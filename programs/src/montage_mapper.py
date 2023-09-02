import sys
import argparse
from collections import defaultdict
from copy import deepcopy
from block_utils import *
import numpy as np
import pandas as pd
import json
import re
import os
import subprocess
import shlex
from concurrent.futures import ThreadPoolExecutor, as_completed
import datetime


def parseAllPafFiles(pafDir):
    pafPaths = [os.path.join(pafDir, f) for f in os.listdir(pafDir) if f.endswith(".paf")]
    alignments = []
    # iterate over all paf files
    for pafPath in pafPaths:
        # parse all alignments from this paf file
        with open(pafPath, "r") as f:
            for line in f:
                alignments.append(Alignment(line))
    return alignments

def runAllAlignments(fastaPathPairs, outDir, minimap2Path, centroalignPath, minimap2Params, centroalignParams,
                     threads=8):
    success = 0

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [
            executor.submit(runMapping, [x, outDir, minimap2Path, centroalignPath, minimap2Params, centroalignParams])
            for x in fastaPathPairs]

        for future in as_completed(futures):
            # retrieve the result
            result = future.result()
            print(result[0])
            if result[1] == 0:
                success += 1
    print(f"[{datetime.datetime.now()}] Successfully mapped [{success}/{len(fastaPathPairs)}] sequence pairs")



def runMapping(args):
    hap1FastaPath = args[0]["hap1_fasta_path"]
    hap2FastaPath = args[0]["hap2_fasta_path"]
    hap1SeqName = args[0]["hap1_seq_name"]
    hap2SeqName = args[0]["hap2_seq_name"]
    mapper = args[0]["mapper"]
    outDir = args[1]
    minimap2Path = args[2]
    centroalignPath = args[3]
    minimap2Params = args[4]
    centroalignParams = args[5]

    outputPafPath = os.path.join(outDir, f"{hap1SeqName}.{hap2SeqName}.{mapper}.paf")

    returncode = 1
    if mapper == "minimap2":
        returncode = runMinimap2(hap1FastaPath,
                                 hap2FastaPath,
                                 outputPafPath,
                                 minimap2Path,
                                 minimap2Params)

    if mapper == "centroalign":
        returncode = runCentroalign(hap1FastaPath,
                                    hap2FastaPath,
                                    outputPafPath,
                                    centroalignPath,
                                    centroalignParams)

    if returncode == 0:
        return f"[{datetime.datetime.now()}] Mapping is done:\t{hap1SeqName}\t{hap2SeqName}\t{mapper}\t{outputPafPath}", 0
    else:
        return f"[{datetime.datetime.now()}] Mapping could not be done:\t{hap2SeqName}\t{hap2SeqName}\t{mapper}", returncode

def runMinimap2(hap1FastaPath, hap2FastaPath, outputPafPath, programPath, params="-x asm5 --eqx"):
    """

    :param params: minimap2 parameters ("-a" is not permitted for producing SAM file)
    :param hap1FastaPath: Path to the fasta file for the first genome
    :param hap2FastaPath: Path to the fasta file for the second genome
    :param outputPafPath: Path to the output Paf file
    """

    cmdString = f"{programPath} -c {params} {hap1FastaPath} {hap2FastaPath} -o {outputPafPath}"
    x = subprocess.run(shlex.split(cmdString), capture_output=True)
    return x.returncode


def runCentroalign(hap1FastaPath, hap2FastaPath, outputPafPath, programPath, params):
    """

    :param params: centroalign parameters
    :param hap1FastaPath: Path to the fasta file for the first genome (should contain only one sequence)
    :param hap2FastaPath: Path to the fasta file for the second genome (should contain only one sequence)
    :param outputPafPath: Path to the output Paf file
    """

    # extract names of the given sequences
    hap1SeqName = None
    hap2SeqName = None
    for name, seq in yieldFasta(hap1FastaPath):
        hap1SeqName = name
        hap1SeqLen = len(seq)
        break

    for name, seq in yieldFasta(hap2FastaPath):
        hap2SeqName = name
        hap2SeqLen = len(seq)
        break

    if hap1SeqName == None or hap2SeqName == None:
        return None

    # concat two sequences in a single fasta file to pass it to centroalign later
    cmdStringConcat = f"cat {hap1FastaPath} {hap2FastaPath}"
    concatFastaPath = os.path.join(os.path.dirname(outputPafPath),
                                   f"{hap1SeqName}_{hap2SeqName}.fasta")
    with open(concatFastaPath, "w") as concatFastaFile:
        x = subprocess.run(shlex.split(cmdStringConcat), stdout=concatFastaFile)
        if x.returncode != 0:
            return x.returncode

    # make a newick file with the two sequence names
    # hap1 name should come first to make it to be used as the reference by centroalign
    treePath = os.path.join(os.path.dirname(outputPafPath),
                            f"{hap1SeqName}_{hap2SeqName}.newick")
    with open(treePath, "w") as f:
        f.write(f"({hap1SeqName}, {hap2SeqName});\n")

    cmdStringAlign = f"{programPath} {params} {concatFastaPath} -T {treePath}"
    x = subprocess.run(shlex.split(cmdStringAlign), capture_output=True)
    if x.returncode != 0:
        return x.returncode

    cigarString = x.stdout.strip().decode('utf-8')
    print(cigarString)
    numberOfMatches, alignmentLength = getNumberOfMatchesAndAlignmentLength(cigarString)
    with open(outputPafPath, "w") as f:
        f.write("\t".join([hap2SeqName, f'{hap2SeqLen}', '0', f'{hap2SeqLen}', '+',
                           hap1SeqName, f'{hap1SeqLen}', '0', f'{hap1SeqLen}',
                           f'{numberOfMatches}', f'{alignmentLength}', '60',
                           'tp:A:P', f'cg:Z:{cigarString}\n']))

    return 0

def shiftAlignmentCoors(alignment, refShift, queryShift):
    """

    :param alignment: Alignment instance parsed from a PAF record
    :param refShift: Amount of shift the start and coordinates of the reference alignment interval should have
    :param queryShift: Amount of shift the start and coordinates of the query alignment interval should have
    """

    alignment.contigStart += queryShift
    alignment.contigEnd += queryShift
    alignment.chromStart += refShift
    alignment.chromEnd += refShift

def getBackToOriginalCoordinates(alignment, contigLengths):
    """
    It is assumed that the ref/query sequences in this alignment are subsets of the original contig/chroms.
    contig/chrom names have the format "{contig_name}_{start}_{end}" illustrating the original coordinates.
    After running this function the names and coordinates of the alignment will be based on the original contig/chrom
    for both ref and query sequences

    :param alignment: Alignment instance
    :param contigLengths: Dictionary of contig lengths
    """

    refSeqName = alignment.chromName
    refOrigName = "_".join(refSeqName.split("_")[0:-2])
    refOrigStart = int(refSeqName.split("_")[-2])

    querySeqName = alignment.contigName
    queryOrigName = "_".join(querySeqName.split("_")[0:-2])
    queryOrigStart = int(querySeqName.split("_")[-2])

    alignment.chromName = refOrigName
    alignment.contigName = queryOrigName

    alignment.chromLength = contigLengths[refOrigName]
    alignment.contigLength = contigLengths[queryOrigName]

    shiftAlignmentCoors(alignment, refOrigStart, queryOrigStart)

def parseIntervalPairsToMap(bedPath):
    """
    :param bedPath: Path to a bed file with 7 columns; hap1_contig, hap1_start, hap1_end, hap2_contig, hap2_start, hap2_end, mapper. mapper can be either "minimap2" or "centroalign"
    :return: List of dictionaries where each dictionary contains the 7 attributes in each bed line
    """

    intervalPairs = []
    with open(bedPath) as f:
        for line in f:
            cols = line.strip().split()
            intervalPairs.append({"hap1_contig": cols[0],
                                  "hap1_start": int(cols[1]),
                                  "hap1_end": int(cols[2]),
                                  "hap2_contig": cols[3],
                                  "hap2_start": int(cols[4]),
                                  "hap2_end": int(cols[5]),
                                  "mapper": cols[6]})
    return intervalPairs

def writeIntervalsToSeparateFastaFiles(intervalPairs, hap1FastaPath, hap2FastaPath, outDir):
    """
    This function creates temporary fasta files for later mappings. Each fasta file will contain the sequence of one
    of the pairs in the given list intervalPairs.

    :param intervalPairs: A list of dictionaries, where each dictionary contains the two seqeuces that has to be mapped against
                          each other
    :param hap1FastaPath: Path to the fasta file for the first genome
    :param hap2FastaPath: Path to the fasta file for the second genome
    :param outDir: output directory
    :return fastaPathPairs: List of dictionaries. Each dictionary has three keys;
                            "hap1_fasta_path", "hap2_fasta_path", "hap1_seq_name", "hap2_seq_name", "mapper"
    """
    intervalsPerContig = defaultdict(list)
    for intervalPair in intervalPairs:
        intervalsPerContig[intervalPair["hap1_contig"]].append(
            (intervalPair["hap1_start"], intervalPair["hap1_end"]))
        intervalsPerContig[intervalPair["hap2_contig"]].append(
            (intervalPair["hap2_start"], intervalPair["hap2_end"]))

    os.makedirs(f"{outDir}", exist_ok=True)

    # for haplotype 1
    for name, seq in yieldFasta(hap1FastaPath):
        for interval in intervalsPerContig[name]:
            s = interval[0]
            e = interval[1]
            newSeq = seq[s:e]
            newSeqName = f"{name}_{s}_{e}"
            with open(f"{outDir}/{newSeqName}.fasta", "w") as f:
                f.write(f">{newSeqName}\n{newSeq}\n")

    # for haplotype 2
    for name, seq in yieldFasta(hap2FastaPath):
        for interval in intervalsPerContig[name]:
            s = interval[0]
            e = interval[1]
            newSeq = seq[s:e]
            newSeqName = f"{name}_{s}_{e}"
            with open(f"{outDir}/{newSeqName}.fasta", "w") as f:
                f.write(f">{newSeqName}\n{newSeq}\n")

    fastaPathPairs = []
    for intervalPair in intervalPairs:
        hap1SeqName = f'{intervalPair["hap1_contig"]}_{intervalPair["hap1_start"]}_{intervalPair["hap1_end"]}'
        hap2SeqName = f'{intervalPair["hap2_contig"]}_{intervalPair["hap2_start"]}_{intervalPair["hap2_end"]}'
        hap1FastaPath = os.path.join(outDir, f"{hap1SeqName}.fasta")
        hap2FastaPath = os.path.join(outDir, f"{hap2SeqName}.fasta")
        fastaPathPairs.append({"hap1_fasta_path": hap1FastaPath,
                               "hap2_fasta_path": hap2FastaPath,
                               "hap1_seq_name": hap1SeqName,
                               "hap2_seq_name": hap2SeqName,
                               "mapper": intervalPair["mapper"]})
    return fastaPathPairs

def main():
    parser = argparse.ArgumentParser(
        description='This program enables mapping different parts of two genomes with different mappers and outputs a paf file')
    parser.add_argument('--bed', type=str,
                        help='It is a bed file with 7 columns; hap1_contig, hap1_start, hap1_end, hap2_contig, hap2_start, hap2_end, mapper. mapper can be either "minimap2" or "centroalign"')
    parser.add_argument('--hap1', type=str,
                        help='hap1 fasta file')
    parser.add_argument('--hap2', type=str,
                        help='hap2 fasta file')
    parser.add_argument('--outDir', type=str,
                        help='output directory')
    parser.add_argument('--prefix', type=str, default="montage_mapper",
                        help='output prefix')
    parser.add_argument('--minimap2Path', type=str,
                        help='path to minimap2 executable')
    parser.add_argument('--centroalignPath', type=str,
                        help='path to centroalign executable')
    parser.add_argument('--minimap2Params', type=str, default="-x asm5 --eqx",
                        help='minimap2 parameters [Default "-x asm5 --eqx"]')
    parser.add_argument('--centroalignParams', type=str, default="",
                        help='centroalign parameters [Default ""]')

    args = parser.parse_args()
    bedPath = args.bed
    hap1FastaPath = args.hap1
    hap2FastaPath = args.hap2
    outDir = args.outDir
    prefix = args.prefix
    minimap2Path = args.minimap2Path
    centroalignPath = args.centroalignPath
    minimap2Params = args.minimap2Params
    centroalignParams = args.centroalignParams

    # save contig lengths
    contigLengths = {}
    for name, seq in yieldFasta(hap1FastaPath):
        contigLengths[name] = len(seq)
    for name, seq in yieldFasta(hap2FastaPath):
        contigLengths[name] = len(seq)

    print(f"[{datetime.datetime.now()}] Parsing all requested intervals that should be mapped against each other.")
    intervalPairs = parseIntervalPairsToMap(bedPath)
    print(f"[{datetime.datetime.now()}] {len(intervalPairs)} interval pairs are parsed.")

    print(f'[{datetime.datetime.now()}] Writing the sequence of each interval in a separate fasta file in {os.path.join(outDir, "fasta_tmp")}')
    fastaPathPairs = writeIntervalsToSeparateFastaFiles(intervalPairs, hap1FastaPath, hap2FastaPath,
                                                        os.path.join(outDir, "fasta_tmp"))


    pafDir = os.path.join(outDir, "paf_temp")
    os.makedirs(pafDir, exist_ok=True)

    print(f"[{datetime.datetime.now()}] Running all requested alignments and save each one in a separate paf file in {pafDir}")
    runAllAlignments(fastaPathPairs,
                     pafDir,
                     minimap2Path, centroalignPath,
                     minimap2Params, centroalignParams,
                     threads=8)

    print(f"[{datetime.datetime.now()}] Parsing all alignments from the paf files in {pafDir}")
    alignments = parseAllPafFiles(pafDir)


    finalPafPath = os.path.join(outDir, f"{prefix}.paf")
    print(f"[{datetime.datetime.now()}] Adjusting the coordiantes of the alignments and saving them in {finalPafPath}")
    for alignment in alignments:
        # adjust the coordinates
        getBackToOriginalCoordinates(alignment, contigLengths)
        # write the new alignment record to the final paf file
        # append should be true to avoid writing from the beginning of the file
        alignment.writeToPaf(finalPafPath, append=True)

    print(f"[{datetime.datetime.now()}] {len(alignments)} alignments are written to {finalPafPath}.")
    print(f"[{datetime.datetime.now()}] Finished!")


if __name__ == "__main__": main()
