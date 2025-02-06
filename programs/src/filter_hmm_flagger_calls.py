import sys
import argparse
from collections import defaultdict
import re
from block_utils import *
import pysam
import copy
import os
import datetime




def parseBlocksFromBed(inputBed):
    # Read the desired blocks. Their start and end positions are converted into the 1-based format
    inputBlocks = defaultdict(list)
    with open(inputBed,"r") as f:
        for line in f:
            if line.startswith("track name"):
                continue
            attrbs = line.strip().split()
            chromName = attrbs[0]
            # start is 0-based in bed format, it gets converted to 1-based here
            start = int(attrbs[1]) + 1
            end = int(attrbs[2])
            info = {"ref_label": attrbs[3], "query_labels":{"asm2ref": [], "ref2asm": []}}
            inputBlocks[chromName].append((start, end, info))
    return inputBlocks

def parseAlignmentsFromSam(inputSamPath, minAlignmentLen, maxDivergence):
    alignments = []
    inputPysamAlignmentFile = pysam.AlignmentFile(inputSamPath, "r")
    for inputPysamRecord in inputPysamAlignmentFile.fetch():
        alignment = Alignment.createFromPysamRecord(inputPysamRecord, inputPysamAlignmentFile.header)
        if (alignment.alignmentLength > minAlignmentLen and
                alignment.divergence is not None and
                alignment.divergence < maxDivergence):
            alignments.append(alignment)
    return  alignments


def saveConservativeBlocksInBed(blockListPerContig, bedPath):
    label2color = {"Err": "162,0,37",
                   "Dup": "250,104,0",
                   "Hap": "0,138,0",
                   "Col": "170,0,255",
                   "Unk": "99, 99, 96",
                   "Msj": "250,200,0",
                   "NNN": "0,0,0"}
    with open(bedPath, "w") as outFile:
        for ctg, blockList in blockListPerContig.items():
            preLabel = None
            preStart = None
            preEnd = None
            for start, end, label in blockList.blocks:
                if preLabel == None:
                    # the first block of this contig
                    # set label and coors
                    preLabel = label
                    preStart = start
                    preEnd = end
                elif preLabel == label:
                    # update end only
                    preEnd = end
                else: # label has changed so write the previous block and label
                    outFile.write("{}\t{}\t{}\t{}\t0\t.\t{}\t{}\t{}\n".format(ctg,
                                                          preStart - 1,
                                                          preEnd,
                                                          preLabel,
                                                          preStart - 1,
                                                          preEnd,
                                                          label2color[preLabel]))
                    # update label and coors
                    preLabel = label
                    preStart = start
                    preEnd = end

            # write the last block of the contig
            if preLabel != None:
                outFile.write("{}\t{}\t{}\t{}\t0\t.\t{}\t{}\t{}\n".format(ctg,
                                                      preStart - 1,
                                                      preEnd,
                                                      preLabel,
                                                      preStart - 1,
                                                      preEnd,
                                                      label2color[preLabel]))

def main():
    parser = argparse.ArgumentParser(description='Filter regions flagged by hmm-flagger using the self homology mappings of the assembly contigs.')
    parser.add_argument('--inputBam', type=str,
                        help='(BAM/SAM format) The self homology mappings of the assembly contigs. This bam file can be created by running minimap2 with the parameters "-D -ax asm5" by using the diploid assembly both as reference and query.')
    parser.add_argument('--inputBed', type=str,
                        help='(BED format) The hmm-flagger output bed file.')
    parser.add_argument('--outputBed', type=str,
                        help='(BAM/SAM format) The filtered hmm-flagger bed file with lower number of false positive calls.')
    parser.add_argument('--threads', type=int, default=4,
                        help='Number of threads for performing label projections between haplotypes.')
    parser.add_argument('--maxDivergence', type=float, default=0.005,
                        help='The alignment records with gap-compressed divergence rate higher than this will be skipped. (Default=0.01)')
    parser.add_argument('--minAlignmentLen', type=int, default=10000,
                        help='The alignment records shorter than this will be skipped. (Default=10000)')


    # Fetch the arguments
    args = parser.parse_args()
    inputSamPath = args.inputBam
    inputBed = args.inputBed
    outputBed = args.outputBed
    threads = args.threads
    maxDivergence = args.maxDivergence
    minAlignmentLen = args.minAlignmentLen

    inputBlocks = parseBlocksFromBed(inputBed)
    print(f"[{datetime.datetime.now()}] Input HMM-Flagger blocks have been parsed.")

    print(f"[{datetime.datetime.now()}] Started parsing alignments.")

    alignments = parseAlignmentsFromSam(inputSamPath, minAlignmentLen, maxDivergence)

    print(f"[{datetime.datetime.now()}] Started running label projections using {threads} threads.")
    projections_asm2ref = runProjectionParallel(alignments,
                                    'asm2ref',
                                    blocks = inputBlocks,
                                    includeEndingIndel = False,
                                    includePostIndel = False,
                                    threads = threads)
    projections_ref2asm = runProjectionParallel(alignments,
                                    'ref2asm',
                                    blocks = inputBlocks,
                                    includeEndingIndel = False,
                                    includePostIndel = False,
                                    threads = threads)


    augmentedBlocks = copy.deepcopy(inputBlocks)

    ############
    # asm to ref
    ############
    for proj in projections_asm2ref:
        ctgRef = proj[0]
        ctgQuery = proj[1]
        orientation = proj[2]
        qBlocks = proj[3]
        rBlocks = proj[4]

        for rBlock, qBlock in zip(rBlocks, qBlocks):
            # for ref block
            rInfo = copy.deepcopy(rBlock[2])
            rInfo["ref_label"] = None
            # the ref label of the query block will be added to the query labels
            rInfo["query_labels"] = {"asm2ref": [qBlock[2]["ref_label"]], "ref2asm": []}
            augmentedBlocks[ctgRef].append((rBlock[0], rBlock[1], rInfo))


    ############
    # ref to asm
    ############
    for proj in projections_ref2asm:
        ctgRef = proj[1]
        ctgQuery = proj[0]
        orientation = proj[2]
        qBlocks = proj[3]
        rBlocks = proj[4]

        for rBlock, qBlock in zip(rBlocks, qBlocks):
            # for ref block
            rInfo = copy.deepcopy(rBlock[2])
            rInfo["ref_label"] = None
            # the ref label of the query block will be added to the query labels
            rInfo["query_labels"] = {"asm2ref": [], "ref2asm": [qBlock[2]["ref_label"]]}
            augmentedBlocks[ctgRef].append((rBlock[0], rBlock[1], rInfo))

    def mergeInfo(info1, info2):
        queryLabels_ref2asm = []
        queryLabels_asm2ref = []
        
        queryLabels_ref2asm.extend(info1["query_labels"]["ref2asm"])
        queryLabels_ref2asm.extend(info2["query_labels"]["ref2asm"])

        queryLabels_asm2ref.extend(info1["query_labels"]["asm2ref"])
        queryLabels_asm2ref.extend(info2["query_labels"]["asm2ref"])

        refLabel = info2["ref_label"] if info1["ref_label"] == None else info1["ref_label"]
        if  info1["ref_label"] != None and info2["ref_label"] != None and info1["ref_label"] != info2["ref_label"]:
            print("ref labels do not match:", info1["ref_label"], info2["ref_label"])
            sys.exit(1)
        return {"ref_label": refLabel, "query_labels": {"ref2asm": queryLabels_ref2asm,  "asm2ref": queryLabels_asm2ref}}



    conservativeBlockListPerContig = {}
    with open(outputBed.replace(".bed", ".debug.bed"), "w") as outFile:
        # wrap block lists with BlockList
        for ctg, blocks in augmentedBlocks.items():
            conservativeBlockListPerContig[ctg] = BlockList()
            blockList = BlockList(blocks)
            blockList.mergeWithCustomFunction(mergeInfo, inplace=True)
            for start, end, info in blockList.blocks:
                ref_label = info["ref_label"]
                query_labels_ref2asm = info["query_labels"]["ref2asm"]
                query_labels_asm2ref = info["query_labels"]["asm2ref"]
                conservative_label = None
                if ref_label == "Err" or ref_label == "Hap" or ref_label == "Unk" or ref_label == "NNN":
                    conservative_label = ref_label
                if ref_label == "Col" or ref_label == "Dup":
                    # change to Hap if there are too many contig mappings (>=3)
                    if len(query_labels_asm2ref) >= 3 or len(query_labels_ref2asm) >= 3:
                        conservative_label = "Hap"
                    else:
                        # if ref label is Dup and there is at least one contig with Dup flag mapped to this block
                        if ref_label == "Dup" and ("Dup" in query_labels_asm2ref or "Dup" in query_labels_ref2asm):
                            conservative_label = "Dup"
                        # if ref label is Col and there is either no query labels or at least one contig with Err flag mapped to this block
                        elif ref_label == "Col" and (("Err" in query_labels_asm2ref) or ("Err" in query_labels_ref2asm) or (len(query_labels_asm2ref) == 0) or (len(query_labels_ref2asm) == 0)):
                            conservative_label = "Col"
                        else: # report Hap in all other cases
                            conservative_label = "Hap"
                outFile.write(f"{ctg}\t{start}\t{end}\t({','.join(query_labels_ref2asm)}),({','.join(query_labels_asm2ref)})\n")
                conservativeBlockListPerContig[ctg].append((start, end, conservative_label))
    # save blocks with conservative labels in the output bed file
    saveConservativeBlocksInBed(conservativeBlockListPerContig, outputBed)

    print(f"[{datetime.datetime.now()}] Done! Output BED file : ", outputBed)

if __name__ == "__main__": main()

