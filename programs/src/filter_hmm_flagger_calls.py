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
            info = {"ref_label": attrbs[3], "query_labels":[]}
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
    parser.add_argument('--maxDivergence', type=float, default=0.01,
                        help='The alignment records with gap-compressed divergence rate higher than this will be skipped. (Default=0.01)')
    parser.add_argument('--minAlignmentLen', type=int, default=10000,
                        help='The alignment records shorter than this will be skipped. (Default=10000)')


    # Fetch the arguments
    args = parser.parse_args()
    inputSamPath = args.inputBam
    inputBed = args.inputBed
    outputBed = args.outputBed
    threads = args.threads
    maxDivergence = parser.maxDivergence
    minAlignmentLen = parser.minAlignmentLen

    inputBlocks = parseBlocksFromBed(inputBed)
    print(f"[{datetime.datetime.now()}] Input HMM-Flagger blocks have been parsed.")

    print(f"[{datetime.datetime.now()}] Started parsing alignments.")

    alignments = parseAlignmentsFromSam(inputSamPath, minAlignmentLen, maxDivergence)

    print(f"[{datetime.datetime.now()}] Started running label projections using {threads} threads.")
    projections = runProjectionParallel(alignments,
                                    'asm2ref',
                                    blocks = inputBlocks,
                                    includeEndingIndel = False,
                                    includePostIndel = False,
                                    threads = threads)


    augmentedBlocks = copy.deepcopy(inputBlocks)
    for proj in projections:
        ctgRef = proj[0]
        ctgQuery = proj[1]
        orientation = proj[2]
        qBlocks = proj[3]
        rBlocks = proj[4]

        for rBlock, qBlock in zip(rBlocks, qBlocks):
            # for ref block
            rInfo = copy.deepcopy(rBlock[2])
            # the ref label of the query block will be added to the query labels
            rInfo["query_labels"] = [qBlock[2]["ref_label"]]
            augmentedBlocks[ctgRef].append((rBlock[0],rBlock[1],rInfo))

            # for query block
            qInfo = copy.deepcopy(qBlock[2])
            # the ref label of the ref block will be added to the query labels
            qInfo["query_labels"] = [rBlock[2]["ref_label"]]
            augmentedBlocks[ctgQuery].append((qBlock[0],qBlock[1],qInfo))


    def mergeInfo(info1, info2):
        queryLabels = []
        queryLabels.extend(info1["query_labels"])
        queryLabels.extend(info2["query_labels"])

        if info1["ref_label"] != info2["ref_label"]:
            print("ref labels do not match:", info1["ref_label"], info2["ref_label"])
            sys.exit(1)
        return {"ref_label": info1["ref_label"], "query_labels": queryLabels}

    with open(outputBed, "w") as outFile:
        # wrap block lists with BlockList
        for ctg, blocks in augmentedBlocks.items():
            blockList = BlockList(blocks)
            blockList.mergeWithCustomFunction(mergeInfo, inplace=True)
            for start, end, info in blockList.blocks:
                ref_label = info["ref_label"]
                query_labels = info["query_labels"]
                outFile.write("{}\t{}\t{}\t{}".format(ctg,
                                                      start - 1,
                                                      end,
                                                      f'ref={ref_label}|query=' + ",".join(query_labels)))
    print(f"[{datetime.datetime.now()}] Done! Output BED file : ", outputBed)

if __name__ == "__main__": main()

