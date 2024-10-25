import sys
import argparse
from collections import defaultdict
import re
from block_utils import *

def main():
    parser = argparse.ArgumentParser(description='Given the alignments find the projection of a set of assembly blocks onto the reference (\'asm2ref\') or vice versa (\'ref2asm\')')
    parser.add_argument('--mode', type=str, default='asm2ref',
                    help='(Str) Default=\'asm2ref\' It can be either {\'ref2asm\' or \'asm2ref\'}. \
                           1.\'asm2ref\': \
                                    In this mode the blocks are given in the coordinates of the assembly and \
                                    the output will be the projections of those blocks onto the reference. \
                           2.\'ref2asm\': \
                                    In this mode the blocks are given in the coordinates of the reference and \
                                    the output will be the projections of those blocks onto the assembly.')
    parser.add_argument('--paf', type=str,
                    help='(PAF format) The alignments of the assembly to the reference. It should include the cigar format.')
    parser.add_argument('--blocks', type=str,
                    help='(BED format) The desired blocks in the assembly (or in the reference if the mode is \'ref2asm\'). It should be sorted and with no overlaps.')
    parser.add_argument('--outputProjectable', type=str,
                    help='(BED format) A path for saving the query blocks that could be projected. The projection of each line is available in the same line of the projection output')
    parser.add_argument('--outputProjection', type=str,
                    help='(BED format) A path for saving the projections of the query blocks.Note that the lines may not be sorted and may have overlaps because of its correspondence with the projected bed file. The 4-th column contains the disimilarity percentage = (indels + mismatches)/(projection block length) * 100. The 5-th column contains other info if present in the input bed file. It is recommended to run bedtools sort (and merge) on this output')
    parser.add_argument('--threads', type=int,
                    help='Number of threads')
    parser.add_argument('--divergence', action='store_true',
                    help='Print divergence percentage (between asm and ref block) as the 4th column in the output bed file')
    parser.add_argument('--flagger', action='store_true',
                    help='Only use when the input bed file in the output of flagger. It will add similar fields to the output bed file.')
    parser.add_argument('--printCigar', action='store_true',
                        help='Add the subset cigar for each projection (in the last column)')
    parser.add_argument('--includeEndingIndel', action='store_true',
                        help='If a projection ended in an indel add that the overlapping indel to the projection')
    parser.add_argument('--includePostIndel', action='store_true',
                        help='If a projection ended right before an indel add that indel will be added to the projection \
                               It will be insertion for ref2asm and deletion for asm2ref \
                             (This option has very limited applications so should not be on usually)')


    # Fetch the arguments
    args = parser.parse_args()
    mode = args.mode
    pafPath = args.paf
    blocksPath = args.blocks
    outputProjectable = args.outputProjectable
    outputProjection = args.outputProjection
    threads = args.threads
    printDiv = args.divergence
    flagger = args.flagger
    printCigar = args.printCigar
    includeEndingIndel = args.includeEndingIndel
    includePostIndel = args.includePostIndel

    # Save the track line if there is one
    trackLine = None
    # Read the desired blocks. Their start and end positions are converted into the 1-based format
    blocks = defaultdict(list)
    with open(blocksPath,"r") as f:
        for line in f:
            if line.startswith("track name"):
                trackLine = line.strip()
                continue
            attrbs = line.strip().split()
            # skip incomplete tracks
            if len(attrbs) < 3:
                continue
            contigName = attrbs[0]
            # start is 0-based in bed format, it gets converted to 1-based here
            start = int(attrbs[1]) + 1
            end = int(attrbs[2])
            info = attrbs[3:] if len(attrbs) > 3 else [""]
            blocks[contigName].append((start, end, info))

    alignments = []
    with open(pafPath,"r") as fPaf:
        for line in fPaf:
            alignments.append(Alignment(line))

    results = runProjectionParallel(alignments, mode, blocks, includeEndingIndel, includePostIndel, threads)

    # Read the alignments one by one and for each of them find the projections by calling findProjections
    with open(outputProjection, "w") as fRef, open(outputProjectable, "w") as fQuery:
        if trackLine != None: # write track line if there was any
                fRef.write(f"{trackLine}\n")
                fQuery.write(f"{trackLine}\n")

        for res in results:
            chromName = res[0]
            contigName = res[1]
            orientation = res[2]
            qBlocks = res[3]
            rBlocks = res[4]
            cigarList = res[5]


            if mode == "asm2ref":
                ctgRef = chromName
                ctgQuery = contigName
            else:
                ctgRef = contigName
                ctgQuery = chromName
            for rBlock, qBlock, cigar in zip(rBlocks, qBlocks, cigarList):
                cigarString = makeCigarString(cigar)
                # Skip if there is no valid projection
                if rBlock[0] == None or qBlock[0] == None: continue
                if flagger:
                    rBlock[2][3] = str(rBlock[0] - 1)
                    rBlock[2][4] = str(rBlock[1])
                if printDiv == True:
                    fRef.write("{}\t{}\t{}\t{:.3f}\t{}".format(ctgRef, rBlock[0] - 1, rBlock[1], rBlock[3], "\t".join(rBlock[2])))
                else:
                    fRef.write("{}\t{}\t{}\t{}".format(ctgRef, rBlock[0] - 1, rBlock[1], "\t".join(rBlock[2])))
                if flagger: 
                    qBlock[2][3] = str(qBlock[0] - 1)
                    qBlock[2][4] = str(qBlock[1])
                if printDiv == True:
                    fQuery.write("{}\t{}\t{}\t{:.3f}\t{}".format(ctgQuery, qBlock[0] - 1, qBlock[1], qBlock[3], "\t".join(qBlock[2])))
                else:
                    fQuery.write("{}\t{}\t{}\t{}".format(ctgQuery, qBlock[0] - 1, qBlock[1], "\t".join(qBlock[2])))

                if printCigar:
                    fQuery.write("\t{}".format(cigarString))
                    fRef.write("\t{}".format(cigarString))
                fQuery.write("\n")
                fRef.write("\n")

if __name__ == "__main__": main()

