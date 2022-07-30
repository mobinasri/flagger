import sys
import argparse
from collections import defaultdict
import re
from block_utils import findProjections, Alignment


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
                    help='(BED format) A path for saving the projections of the query blocks.Note that the lines may not be sorted and may have overlaps because of its correspondence with the projected bed file. It is recommended to run bedtools sort (and merge) on this output')
    
    # Fetch the arguments
    args = parser.parse_args()
    mode = args.mode
    pafPath = args.paf
    blocksPath = args.blocks
    outputProjectable = args.outputProjectable
    outputProjection = args.outputProjection

    # Read the desired blocks. Their start and end positions are converted into the 1-based format
    blocks = defaultdict(list)
    with open(blocksPath,"r") as f:
        for line in f:
            attrbs = line.strip().split()
            contigName = attrbs[0]
            # start is 0-based in bed format, it gets converted to 1-based here
            start = int(attrbs[1]) + 1
            end = int(attrbs[2])
            blocks[contigName].append((start, end))

    # Read the alignments one by one and for each of them find the projections by calling findProjections
    with open(pafPath,"r") as fPaf, open(outputProjection, "w") as fRef, open(outputProjectable, "w") as fQuery:
        for line in fPaf:
            # Extract the alignment attributes like the contig name, alignment boundaries, orientation and cigar 
            alignment = Alignment(line);
            chromName = alignment.chromName
            contigName = alignment.contigName
            # rBlocks contains the projections and 
            # qBlocks contains the projectable blocks
            if mode == "asm2ref":
                if len(blocks[contigName]) == 0: # Continue if there is no block in the contig
                    continue
                #print(blocks[contigName], contigStart, contigEnd, chrom, chromStart, chromEnd)
                qBlocks, rBlocks = findProjections(mode, 
                                                   alignment.cigarList, 
                                                   blocks[contigName],
                                                   alignment.chromLength, 
                                                   alignment.chromStart + 1, alignment.chromEnd, # make 1-based start
                                                   alignment.contigLength, 
                                                   alignment.contigStart + 1, alignment.contigEnd, # make 1-based start
                                                   alignment.orientation)
            else:
                if len(blocks[chromName]) == 0: # Continue if there is no block in the chrom
                    continue
                qBlocks, rBlocks = findProjections(mode, 
                                                   alignment.cigarList, 
                                                   blocks[chromName],
                                                   alignment.chromLength, 
                                                   alignment.chromStart + 1, alignment.chromEnd, # make 1-based start
                                                   alignment.contigLength, 
                                                   alignment.contigStart + 1, alignment.contigEnd, # make 1-based start
                                                   alignment.orientation)
            if mode == "asm2ref":
                for rBlock in rBlocks:
                    fRef.write("{}\t{}\t{}\n".format(chromName, rBlock[0] - 1, rBlock[1]))
                for qBlock in qBlocks:
                    fQuery.write("{}\t{}\t{}\n".format(contigName, qBlock[0] - 1, qBlock[1]))
            else: # mode = "ref2asm"
                for rBlock in rBlocks:
                    fRef.write("{}\t{}\t{}\n".format(contigName, rBlock[0] - 1, rBlock[1]))
                for qBlock in qBlocks:
                    fQuery.write("{}\t{}\t{}\n".format(chromName, qBlock[0] - 1, qBlock[1]))
main()

