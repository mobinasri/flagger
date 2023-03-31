import sys
import argparse
from collections import defaultdict
import re
from block_utils import findProjections, Alignment
from multiprocessing import Pool

def runProjection(line, mode, blocks):
    # Extract the alignment attributes like the contig name, alignment boundaries, orientation and cigar 
    alignment = Alignment(line);
    chromName = alignment.chromName
    contigName = alignment.contigName
    if alignment.isPrimary == False:
        return [chromName, contigName, [], []]
    # rBlocks contains the projections and 
    # qBlocks contains the projectable blocks
    if mode == "asm2ref":
        if len(blocks[contigName]) == 0: # Continue if there is no block in the contig
            return [chromName, contigName, [], []]
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
            return [chromName, contigName, [], []]
        qBlocks, rBlocks = findProjections(mode,
                                            alignment.cigarList,
                                            blocks[chromName],
                                            alignment.chromLength,
                                            alignment.chromStart + 1, alignment.chromEnd, # make 1-based start
                                            alignment.contigLength,
                                            alignment.contigStart + 1, alignment.contigEnd, # make 1-based start
                                            alignment.orientation)

    return [chromName, contigName, qBlocks, rBlocks]

def runProjectionParallel(pafPath, mode, blocks, threads):
    allPafLines = []
    with open(pafPath,"r") as fPaf:
        for line in fPaf:
            allPafLines.append(line)
    pool = Pool(threads)
    print("Started projecting")
    results = pool.starmap(runProjection, [(line, mode, blocks) for line in allPafLines])
    pool.close()
    return results



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
    parser.add_argument('--divergence',
                    help='Print divergence percentage (between asm and ref block) as the 4th column in the output bed file')
    
    # Fetch the arguments
    args = parser.parse_args()
    mode = args.mode
    pafPath = args.paf
    blocksPath = args.blocks
    outputProjectable = args.outputProjectable
    outputProjection = args.outputProjection
    threads = args.threads
    printDiv = args.divergence

    # Save the track line if there is one
    trackLine = None
    # Read the desired blocks. Their start and end positions are converted into the 1-based format
    blocks = defaultdict(list)
    with open(blocksPath,"r") as f:
        for line in f:
            if line.startswith("track name"):
                trackLine = line
                continue
            attrbs = line.strip().split()
            contigName = attrbs[0]
            # start is 0-based in bed format, it gets converted to 1-based here
            start = int(attrbs[1]) + 1
            end = int(attrbs[2])
            info = attrbs[3:] if len(attrbs) > 3 else [""]
            blocks[contigName].append((start, end, info))

    results = runProjectionParallel(pafPath, mode, blocks, threads)

    # Read the alignments one by one and for each of them find the projections by calling findProjections
    with open(outputProjection, "w") as fRef, open(outputProjectable, "w") as fQuery:
        for res in results:
            chromName = res[0]
            contigName = res[1]
            qBlocks = res[2]
            rBlocks = res[3]
            if trackLine != None: # write track line if there was any
                fRef.write(f"{trackLine}\n")
                fQuery.write(f"{trackLine}\n")

            if mode == "asm2ref":
                if printDiv == True: # print divergence percentage in the 4th column
                    for rBlock in rBlocks:
                        fRef.write("{}\t{}\t{}\t{:.3f}\t{}\n".format(chromName, rBlock[0] - 1, rBlock[1], rBlock[3], rBlock[2]))
                    for qBlock in qBlocks:
                        fQuery.write("{}\t{}\t{}\t{:.3f}\t{}\n".format(contigName, qBlock[0] - 1, qBlock[1], qBlock[3], qBlock[2]))
                else:  # skip divergence ratio
                    for rBlock in rBlocks:
                        fRef.write("{}\t{}\t{}\t{}\n".format(chromName, rBlock[0] - 1, rBlock[1], "\t".join(rBlock[2])))
                    for qBlock in qBlocks:
                        fQuery.write("{}\t{}\t{}\t{}\n".format(contigName, qBlock[0] - 1, qBlock[1], "\t".join(qBlock[2])))
            else: # mode = "ref2asm"
                if printDiv == True: # print divergence percentage in the 4th column
                    for rBlock in rBlocks:
                        fRef.write("{}\t{}\t{}\t{:.3f}\t{}\n".format(contigName, rBlock[0] - 1, rBlock[1], rBlock[3], rBlock[2]))
                    for qBlock in qBlocks:
                        fQuery.write("{}\t{}\t{}\t{:.3f}\t{}\n".format(chromName, qBlock[0] - 1, qBlock[1], rBlock[3], qBlock[2]))
                else: # skip divergence ratio
                    for rBlock in rBlocks:
                        fRef.write("{}\t{}\t{}\t{}\n".format(contigName, rBlock[0] - 1, rBlock[1], "\t".join(rBlock[2])))
                    for qBlock in qBlocks:
                        fQuery.write("{}\t{}\t{}\t{}\n".format(chromName, qBlock[0] - 1, qBlock[1], "\t".join(qBlock[2])))
main()

