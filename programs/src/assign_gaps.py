from collections import defaultdict
import sys
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description='Find and annotate gaps (gaps are the blocks that are not present in any of the given bed files). It is assumed that there is no common block between any two given bed files. Each given bed file will produce an output bed file containing the assigned gaps.')
    parser.add_argument('--fai', type=str, help='fai file')
    parser.add_argument('files', nargs='+')
    args = parser.parse_args()
    bedFilepaths = args.files
    faiFilepath = args.fai

    # Get the lengths of all contigs from fai file
    contigLengths = {}
    with open(faiFilepath, "r") as f:
            for line in f:
                attrbs = line.strip().split()
                contig = attrbs[0]
                cLength = int(attrbs[1])
                contigLengths[contig] = cLength
    
    # Get the blocks from the given bed files
    blocks = []
    for i in range(len(bedFilepaths)):
        filepath = bedFilepaths[i]
        with open(filepath, "r") as f:
            for line in f:
                attrbs = line.strip().split()
                contig = attrbs[0]
                start = int(attrbs[1])
                end = int(attrbs[2])
                blocks.append((contig, start, end, i)) # "i" is the index of the bed file

    # Each given bed file will produce an output bed file containing the assigned gaps
    outFiles = []
    for i in range(len(bedFilepaths)):
        filepath = bedFilepaths[i]
        outFilename = os.path.splitext(os.path.basename(filepath))[0]
        outFiles.append(open("{}.gaps.bed".format(outFilename),"w"))

    blocks.sort()
    for i in range(1, len(blocks)):
        prevBlock = blocks[i-1]
        currBlock = blocks[i]
        if currBlock[0] != prevBlock[0]: # means that a new contig starts here, so currBlock is the first block of the new contig and prevBlock is the last block of the previous contig
            if prevBlock[2] < contigLengths[prevBlock[0]]: # if the last block does not cover the last bases of the contig
                outFiles[currBlock[3]].write("{}\t{}\t{}\n".format(prevBlock[0], prevBlock[2], contigLengths[prevBlock[0]]))
            if currBlock[1] > 0: # if the first block does not cover the first bases of the contig
                outFiles[currBlock[3]].write("{}\t{}\t{}\n".format(currBlock[0], 0, currBlock[1]))
        else: # if we are still in the previous contig
            gapStart = prevBlock[2]
            gapEnd = currBlock[1]
            gapLength = gapEnd - gapStart
            if currBlock[3] == prevBlock[3]: # if the type of block is same as the previous one
                outFiles[currBlock[3]].write("{}\t{}\t{}\n".format(currBlock[0], gapStart, gapEnd))
            elif gapLength > 0: # if the type of block is different from the previous one and there is a gap
                prevBlockLength = prevBlock[2] - prevBlock[1]
                currBlockLength = currBlock[2] - currBlock[1]
                # Here I assume that the larger block takes larger portion of the unannotated gap (linearly)
                prevW = prevBlockLength / (prevBlockLength + currBlockLength)
                gapMidpoint = int(gapLength * prevW) + gapStart
                outFiles[prevBlock[3]].write("{}\t{}\t{}\n".format(currBlock[0], gapStart, gapMidpoint))
                outFiles[currBlock[3]].write("{}\t{}\t{}\n".format(currBlock[0], gapMidpoint, gapEnd))
    for outFile in outFiles:
        outFile.flush()
        outFile.close()



if __name__ == "__main__":
    main()

