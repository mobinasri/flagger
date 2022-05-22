import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description='Split a BED file into multiple ones')
    parser.add_argument('--bed', type=str, help='Input BED file')
    parser.add_argument('--n', type=str, help='Number of output bed files')
    parser.add_argument('--dir', type=str, help='Directory for saving split bed files')
    parser.add_argument('--prefix', type=str, help='Prefix of the output bed files')
    args = parser.parse_args()
    bedPath = args.bed
    prefix = args.prefix
    n = int(args.n)
    outputDir = args.dir
    BedBlocks = []
    totalSize = 0
    with open(bedPath,"r") as f:
        for line in f:
            elems = line.strip().split()
            contig = elems[0]
            start = int(elems[1])
            end = int(elems[2])
            BedBlocks.append((contig,start,end))
            totalSize += end - start
    intervalSize = int(totalSize / n) + 1
    currBlock = BedBlocks[0]
    currInterval = 0
    startPoint = (currBlock[0],0)
    blockIdx = 0
    for i in range(1,n+1):
        currInterval = 0
        with open("{}/{}_{}.bed".format(outputDir,prefix,i),"w") as f:
            while (startPoint[1] + intervalSize - currInterval) >= currBlock[2]:
                f.write("{}\t{}\t{}\n".format(startPoint[0],startPoint[1],currBlock[2]))
                currInterval += currBlock[2] - startPoint[1]
                blockIdx += 1
                if blockIdx >= len(BedBlocks):
                    break
                currBlock = BedBlocks[blockIdx]
                startPoint = (currBlock[0],currBlock[1])
            if (blockIdx >= len(BedBlocks)):
                    break
            if (currInterval < intervalSize):
                f.write("{}\t{}\t{}\n".format(startPoint[0],startPoint[1],startPoint[1] + intervalSize - currInterval))
                startPoint = (startPoint[0], startPoint[1] + intervalSize - currInterval)

main()
