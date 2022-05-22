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
    nContigs = 0
    with open(bedPath,"r") as f:
        for line in f:
            nContigs += 1
            elems = line.strip().split()
            contig = elems[0]
            start = int(elems[1])
            end = int(elems[2])
            BedBlocks.append((contig,start,end))
            totalSize += end - start
    n = min(n, nContigs)
    BedBlocks.sort(key= lambda x : x[1] - x[2])
    intervalSize = int(totalSize / (n+1)) + 1
    groupContigs = [[BedBlocks[i]] for i in range(n)] # initialize groups with the longest contigs
    groupSizes = [BedBlocks[i][2] - BedBlocks[i][1] for i in range(n)]
    contigIdx = n
    print("Interval: {} Mb".format(intervalSize/1e6))
    for i in range(n):
        while groupSizes[i] < intervalSize and contigIdx < len(BedBlocks):
            groupContigs[i].append(BedBlocks[contigIdx])
            groupSizes[i] += BedBlocks[contigIdx][2] - BedBlocks[contigIdx][1]
            contigIdx += 1
    # Distribute the remaining small contigs into different groups
    for i in range(contigIdx, len(BedBlocks)):
            groupContigs[i%n].append(BedBlocks[contigIdx])
            groupSizes[i%n] += BedBlocks[contigIdx][2] - BedBlocks[contigIdx][1]
            contigIdx += 1

    total = 0 
    for i in range(n):
        print("{}:\t{:.2f} Mb".format(i, groupSizes[i]/1e6))
        with open("{}/{}_{}.bed".format(outputDir, prefix, i+1),"w") as f:
            for contig, start, end in groupContigs[i]:
                total += end - start
                f.write("{}\t{}\t{}\n".format(contig, start, end))
    # check if the total size is correct
    assert(total == totalSize)
main()
