import argparse

# This script takes a fai file and prints out a bed file
# covering the whole assembly with each track not being 
# longer than 2 times the given split size

def main():
    parser = argparse.ArgumentParser(description='# This script takes a fai file and prints out a bed file covering the whole assembly with each track not being longer than 2 times the given split size. It can be useful for splitting an assembly into shorter contigs')
    parser.add_argument('--fai', type=str, help='path to fasta index file')
    parser.add_argument('--splitSize', type=int, help='split size', default=None)
    args = parser.parse_args()
    faiPath = args.fai
    splitSize = args.splitSize

    ctg2len = {}
    with open(faiPath, "r") as f:
        for line in f:
            cols = line.split()
            ctg = cols[0]
            length = int(cols[1])
            ctg2len[ctg] = length
    for ctg, length in ctg2len.items():
        n = int(length / splitSize)
        if 1 < n:
            for i in range(1,n):
                s = (i-1) * splitSize
                e = i * splitSize
                print(f'{ctg}\t{s}\t{e}')
            s = (n-1) * splitSize
            e = length
            print(f'{ctg}\t{s}\t{e}')
        else:
            s = 0
            e = length
            print(f'{ctg}\t{s}\t{e}')

main()
