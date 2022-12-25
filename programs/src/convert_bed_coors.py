import sys

# If the assembly contigs are split and their contig names have this format contigname:start-end
# this short script can convert the coordinates back to the original coordinates

def main():
    with open(sys.argv[1]) as f:
        for line in f:
            cols = line.strip().split()
            start = int(cols[1])
            end = int(cols[2])
            window = cols[0]
            ctg, interval = window.split(':')
            offset = int(interval.split("-")[0])
            print(f'{ctg}\t{offset + start}\t{offset + end}')

main()
