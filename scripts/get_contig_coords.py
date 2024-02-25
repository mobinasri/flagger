import sys,re
import argparse


def print_coords(seq_name, seq):
    for match in re.finditer(r'[ATGCatgc]+', seq):
        print(f'{seq_name}\t{match.start() + 1}\t{match.end()}')



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputFasta')
    args = parser.parse_args()

    seq_name = None
    with open(args.inputFasta, "r") as f:
        for line in f:
            if line.strip()[0] == ">":
                if seq_name != None:
                    print_coords(seq_name, seq)
                seq_name = line.strip()[1:]
                seq = ""
                continue
            seq += line.strip()
    print_coords(seq_name, seq)

if __name__ == "__main__":
    main()
