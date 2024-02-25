import sys,re
import argparse


def break_into_contigs(seq_name, seq):
    ctg_seqs = re.sub('[nN]+','\n',seq).split('\n')
    if len(ctg_seqs) == 1:
        print(f'>{seq_name}')
        print(ctg_seqs[0])
    else:
        for j, ctg_seq in enumerate(ctg_seqs):
            print(f'>{seq_name}_{j}')
            print(ctg_seq)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputFasta')
    args = parser.parse_args()

    seq_name = None
    with open(args.inputFasta, "r") as f:
        for line in f:
            if line.strip()[0] == ">":
                if seq_name != None:
                    break_into_contigs(seq_name, seq)
                seq_name = line.strip()[1:]
                seq = ""
                continue
            seq += line.strip()
    break_into_contigs(seq_name, seq)

if __name__ == "__main__":
    main()
