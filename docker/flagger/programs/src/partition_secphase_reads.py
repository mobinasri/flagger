import sys
import argparse
import math


def main():
    parser = argparse.ArgumentParser(description='Given the secphase output, generate two lists of reads; one list for the reads that should be transferred from hap1 to hap2 and the other one for the hap2-to-hap1 reads')
    parser.add_argument('--secphase', type=str,
                    help='the output of secphase')
    parser.add_argument('--output1', type=str,
                    help='File path to output hap1-to-hap2 read names')
    parser.add_argument('--output2', type=str,
                    help='File path to output hap2-to-hap1 read names')
    args = parser.parse_args()
    secphaseOutputFile = args.secphase
    outputFile1 = args.output1
    outputFile2 = args.output2

    with open(secphaseOutputFile, "r") as fin, open(outputFile1, "w") as fout1, open(outputFile2, "w") as fout2:
        readname = None
        prevContig = None
        nextContig = None
        for line in fin:
            if line.strip() == "":
                continue
            cols = line.strip().split()

            if cols[0] == "$" and readname != None and prevContig != None and nextContig != None:
                #print(readname, prevContig, nextContig)
                # Check if it is hap1-to-hap2 or the reverse
                if "#1" in prevContig and "#2" in nextContig:
                    fout1.write("{}\n".format(readname))
                elif "#2" in prevContig and "#1" in nextContig:
                    fout2.write("{}\n".format(readname))

            if cols[0] == "$":
                if "RECORDS" in line:
                    readname = cols[-1][:-1] # remove ":"
                else:
                    readname = cols[-1]
            elif cols[0] == "*":
                prevContig = cols[2]
            elif cols[0] == "@":
                nextContig = cols[2]
        if readname != None:
            # Check if it is hap1-to-hap2 or the reverse
            if "#1" in prevContig and "#2" in nextContig:
                fout1.write("{}\n".format(readname))
            elif "#2" in prevContig and "#1" in nextContig:
                fout2.write("{}\n".format(readname))


if __name__ == "__main__":
    main()
