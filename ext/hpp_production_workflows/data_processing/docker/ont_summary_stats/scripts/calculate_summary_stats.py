#/usr/bin/env python
# Miten Jain (mjain3@ucsc.edu)
# calculate_summary_stats.py
# updated to use ONT files with multiple headers

import os, sys, time
import numpy
import gzip
from optparse import OptionParser

########################################################################
# Main
# Here is the main program
########################################################################

def main(myCommandLine=None):

    t0 = time.time()

    #Parse the inputs args/options
    usageStr = 'python calculate_summary_stats.py sequencing_summary.txt sequencing_summary2.txt'
    parser = OptionParser(usage=usageStr, version='%prog 0.0.3')

    #Parse the options/arguments
    options, args = parser.parse_args()

    #Print help message if no input
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    inFile_list = sys.argv[1:]
    read_length = []
    bases = 0
    for inFile in inFile_list:
        if inFile.endswith('.gz'):
            file = gzip.open(inFile,'rt')
        else:
            file = open(inFile, 'r')

        header = file.readline()
        header = header.strip().split()
        length_index = header.index("sequence_length_template")

        for line in file:
            if 'file' in line:  
                continue
            else:
                line = line.strip().split()
                length = int(line[length_index])
                read_length.append(length)
                bases += length
        file.close()

    read_length = sorted(read_length, reverse=True)

    total_bases = numpy.sum(bases)
    total_gigabases = round(total_bases / 1E9, 2)
    target = total_bases / 2.0


    print('File\tread_N50\tGb\tcoverage\t100kb+\t200kb+\t300kb+\t400kb+\t500kb+\t1Mb+\twhales')

    # stats based on human genome (3.3E9 bases per genome)
    coverage = round(total_bases / 3.3E9, 2)
    lt100 = round(sum([i for i in read_length if i >= 100000]) / 3.3E9, 2)
    lt200 = round(sum([i for i in read_length if i >= 200000]) / 3.3E9, 2)
    lt300 = round(sum([i for i in read_length if i >= 300000]) / 3.3E9, 2)
    lt400 = round(sum([i for i in read_length if i >= 400000]) / 3.3E9, 2)
    lt500 = round(sum([i for i in read_length if i >= 500000]) / 3.3E9, 2)
    lt1000 = round(sum([i for i in read_length if i >= 1000000]) / 3.3E9, 2)

    # number of 1 Mb+ reads
    num1000 = len([i for i in read_length if i >= 1000000])

    n50 = 0
    counter = 0
    for item in read_length:
        n50 += item
        counter += 1
        if n50 >= target:
            print(sys.argv[1:], '\t', read_length[counter], '\t', total_gigabases, '\t', coverage, '\t', lt100, '\t', lt200, '\t', lt300, '\t', lt400, '\t', lt500, '\t', lt1000, '\t', num1000, file=sys.stdout)
            break

    print('\ntotal time for the program %.3f' % (time.time()-t0), file=sys.stderr)

if (__name__ == '__main__'):
    main()
    raise SystemExit
