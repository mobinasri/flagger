'''
Purpose: parse paf file to identify regions of 100% homozygosity bigger than a window size of m, output bedfile
Author: Mira Mastoras, mmastora@ucsc.edu
Usage: python3 find_homozygous_regions.py -p paf_file -m min_length -e extend_windows -o output_prefix
'''

import argparse
from block_utils import Alignment

def arg_parser():
    '''
    Parses command line arguments with argparse
    '''
    parser = argparse.ArgumentParser(
        prog='',
        description="")

    parser.add_argument("-p", "--paf_file",
                        required=True,
                        help="paf file of aligments between two haplotype assemblies. --eqx flag must be used in alignment")
    parser.add_argument("-m", "--min_length",
                        required=False,
                        default=int(1000),
                        help="Minimum window size in bp for homozygous regions")
    parser.add_argument("-e", "--extend_windows",
                        required=False,
                        help="(Optional) Number of bp to extend windows by, to capture surrounding heterozygosity")
    parser.add_argument("-o", "--out_bed",
                        required=True,
                        help="Output prefix")

    return parser.parse_args()

def find_homozygous_regions(paf_line, min_length):
    '''

    :param pafLine: line from paf file
    :param minLength: Minimum window size for homozygous regions
    :return: list of 100% homozygous regions with each element in list as [(ref chr, ref start, ref end) , (query chr, query start, query end)]
    '''
    # create alignment object from paf line
    alignment = Alignment(paf_line)

    # initialize first block with position of alignment start on ref (chrom) and query (contig)
    block_chrom_start= alignment.chromStart

    # check orientation of contig
    # if it is positive, alignment starts from contig start
    if alignment.orientation == "+":
        multiplier = 1
        block_contig_start = alignment.contigStart

    # if it is negative, alignment starts from contig end, and we need a -1 multiplier to subtract cigar[1] (number of bases) from contig end
    else:
        multiplier = -1
        block_contig_start = alignment.contigEnd

    hom_regions=[]

    for cigar in alignment.cigarList: # cigar = ("M | D | I | X | = ", num bases)
        if cigar[0]=='=': # query exactly matches ref
            # if there's a match, check if we want to record this region
            if cigar[1] >= min_length:

                if alignment.orientation=="+":
                    hom_block_contig = (alignment.contigName, block_contig_start, block_contig_start + (multiplier * cigar[1]))

                # for - orientation, need to swap start and end coords so bed format is correct (start < end )
                else:
                    hom_block_contig = (alignment.contigName, block_contig_start + (multiplier * cigar[1]), block_contig_start)

                hom_block_chrom = (alignment.chromName, block_chrom_start, block_chrom_start +  cigar[1])
                hom_regions.append([hom_block_chrom, hom_block_contig])

            # move start by # bases mismatch, for both ref and query
            block_chrom_start = block_chrom_start + cigar[1]
            block_contig_start = block_contig_start + (multiplier * cigar[1])

        if cigar[0] =="D": # query has deletion relative to ref
            block_chrom_start = block_chrom_start +  cigar[1]

        if cigar[0] == "I":  # query has insertion relative to ref
            block_contig_start = block_contig_start + (multiplier * cigar[1])

        if cigar[0] =="X": # query has a mismatch (snp)
            block_chrom_start = block_chrom_start + cigar[1]
            block_contig_start = block_contig_start + (multiplier * cigar[1])

    return(hom_regions)

def main():
    args = arg_parser()

    # get arguments
    paf_path=args.paf_file
    min_length=int(args.min_length)

    # for each line in paf file, find homozygous regions
    hom_regions=[]
    with open(paf_path, "r") as fpaf:
        for line in fpaf:
            hom_regions.append(find_homozygous_regions(line, min_length))

    out_bed= args.out_bed + ".bed"
    out_flanking_bed = args.out_bed + "flanking_" + args.extend_windows+ ".bed"

    # write output bed file
    with open(out_bed, 'w') as out:
        for line in hom_regions:
            for reg in line:
                print(reg[0][0], max(0,reg[0][1]), reg[0][2] , reg[1][0], max(0,reg[1][1]), reg[1][2], sep="\t", file=out)

    # write optional bedfile with windows extended by -e bp into flanking regions
    if args.extend_windows:
        with open(out_flanking_bed, 'w') as out:
            for line in hom_regions:
                for reg in line:
                    print(reg[0][0], max(0,reg[0][1] - int(args.extend_windows)), reg[0][2] + int(args.extend_windows), reg[1][0], max(0,reg[1][1] - int(args.extend_windows)), reg[1][2] + int(args.extend_windows), sep="\t", file=out)


if __name__ == '__main__':
    main()