import sys
import argparse
from block_utils import findProjections, Alignment, parseAssemblyIntervals, subtractInterval, getLongDeletionBlocks

def getColCandidates(pafPath: str , faiPath: str, indelThreshold: int):
    deletionBlocks = []
    refIntervals = parseAssemblyIntervals(faiPath)
    with open(pafPath) as f:
        for line in f:
            blocks = []
            alignment = Alignment(line)
            if alignment.isPrimary == False:
                continue

            # Deletion-derived Collapsed Candidates
            #
            # Find long insertion blocks; in the coordinates of the assembly
            # and add them to "insertionBlocks"
            deletionBlocks.extend(getLongDeletionBlocks(alignment, indelThreshold))

            # Loss-of-mapping-derived Collapsed Candidates
            #
            # Subtract the alignment from the whole assembly intervals
            # This is for finding the blocks with no alignment to the reference
            refIntervals[alignment.chromName] = subtractInterval(refIntervals[alignment.chromName], 
                                                                 (alignment.chromStart, alignment.chromEnd))
 
    
    allCandidates = [] # candidates may have redundant overlaps should be merged by bedtools merge

    for chrom in refIntervals:
        intervals = refIntervals[chrom]
        if len(intervals) > 0:
            for start, end in intervals:
                allCandidates.append((chrom, start, end))
    allCandidates.extend(deletionBlocks)
    return allCandidates



def main():
    parser = argparse.ArgumentParser(description='A program for extracting collapsed candidates in a draft assembly. It needs the alignments of the assembly contigs to the reference (or high quality assembly). The alignment file should be sorted and be in the PAF format. The output is a BED file needs to be sorted and merged.')
    parser.add_argument('--paf', type=str,
                    help='(PAF format) The alignments of the assembly to the reference. It should include the cigar format.')
    parser.add_argument('--fai', type=str,
                    help='(Fasta index) The fasta index of the reference (not the assembly). It will be used for extracting the parts of the chromosomes with no mapping')
    parser.add_argument('--output', type=str,
                    help='(BED format) A path for saving the candidate collapsed blocks in the assembly. Note that the lines may not be sorted and may have overlaps.')
    parser.add_argument('--indel', type=int,
                    help='The threshold for indel-based candidates (default = 100)', default=100)
    
    # Fetch the arguments
    args = parser.parse_args()
    pafPath = args.paf
    faiPath = args.fai
    outPath = args.output
    indelThreshold = args.indel

    colCandidates = getColCandidates(pafPath, faiPath, indelThreshold)
    with open(outPath, "w+") as f:
        for chrom, start, end in colCandidates:
            f.write("{}\t{}\t{}\n".format(chrom, start, end))

main()
