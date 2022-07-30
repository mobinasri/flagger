import sys
import argparse
from block_utils import findProjections, Alignment, parseAssemblyIntervals, subtractInterval, getLongInsertionBlocks

def getDupCandidates(pafPath: str , faiPath: str, indelThreshold: int):
    preAlignment = None
    insertionBlocks = []
    overlapBlocks = []
    assemblyIntervals = parseAssemblyIntervals(faiPath)
    with open(pafPath) as f:
        for line in f:
            blocks = []
            alignment = Alignment(line)
            if alignment.isPrimary == False:
                continue

            # Insertion-derived Duplicated Candidates
            #
            # Find long insertion blocks; in the coordinates of the assembly
            # and add them to "insertionBlocks"
            insertionBlocks.extend(getLongInsertionBlocks(alignment, indelThreshold))

            # Loss-of-mapping-derived Duplicated Candidates
            #
            # Subtract the alignment from the whole assembly intervals
            # This is for finding the blocks with no alignment to the reference
            assemblyIntervals[alignment.contigName] = subtractInterval(assemblyIntervals[alignment.contigName], 
                                                                       (alignment.contigStart, alignment.contigEnd))

            # Overlap-derived Duplicated Candidates
            #
            # Check if the current alignment has overlap with the previous alignment
            # If there is an overlap extract the overlap in the coordinates of the assembly
            if preAlignment == None:
                preAlignment = alignment
                continue
            else:
                # If the reference chromosome didn't change and
                # the start of the current alignmet is before the 
                # end of the previous alignment
                if alignment.chromName == preAlignment.chromName and \
                   alignment.chromStart < preAlignment.chromEnd:
                        # Get the overlap in the reference coordiantes
                        refOverlapInterval = (alignment.chromStart + 1, min(alignment.chromEnd, preAlignment.chromEnd))

                        # Find the projection of the overlap in the assembly coordinates
                        # For the current alignment
                        projectables, projections = findProjections('ref2asm', 
                                                           alignment.cigarList, 
                                                           [refOverlapInterval],
                                                           alignment.chromLength, 
                                                           alignment.chromStart + 1, alignment.chromEnd, # make 1-based start
                                                           alignment.contigLength, 
                                                           alignment.contigStart + 1, alignment.contigEnd, # make 1-based start
                                                           alignment.orientation)
                        # There is only one projection; projections[0]
                        # that is actually the projection of [refOverlapInterval]
                        overlapBlocks.append((alignment.contigName, projections[0][0] - 1, projections[0][1]))

                        # Find the projection of the overlap in the assembly coordinates
                        # For the previous alignment
                        projectables, projections = findProjections('ref2asm', 
                                                           preAlignment.cigarList, 
                                                           [refOverlapInterval],
                                                           preAlignment.chromLength, 
                                                           preAlignment.chromStart + 1, preAlignment.chromEnd, # make 1-based start
                                                           preAlignment.contigLength, 
                                                           preAlignment.contigStart + 1, preAlignment.contigEnd, # make 1-based start
                                                           preAlignment.orientation)
                        # There is only one projection; projections[0]
                        # that is actually the projection of [refOverlapInterval]
                        overlapBlocks.append((preAlignment.contigName, projections[0][0] - 1, projections[0][1]))
                # Update the previous alignment
                preAlignment = alignment
    
    
    allCandidates = [] # candidates may have redundant overlaps should be merged by bedtools merge

    for contig in assemblyIntervals:
        intervals = assemblyIntervals[contig]
        if len(intervals) > 0:
            for start, end in intervals:
                allCandidates.append((contig, start, end))
    allCandidates.extend(overlapBlocks)
    allCandidates.extend(insertionBlocks)
    return allCandidates



def main():
    parser = argparse.ArgumentParser(description='A program for extracting collapsed candidates in a draft assembly. It needs the alignments of the assembly contigs to the reference (or high quality assembly). The alignment file should be sorted and be in the PAF format. The output is a BED file needs to be sorted and merged.')
    parser.add_argument('--paf', type=str,
                    help='(PAF format) The alignments of the assembly to the reference. It should include the cigar format.')
    parser.add_argument('--fai', type=str,
                    help='(Fasta index) The fasta index of the assembly (not the reference). It will be used for extracting the parts of the contigs with no mapping')
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

    dupCandidates = getDupCandidates(pafPath, faiPath, indelThreshold)
    with open(outPath, "w+") as f:
        for contig, start, end in dupCandidates:
            f.write("{}\t{}\t{}\n".format(contig, start, end))

main()
