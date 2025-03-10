#!/usr/bin/env python3
import sys
import os
import argparse
from block_utils import Alignment, parseAssemblyIntervals, subtractInterval, getLongDeletionBlocks, runProjectionParallel, splitCandidatesByHaplotype, mergeProjectionDictionaries, BlockList

def match_intervals(interval1, interval2):
    """
    Compare two intervals (each a tuple: (start, end)) and return a label indicating their relationship.
    
    Returns:
      - "no_overlap" if there is no overlap.
      - "perfect" if the intervals are exactly equal.
      - "contained" if the first interval is fully contained in the second.
      - "container" if the second interval is fully contained in the first.
      - "partial" if they overlap but none of the above conditions hold.
    """
    s1, e1 = interval1
    s2, e2 = interval2
    if s1 > e2 or s2 > e1:
        return "no_overlap"
    if s1 == s2 and e1 == e2:
        return "perfect"
    if s1 >= s2 and e1 <= e2:
        return "contained"
    if s2 >= s1 and e2 <= e1:
        return "container"
    return "partial"

def getColCandidates(pafPath: str, faiPath: str, refpafPath: str, indelThreshold: int):
    """
    Process the PAF file and FASTA index and reference PAF file to extract collapsed candidates.
    Returns:
      - allCandidates: list of (chrom, start, end) candidate intervals from final projection.
      - combinedMatches: list of tuples with combined projection information.
          Each tuple has 13 values:
            (p1_chrom, p1_projectable_start, p1_projectable_end, p1_contig, p1_projblock_start, p1_projblock_end,
             p2_chrom, p2_projectable_start, p2_projectable_end, p2_contig, p2_projblock_start, p2_projblock_end,
             match_type)
      - proj1_list: list of dictionaries with projection 1 details.
      - proj2_list: list of dictionaries with projection 2 details.
    """
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
            # Find long deletion blocks; in the coordinates of the assembly
            # and add them to "deletionBlocks"
            deletionBlocks.extend(getLongDeletionBlocks(alignment, indelThreshold))

            # Loss-of-mapping-derived Collapsed Candidates
            #
            # Subtract the alignment from the whole assembly intervals
            # This is for finding the blocks with no alignment to the reference
            refIntervals[alignment.chromName] = subtractInterval(refIntervals[alignment.chromName],
                                                                 (alignment.chromStart, alignment.chromEnd))
            
    # Organize preliminaryCandidates into a dictionary per chromosome
    preliminaryCandidatesDict = {}
    for chrom in refIntervals:
        if chrom not in preliminaryCandidatesDict:
            preliminaryCandidatesDict[chrom] = BlockList()
        for start, end in refIntervals[chrom]:
            preliminaryCandidatesDict[chrom].append((start, end))

    for chrom, start, end in deletionBlocks:
        if chrom not in preliminaryCandidatesDict:
            preliminaryCandidatesDict[chrom] = BlockList()
        preliminaryCandidatesDict[chrom].append((start, end))

    # Merge overlapping blocks in each chromosome
    for chrom, blockList in preliminaryCandidatesDict.items():
        blockList.mergeWithOverlapCount()

    # Convert preliminaryCandidatesDict back to a list of (chrom, start, end, info)
    preliminaryCandidates = [
        (chrom, start, end, info)
        for chrom, blockList in preliminaryCandidatesDict.items()
        for start, end, info in blockList.blocks
    ]

    # Save the list of preliminaryCandidates in a bed file
    print(f"Number of valid preliminaryCandidates after merging: {len(preliminaryCandidates)}")

    # Generate an alignments list from the reference haplotypes alignment. It will be used to project from one reference haplotype to the other 
    alignmentsfirst = []
    with open(refpafPath) as f:
        for line in f:
            blocks = []
            alignment = Alignment(line)
            if alignment.isPrimary == False:
                continue
            alignmentsfirst.append(alignment)
        # Split candidates and create dictionaries to obtain projections.
        # The candidates have to be split so that runProjectionParallel can be run in the two modes based on the haplotype of the candidate (all candidates are from the reference coordinates)
        refCandidates, queryCandidates = splitCandidatesByHaplotype(preliminaryCandidates, alignmentsfirst)

        # Project candidates that have been found on the reference genome/assembly haplotype used as reference in the alignment onto the other haplotype
        projectedrefBlocks = runProjectionParallel(
            alignmentsfirst,
            mode="ref2asm",  # Adjust based on your alignment direction
            blocks=refCandidates,
            includeEndingIndel=False,
            includePostIndel=False,
            threads=4  # Adjust thread count as needed
        )

        # Project candidates that have been found on the genome/assembly haplotype used as query in the alignment onto the other haplotype
        projectedqueryBlocks = runProjectionParallel(
            alignmentsfirst,
            mode="asm2ref",  # Adjust based on your alignment direction
            blocks=queryCandidates,
            includeEndingIndel=False,
            includePostIndel=False,
            threads=4  # Adjust thread count as needed
        )

    print(f"Number of projectedrefBlocks: {sum(len(result[4]) for result in projectedrefBlocks)}")
    print(f"Number of projectedqueryBlocks: {sum(len(result[4]) for result in projectedqueryBlocks)}")
    print(f"projectedrefBlocks: {projectedrefBlocks}")
    print(f"projectedqueryBlocks: {projectedqueryBlocks}")

    # Merge blocks to find projections a second time 
    mergedBlocks = mergeProjectionDictionaries(projectedrefBlocks, projectedqueryBlocks)
    # Save the mergedBlocks output of the first projection in a bed file
    print(f"Number of chromosomes in mergedBlocks: {len(mergedBlocks.keys())}")
    print(f"Total number of mergedBlocks: {sum(len(blocks) for blocks in mergedBlocks.values())}")
    print(f"mergedBlocks: {mergedBlocks}")

    # Convert each key's list in mergedBlocks into a separate BlockList
    mergedBlocksList = {
        key: BlockList([(start, end, info) for start, end, _, _, info in blocks])
        for key, blocks in mergedBlocks.items()
    }
    
    # Merge overlapping blocks in each chromosome
    for chrom, blockList in mergedBlocksList.items():
        blockList.mergeWithOverlapCount()

    # Convert mergedBlocksList back to a dictionary with chrom as keys and lists of (start, end, info)
    mergedDict = {
        chrom: [(start, end, info) for start, end, info in blockList.blocks]
        for chrom, blockList in mergedBlocksList.items()
    }

    # Project merged blocks to the query coordinate system
    alignments = []
    with open(pafPath) as f:
        for line in f:
            alignment = Alignment(line)
            if alignment.isPrimary == False:
                continue
            alignments.append(alignment)

    # Final projection 
    finalProjections = runProjectionParallel(
        alignments,
        mode="ref2asm",
        blocks=mergedDict,
        includeEndingIndel=False,
        includePostIndel=False,
        threads=4
    )

    allcandidates = {}
    for result in finalProjections:
        chromName, contigName, _, projectableBlocks, projectionBlocks, _ = result
        if contigName not in allcandidates:
            allcandidates[contigName] = []
        # Add the "collapsed" info field to each block
        for projectionBlock, projectableBlock in zip(projectionBlocks, projectableBlocks):
            allcandidates[contigName].append((
                projectionBlock[0],  # Start of projected block
                projectionBlock[1],  # End of projected block
                chromName,           # Chromosome name
                projectableBlock[0], # Start of projectable block
                projectableBlock[1]  # End of projectable block
            ))

    # Extract the projections from finalProjections and make them a list that can be converted to a bed file
    allCandidates = []
    for contig, intervals in allcandidates.items():
        if len(intervals) > 0:
            for start, end, chromName, startref, endref in intervals:
                allCandidates.append((contig, start, end))  # Query-only

    # Build projection lists for combining outputs from projection 1 and projection 2
    # Process first projection results (from both projectedrefBlocks and projectedqueryBlocks)
    firstProjectionResults = projectedrefBlocks + projectedqueryBlocks
    proj1_list = []
    for res in firstProjectionResults:
        p1_chrom = res[0]
        p1_contig = res[1]
        p1_projectable = res[3]  # List of tuples: (start, end, info)
        p1_projection = res[4]   # List of tuples: (start, end)
        for pb, pr in zip(p1_projectable, p1_projection):
            proj1_list.append({
                "p1_chrom": p1_chrom,
                "p1_projectable_start": pb[0],
                "p1_projectable_end": pb[1],
                "p1_contig": p1_contig,
                "p1_projblock_start": pr[0],
                "p1_projblock_end": pr[1]
            })

    # Process final projection results
    proj2_list = []
    for res in finalProjections:
        p2_chrom = res[0]
        p2_contig = res[1]
        p2_projectable = res[3]
        p2_projection = res[4]
        for pb, pr in zip(p2_projectable, p2_projection):
            proj2_list.append({
                "p2_chrom": p2_chrom,
                "p2_projectable_start": pb[0],
                "p2_projectable_end": pb[1],
                "p2_contig": p2_contig,
                "p2_projblock_start": pr[0],
                "p2_projblock_end": pr[1]
            })

    # Match entries from projection 1 and projection 2 based on contig and overlapping coordinates.
    # A match is accepted if:
    #   - The candidate contig from projection 1 (p1_contig) equals the chromosome in projection 2 (p2_chrom)
    #   - And the projection block from projection 1 compared with the projectable block from projection 2 are either perfectly equal,
    #     one is contained within the other, or they partially overlap.
    combinedMatches = []
    for entry1 in proj1_list:
        for entry2 in proj2_list:
            if entry1["p1_contig"] == entry2["p2_chrom"]:
                match_type = match_intervals(
                    (entry1["p1_projblock_start"], entry1["p1_projblock_end"]),
                    (entry2["p2_projectable_start"], entry2["p2_projectable_end"])
                )
                if match_type != "no_overlap":
                    combinedMatches.append((
                        entry1["p1_chrom"],
                        entry1["p1_projectable_start"],
                        entry1["p1_projectable_end"],
                        entry1["p1_contig"],
                        entry1["p1_projblock_start"],
                        entry1["p1_projblock_end"],
                        entry2["p2_chrom"],
                        entry2["p2_projectable_start"],
                        entry2["p2_projectable_end"],
                        entry2["p2_contig"],
                        entry2["p2_projblock_start"],
                        entry2["p2_projblock_end"],
                        match_type
                    ))

    return allCandidates, combinedMatches, proj1_list, proj2_list

def main():
    parser = argparse.ArgumentParser(description='A program for extracting collapsed candidates in a draft assembly. It needs the alignments of the assembly contigs to the reference (or high quality assembly) and the alignment of the reference haplotypes. The alignment files should be sorted and be in the PAF format. The output is a BED file that needs to be sorted and merged, along with additional CSV files containing projection information.')
    parser.add_argument('--paf', type=str,
                        help='(PAF format) The alignments of the assembly to the reference. It should include the cigar format.')
    parser.add_argument('--fai', type=str,
                        help='(Fasta index) The fasta index of the reference (not the assembly). It will be used for extracting the parts of the chromosomes with no mapping')
    parser.add_argument('--refpaf', type=str,
                        help='(PAF format) The alignments of the two reference haplotypes. It should include the cigar format. Needed to map deletions found in one of the query haplotypes, projecting first on the other reference haplotype.')
    parser.add_argument('--output', type=str,
                        help=('Output directory where all files will be saved. The main candidate BED file '
                              'will be saved in this directory and additional files in a subdirectory.'))
    parser.add_argument('--indel', type=int,
                        help='The threshold for indel-based candidates (default = 100)', default=100)
    
    # Fetch the arguments
    args = parser.parse_args()
    pafPath = args.paf
    faiPath = args.fai
    refpafPath = args.refpaf
    outDir = args.output
    indelThreshold = args.indel

    # Create the main output directory if it does not exist.
    if not os.path.exists(outDir):
        os.makedirs(outDir)

    # Use the basename of the output directory as the file name prefix.
    main_prefix = os.path.basename(os.path.normpath(outDir))

    # Create a subdirectory for additional files.
    additional_dir = os.path.join(outDir, "additional_files")
    if not os.path.exists(additional_dir):
        os.makedirs(additional_dir)

    # Run the candidate extraction.
    colCandidates, combinedMatches, proj1_list, proj2_list = getColCandidates(pafPath, faiPath, refpafPath, indelThreshold)
    
    # Write the main candidate BED file (colCandidates) in the main output directory.
    main_bed_path = os.path.join(outDir, f"{main_prefix}.bed")
    with open(main_bed_path, "w") as f:
        for chrom, start, end in colCandidates:
            f.write("{}\t{}\t{}\n".format(chrom, start, end))
    
    # Save the combined projection matches to a CSV file
    combined_csv_path = os.path.join(additional_dir, f"{main_prefix}_combined_projections.csv")
    with open(combined_csv_path, "w") as f:
        # Write header
        f.write("p1_chrom,p1_projectable_start,p1_projectable_end,p1_contig,p1_projblock_start,p1_projblock_end,"
                "p2_chrom,p2_projectable_start,p2_projectable_end,p2_contig,p2_projblock_start,p2_projblock_end,match_type\n")
        for row in combinedMatches:
            f.write(",".join(str(x) for x in row) + "\n")
    
    # Save projection 1 details to a CSV file
    proj1_csv_path = os.path.join(additional_dir, f"{main_prefix}_projection1.csv")
    with open(proj1_csv_path, "w") as f:
        f.write("p1_chrom,p1_projectable_start,p1_projectable_end,p1_contig,p1_projblock_start,p1_projblock_end\n")
        for entry in proj1_list:
            f.write("{},{},{},{},{},{}\n".format(entry["p1_chrom"], entry["p1_projectable_start"], entry["p1_projectable_end"],
                                                 entry["p1_contig"], entry["p1_projblock_start"], entry["p1_projblock_end"]))
    # Save projection 2 details to a CSV file
    proj2_csv_path = os.path.join(additional_dir, f"{main_prefix}_projection2.csv")
    with open(proj2_csv_path, "w") as f:
        f.write("p2_chrom,p2_projectable_start,p2_projectable_end,p2_contig,p2_projblock_start,p2_projblock_end\n")
        for entry in proj2_list:
            f.write("{},{},{},{},{},{}\n".format(entry["p2_chrom"], entry["p2_projectable_start"], entry["p2_projectable_end"],
                                                 entry["p2_contig"], entry["p2_projblock_start"], entry["p2_projblock_end"]))

if __name__ == "__main__":
    main()
