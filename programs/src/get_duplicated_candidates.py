#!/usr/bin/env python3
import sys
import argparse
import os
from block_utils import findProjections, Alignment, parseAssemblyIntervals, subtractInterval, getLongInsertionBlocks

def getDupCandidates(pafPath: str, faiPath: str, indelThreshold: int):
    """
    Process the PAF file and assembly FASTA index to extract candidate duplicated regions.
    Returns:
      - allCandidates: list of (contig, start, end) from the union of loss-of-mapping, overlap- and insertion-derived candidates.
      - projectionDetails: list of tuples with projection information.
      - overlapBlocks: list of (contig, start, end) from overlap-derived candidates.
      - assemblyIntervals_list: list of (contig, start, end) corresponding to intervals with no mapping.
      - insertionBlocks: list of (contig, start, end) from insertion-derived candidates.
      - problematicBlocks: list of tuples (alignment_type, chrom, ref_overlap_start, ref_overlap_end, contig, contig_start, contig_end)
         for which the projection was not present.
    """
    preAlignment = None
    insertionBlocks = []
    overlapBlocks = []
    problematicBlocks = []
    projectionDetails = []  # Each record: (alignment_type, contig, chrom, ref_overlap_start, ref_overlap_end, projectable, projected, cigar)

    # Parse the assembly intervals from the .fai file.
    assemblyIntervals = parseAssemblyIntervals(faiPath)

    with open(pafPath) as f:
        for line in f:
            alignment = Alignment(line)
            if not alignment.isPrimary:
                continue

            # Insertion-derived Duplicated candidates:
            # Find long insertion blocks (in assembly coordinates) and add them to "insertionBlocks"
            insertionBlocks.extend(getLongInsertionBlocks(alignment, indelThreshold))

            # Loss-of-mapping-derived Duplicated candidates:
            # Subtract the alignment from the whole assembly intervals to get regions with no mapping.
            assemblyIntervals[alignment.contigName] = subtractInterval(
                assemblyIntervals[alignment.contigName],
                (alignment.contigStart, alignment.contigEnd)
            )

            # Overlap-derived Duplicated candidates:
            # Check if the current alignment has overlap with the previous alignment.
            # If there is an overlap, extract the overlap (in reference coordinates) and try to project it to the assembly.
            if preAlignment is None:
                preAlignment = alignment
                continue
            else:
                # If the reference chromosome didn't change and
                # the start of the current alignment is before the 
                # end of the previous alignment
                if alignment.chromName == preAlignment.chromName and \
                   alignment.chromStart < preAlignment.chromEnd:
                    # Define the overlap interval in the reference coordinates (1-based start).
                    refOverlapInterval = (alignment.chromStart + 1, min(alignment.chromEnd, preAlignment.chromEnd), "NA")
                    
                    # Process the current alignment:
                    curr_proj, curr_projections, curr_cigar = findProjections(
                        'ref2asm',
                        alignment.cigarList,
                        [refOverlapInterval],
                        alignment.chromLength,
                        alignment.chromStart + 1, alignment.chromEnd,
                        alignment.contigLength,
                        alignment.contigStart + 1, alignment.contigEnd,
                        alignment.orientation,
                        False, False
                    )
                    # There is only one projection; projections[0]
                    # that is actually the projection of [refOverlapInterval]
                    curr_projectable = curr_proj[0] if curr_proj else None
                    curr_projected = curr_projections[0] if curr_projections else None
                    curr_cigar_item = curr_cigar[0] if curr_cigar else None
                    projectionDetails.append((
                        "current",
                        alignment.contigName,
                        alignment.chromName,
                        refOverlapInterval[0],
                        refOverlapInterval[1],
                        curr_projectable,
                        curr_projected,
                        curr_cigar_item
                    ))
                    if curr_projections:
                        # Save the projected interval (convert from 1-based to 0-based start).
                        overlapBlocks.append((alignment.contigName, curr_projected[0] - 1, curr_projected[1]))
                    else:
                        # Save the overlap interval and alignment coordinates if projection is missing.
                        problematicBlocks.append((
                            alignment.chromName,
                            refOverlapInterval[0],
                            refOverlapInterval[1],
                            alignment.contigName,
                            alignment.contigStart,
                            alignment.contigEnd,
                            preAlignment.contigName,
                            preAlignment.contigStart,
                            preAlignment.contigEnd
                        ))
                    
                    # Process the previous alignment:
                    pre_proj, pre_projections, pre_cigar = findProjections(
                        'ref2asm',
                        preAlignment.cigarList,
                        [refOverlapInterval],
                        preAlignment.chromLength,
                        preAlignment.chromStart + 1, preAlignment.chromEnd,
                        preAlignment.contigLength,
                        preAlignment.contigStart + 1, preAlignment.contigEnd,
                        preAlignment.orientation,
                        False, False
                    )
                    pre_projectable = pre_proj[0] if pre_proj else None
                    pre_projected = pre_projections[0] if pre_projections else None
                    pre_cigar_item = pre_cigar[0] if pre_cigar else None
                    projectionDetails.append((
                        "previous",
                        preAlignment.contigName,
                        preAlignment.chromName,
                        refOverlapInterval[0],
                        refOverlapInterval[1],
                        pre_projectable,
                        pre_projected,
                        pre_cigar_item
                    ))
                    if pre_projections:
                        overlapBlocks.append((preAlignment.contigName, pre_projected[0] - 1, pre_projected[1]))
                    else:
                        problematicBlocks.append((
                            alignment.chromName,
                            refOverlapInterval[0],
                            refOverlapInterval[1],
                            alignment.contigName,
                            alignment.contigStart,
                            alignment.contigEnd,
                            preAlignment.contigName,
                            preAlignment.contigStart,
                            preAlignment.contigEnd
                        ))
                # Update the previous alignment to the current one.
                preAlignment = alignment

    # Convert assemblyIntervals (a dict) to a flat list of intervals.
    assemblyIntervals_list = []
    for contig, intervals in assemblyIntervals.items():
        for start, end in intervals:
            assemblyIntervals_list.append((contig, start, end))

    # Combine all candidate intervals.
    allCandidates = []
    allCandidates.extend(assemblyIntervals_list)
    allCandidates.extend(overlapBlocks)
    allCandidates.extend(insertionBlocks)

    return allCandidates, projectionDetails, overlapBlocks, assemblyIntervals_list, insertionBlocks, problematicBlocks

def main():
    parser = argparse.ArgumentParser(
        description=('A program for extracting duplicated candidates in a draft assembly. '
                     'It requires the alignments (in PAF format, including CIGAR strings) of the assembly contigs '
                     'to a reference (or high-quality assembly), and the assembly's FASTA index (.fai). '
                     'The main candidate BED file and additional files (projection details, overlap, assembly, insertion, '
                     'and missing projection intervals) will be saved in the specified output directory.')
    )
    parser.add_argument('--paf', type=str,
                        help='(PAF format) The alignments of the assembly to the reference. Should include the CIGAR format.')
    parser.add_argument('--fai', type=str,
                        help='(Fasta index) The FASTA index of the assembly (not the reference).')
    parser.add_argument('--output', type=str,
                        help=('Output directory where all files will be saved. The main candidate BED file '
                              'will be saved in this directory and additional files in a subdirectory.'))
    parser.add_argument('--indel', type=int,
                        help='Threshold for indel-based candidates (default = 100)', default=100)
    
    args = parser.parse_args()
    pafPath = args.paf
    faiPath = args.fai
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
    (allCandidates, projectionDetails, overlapBlocks,
     assemblyIntervals_list, insertionBlocks, problematicBlocks) = getDupCandidates(pafPath, faiPath, indelThreshold)

    # Write the main candidate BED file (allCandidates) in the main output directory.
    main_bed_path = os.path.join(outDir, f"{main_prefix}.bed")
    with open(main_bed_path, "w") as f:
        for contig, start, end in allCandidates:
            f.write(f"{contig}\t{start}\t{end}\n")

    # Write the projection details CSV file.
    projection_csv_path = os.path.join(additional_dir, f"{main_prefix}_projection_blocks.csv")
    with open(projection_csv_path, "w") as f:
        # Write a header line.
        f.write("alignment_type,contig,chrom,ref_overlap_start,ref_overlap_end,projectable,projected,cigar\n")
        for record in projectionDetails:
            # Each record is a tuple:
            # (alignment_type, contig, chrom, ref_overlap_start, ref_overlap_end, projectable, projected, cigar)
            alignment_type, contig, chrom, ref_start, ref_end, projectable, projected, cigar = record
            proj_str = f"{projectable[0]}-{projectable[1]}" if projectable else ""
            projd_str = f"{projected[0]}-{projected[1]}" if projected else ""
            f.write(f"{alignment_type},{contig},{chrom},{ref_start},{ref_end},{proj_str},{projd_str},{cigar if cigar else ''}\n")

    # Write the overlapBlocks BED file.
    overlap_bed_path = os.path.join(additional_dir, f"{main_prefix}_overlapBlocks.bed")
    with open(overlap_bed_path, "w") as f:
        for contig, start, end in overlapBlocks:
            f.write(f"{contig}\t{start}\t{end}\n")

    # Write the assembly intervals (loss-of-mapping) BED file.
    assembly_bed_path = os.path.join(additional_dir, f"{main_prefix}_assemblyIntervals.bed")
    with open(assembly_bed_path, "w") as f:
        for contig, start, end in assemblyIntervals_list:
            f.write(f"{contig}\t{start}\t{end}\n")

    # Write the insertionBlocks BED file.
    insertion_bed_path = os.path.join(additional_dir, f"{main_prefix}_insertionBlocks.bed")
    with open(insertion_bed_path, "w") as f:
        for contig, start, end in insertionBlocks:
            f.write(f"{contig}\t{start}\t{end}\n")

    # Write the problematic (missing projection) CSV file.
    missing_projection_csv_path = os.path.join(additional_dir, f"{main_prefix}_missing_projection_blocks.csv")
    with open(missing_projection_csv_path, "w") as f:
        # Header: alignment_type, chrom, ref_overlap_start, ref_overlap_end, contig, contig_start, contig_end
        f.write("overlap_ref_chrom,ref_overlap_start,ref_overlap_end,contig_align,contig_align_start,contig_align_end,contig_prealign,contig_prealign_start,contig_prealign_end\n")
        for pb in problematicBlocks:
            # pb is a tuple: (overlap_ref_chrom,ref_overlap_start,ref_overlap_end,contig_align,contig_align_start,contig_align_end,contig_prealign,contig_prealign_start,contig_prealign_end)
            f.write("{},{},{},{},{},{},{},{},{}\n".format(pb[0], pb[1], pb[2], pb[3], pb[4], pb[5], pb[6], pb[7], pb[8]))

if __name__ == "__main__":
    main()
