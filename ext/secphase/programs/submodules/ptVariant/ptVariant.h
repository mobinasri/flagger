#ifndef PT_VARIANT_H
#define PT_VARIANT_H

#include "sam.h"
#include "faidx.h"
#include <time.h>
#include "bgzf.h"
#include <regex.h>
#include "sonLib.h"
#include <assert.h>
#include <math.h>
#include <float.h>
#include "cigar_it.h"
#include "common.h"
#include "vcf.h"
#include "edlib.h"
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "ptBlock.h"
#include "ptAlignment.h"


/*! @typedef
 @abstract Structure for a variant record
 @field contig		Name of the contig where the variant is located
 @field pos		    Position of the variant in contig
 @field type		Type of the variant;
                    (Look at https://github.com/samtools/htslib/blob/develop/htslib/vcf.h)
 @field vaf		    Variant allele frequency
 @field gq		    Genotype quality
 @field gt          Genotype list
 @field gt_len		Length of gt
 @field alleles		allele sequences (starts with the REF allele)
 @field longest_allele_len	Length of the longest allele (could be insertion, deletion or snp)
 @field ps              Phase block number
 */
typedef struct {
    char contig[50];
    int32_t pos;
    int8_t type;
    float vaf;
    int8_t gq;
    int8_t *gt;
    int32_t gt_len;
    stList *alleles;
    int32_t longest_allele_len;
    int64_t ps;
} ptVariant;


/// Create a ptVariant structure
/**
 * Parameters are listed in the definition of the structure
 * @return 	Pointer to a new variant struct on success; NULL on failure
 * 		The ptMarker struct returned by a successful call should be freed
 * 		via ptMarker_destruct() when it is no longer needed.
 */
ptVariant *ptVariant_construct(char *contig, int32_t pos, int8_t type, float vaf, int8_t gq, int64_t ps);


/// Free a ptVariant structure
/*
 * @param variant       Pointer to the variant record
 */
void ptVariant_destruct(ptVariant *variant);


int ptVariant_cmp(const void *a, const void *b);


/// Add an allele sequence to a variant record
/**
   @param variant	Pointer to the variant record
   @param allele	Allele sequence to be added

   Note that alleles should be added in the correct order
   starting from the REF allele to be consistent with the
   genotype indices
 */
void ptVariant_append_allele(ptVariant *variant, char *allele);


/// Add a genotype to a variant record
/**
   @param variant       Pointer to the variant record
   @param gt        	genotype index to be added
 */
void ptVariant_append_gt(ptVariant *variant, int8_t gt);


ptVariant *ptVariant_copy(ptVariant *src);


// Returns true if two variants are at the same location (and of course same contig/chrom) otherwise false
bool ptVariant_is_equal(ptVariant *var1, ptVariant *var2);


//Returns true if variant exists in the list otherwise false
bool ptVariant_exist_in_list(ptVariant *var, stList *var_list);


void ptVariant_print(ptVariant *variant);

/// Swap the genotypes of the variants in a phase block
// It is assumed that the variants are sorted by phase block
/**
   @param variants	A list of variant records
   @param start		The index of the first variant in the phase block
   @param end		The index of the last variant in the phase block
*/
void ptVariant_swap_gt(stList *variants, int start, int end);


/// Parse variants from a vcf file and save the phased variants
// in a list. Each variant is saved in a ptVariant structure.
// If consistent_gt is True this function will check if the 
// first genotypes of all the phased variants in a phase block 
// are more consistent with the reference allele than the second
// genotypes.If it was not consistent then all the first and 
// second genotypes will be swapped.
// It can be useful for polishing purposes when we assume that 
// the first genotypes show the changes that has to be made on 
// the assembly sequence.
// Only the variants with two genotypes (diploid genome) are included.
/**
   @param vcf_path	The path to a vcf file
*/
stList *read_phased_variants(char *vcf_path, bool consistent_gt, int min_gq);


// 
// Make a new list of variants after removing the variants
// whose first genotypes are 0, which means that their first 
// genotypes are same as reference
//
/*
 * @param variants	stList of variants (can be the output of read_phased_variants())
 *
 */
stList *filter_ref_variants(stList *variants);


// 
// Make a table of blocks that surround the given variants.
// The keys of the table are contig names and the values
// are related lists of variant blocks.
// 
// This function first makes a symmetric window around each variant
// The length of window is (2 * min_margin + 1) initially
// If min_margin is smaller than the length of variant then the
// length of the window will increase to (2 * longest_allele_len + 1)
// This is to make sure that we take enough flanking sequence to
// have a reliable alignment between reference and read sequence
// (which will be performed later by edlib functions)
// 
// If the blocks of close variants have overlap they will be merged
// into one block. 
// The data attribute of each block contains a stList of variants 
// that it spans over. This information is neccessary for polishing the
// reference sequence later.
//
// Note that variants should be sorted by contig name and then by position
// ptVariant_cmp and stList_sort can be used for doing so
//
/*
 * @param variants      	stList of variants (can be the output of filter_ref_variants())
 * @param fai			Fasta index (needed for not exceeding the length of contigs)
 * @param min_margin		Minimum margin length from each side of a variant
 * @return variant_blocks	Table of variant blocks (keys=contig names)
 */
stHash *ptVariant_extract_variant_blocks(stList *variants, const faidx_t *fai, int min_margin);


//
// Fetch the read sequence of a given block
// The attributes sqs and sqe should have been set in the block
//
/*
 * @param alignment	Alignment saved as a ptAlignment struct
 * @param block		Block saved as a ptBlock struct
 * @return block_seq	Read sequence of the given block
 */
char *fetch_read_seq(ptAlignment *alignment, ptBlock *block);


//
// Fetch the corrected reference sequence of a given block
// The variants in the block will be applied to the reference sequence.
// The attributes rfs and rfe should have been set in the block.
//
/*
 * @param fai			Fasta index
 * @param block         	Block saved as a ptBlock struct, data attribute should contain 
 * 				the list of variants in that block
 * @param conting_name		Contig name
 * @return seq_corrected	Corrected reference sequence of the block    
 */

char *fetch_corrected_ref_seq(const faidx_t *fai, ptBlock *block, char *contig_name);


//
// Project blocks onto read coordinates
// Make a new list of blocks. For each block rds_f and rde_f attributes are set after 
// projecting the reference coordinates to read by iterating over cigar operations.
// The new blocks contain the variant information
// Note that the output blocks of this function may have overlaps on read coordinates.
// To merge the overlapping blocks merge_variant_read_blocks() can be called.
//
/*
 * @param alignment		Alignment saved as a ptAlignment struct
 * @param ref_variant_blocks	stList of blocks with rfs and rfe
 * @return read_blocks		stList of blocks with rds_f and rde_f
 */
// ref_variant_blocks should be sorted by rfs with no overlap
stList *project_blocks_to_read(ptAlignment *alignment, stList *ref_variant_blocks);


/**
 * Merge overlapping blocks in the given list of variant blocks
 * When two overlapping blocks are merged all the variants in the
 * two blocks will be added to the final merged block
   @param blocks     a list of blocks sorted by rds_f
   @return      A list of merged blocks
 */
stList *merge_variant_read_blocks(stList *blocks);


//
// Extract read sequence and corrected reference sequence 
// for the given block and then compute the edit distance 
// between the two sequences.
//
/*
 * @param alignment     Alignment saved as a ptAlignment struct
 * @param fai		Fasta index
 * @param block         Block saved as a ptBlock struct
 * @param contig_name	Contig name
 * @return edit_distance	Edit distance between read sequence and corrected ref sequence
 */
int get_edit_distance(ptAlignment *alignment, faidx_t *fai, ptBlock *block);



// Project read blocks to reference by iterating over the alignment 
// cigar operations.
// Extract read sequence and corrected reference sequence for each 
// block and then compute the edit distance between the two sequences. 
// The output of this function would be the summation of edit distances 
// between read and reference for all blocks.
// Note that the input blocks should be in read coordinates (recommended to be 
// the output of merge_variant_read_blocks())
//
/*
 * @param alignment     	Alignment saved as a ptAlignment struct
 * @param fai           	Fasta index
 * @param variant_read_blocks	stList of blocks each of which saved as a ptBlock struct
 * @param contig_name   	Contig name
 * @return edit_distance	Summation of edit distances between read and reference for all blocks
 */
int get_total_edit_distance(ptAlignment *alignment, const faidx_t *fai, char *contig_name,
                            stList *variant_read_blocks, stList *projected_blocks);


// This function performs three main tasks:
// 1. Parse all phased variants from the vcf file
// 2. Keep the variants whose first genotypes are different from reference (So can be used for polishing)
// 3. Construct blocks around the selected variants and save them in a stHash table (keys are contig names and values
//    are lists of blocks saved as ptBlock. The data attribute in each block is a stList of the variants in that block)
//
/*
 * @param vcf_path		Path to VCF file with phased variants
 * @param fai                   Fasta index
 * @param min_margin		Minimum size of the block surrounding each variant
 * @param min_gq		Minimum genotype quality
 *
 */
stHash *ptVariant_parse_variants_and_extract_blocks(char *vcf_path, char *bed_path, faidx_t *fai, int min_margin, int min_gq);


//
// Save the variant reference blocks in a bed file
//
/*
 * @param variant_ref_blocks	Output of parse_variants_and_extract_blocks()
 * @param bed_path		Path to the bed file for saving the blocks
 *
 */
void ptVariant_save_variant_ref_blocks(stHash *variant_ref_blocks, char *bed_path);



// 
// 1. Project variant reference blocks to the read coordinates per alignment
// 2. Intersect and merge blocks in the read coordinates
//
//
/*
 * @param variant_ref_blocks    Output of parse_variants_and_extract_blocks()
 * @param alignments            An array of alignments each of which saved as ptAlignment struct
 * @param alignments_len        Length of alignments array
 * @param sam_hdr		SAM header 
 *
 */
stList *
ptVariant_get_merged_variant_read_blocks(stHash *variant_ref_blocks, ptAlignment **alignments, int alignments_len);


/**
 * This function does these main tasks
 *  1. Project read blocks to reference
 *  2. Calculate the edit distance between corrected reference and read sequence
 *  3. Set the alignment scores as negative values of edit distances (The best alignment should have the highest score)
 *
 * @param variant_ref_blocks    Output of parse_variants_and_extract_blocks()
 * @param alignments            An array of alignments each of which saved as ptAlignment struct
 * @param alignments_len        Length of alignments array
 * @param sam_hdr		SAM header
 * @return Array of stList. Each element is a list of proejected blocks for each alignment containing the start and end coordinates on both ref and read
 *
 */
stList **
set_scores_as_edit_distances(stList *read_blocks_merged, ptAlignment **alignments, int alignments_len, faidx_t *fai);

/**
 * If there is at least one alignment that encompasses a variant block
 * completely then return true otherwise false
 *
 *
 * @param alignments                        An array of alignments each of which saved as ptAlignment struct
 * @param alignments_len                    Length of alignments array
 * @param variant_ref_blocks_per_contig     Output of parse_variants_and_extract_blocks()
 */
bool overlap_variant_ref_blocks(stHash *variant_ref_blocks_per_contig, ptAlignment **alignments, int alignments_len);

// defined for ptBlock->extend_data
void ptVariant_extend_stList(void *curr_vars_, void *new_vars_);

// defined for ptBlock->copy_data
void *ptVariant_copy_stList(void *vars_);

// defined for ptBlock->destruct_data
void *ptVariant_destruct_stList(void *vars_);

stList *ptVariant_subset_stList(stList *variants, stHash *blocks_per_contig);

#endif /* PT_VARIANT_H */
