#ifndef PT_MARKER_H
#define PT_MARKER_H

#include <getopt.h>
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
#include "ptAlignment.h"
#include "ptBlock.h"


/*! @typedef
 * @abstract Structure for saving a marker's information.
 * @field alignment_idx	    Index of the alignment (could be secondary or primary) from where this marker is taken
 * @field read_pos_f	    Position of the marker in the read (All positions should be computed w.r.t
 * 				            the forward strand for the inconsistency between negative and positive alignments)
 * @field base_idx		    Index of the base in the query seq of the corresponding alignment. It may be different
 * 				            from read_pos_f because of two reasons; negative orientation and hard clipping
 * @field base_q		    Base quality of the marker (could be the raw quality or any processed form e.g BAQ)
 * @field is_match		    True if the marker base is equal to the reference and False otherwise
 */
typedef struct {
	int32_t alignment_idx;
	int32_t read_pos_f;
	int32_t base_idx;
	int32_t base_q;
	bool is_match;
	int32_t ref_pos;
}ptMarker;


/**
 * Create a ptMarker struct
 *
 * @param alignment_idx		The index of the alignment (could be secondary or primary) from where this marker is taken
 * @param read_pos_f		The position of the marker in the read (All positions should be computed w.r.t to
 * 				            the forward strand for the inconsistency between negative and positive alignments)
 * @param base_q    		The base quality of the marker (could be the raw quality or any processed form e.g BAQ)
 * @param is_match		    True if the marker base is equal to the reference and False otherwise
 * @return 			        Pointer to a new marker struct on success; NULL on failure
 * 			                The ptMarker struct returned by a successful call should be freed
 * 				            via free() when it is no longer needed.
 */
ptMarker* ptMarker_construct(int32_t alignment_idx, int32_t base_idx, int32_t read_pos_f, uint8_t base_q, bool is_match, int32_t ref_pos);



/**
 * Make a copy of a marker saved as a ptMarker struct
 *
 * @param src		Marker to be copied
 * @return dest		Copied marker
 *
 */
ptMarker* ptMarker_copy(ptMarker* src);



/**
 * Compare two ptMarker structs.
 * The comparison is first based on their positions in the read and if both given markers
 * are located at the same position then it looks at alignment_idx.
 *
 * @param a	Pointer to the first ptMarker structure
 * @param b	Pointer to the second ptMarker structure
 * @return 	Negative when a should be before b and positive otherwise
 *
 */
int ptMarker_cmp(const void *a, const void *b);


/**
 * Get the initial list of markers that only contains mismatch bases
 *
 * @param alignments            An array of alignments each of which saved as ptAlignment struct
 * @param alignments_len        Length of alignments array
 * @param min_q                 Minimum raw base quality of a marker
 *
 */
stList* ptMarker_get_initial_markers(ptAlignment** alignments, int alignments_len, int min_q);


/**
 * This function will create a marker with its attribute is_match set to true
 * It iterates over the cigar string to find the related seq and ref coordinates too
 * and put this information into the created marker
 *
 * @param alignments        An array of alignments each of which saved as ptAlignment struct
 * @param alignment_idx		Index of the alignment that should be used for making a match marker
 * @param read_pos_f		Read position (not seq position) for which a match marker should be created
 * @return match		    Match marker created based on the given coordinate
 *
 */
ptMarker* ptMarker_construct_match(ptAlignment** alignments, int32_t alignment_idx, int32_t read_pos_f);


/**
 * If there exists a marker with a quality lower than the given threshold
 * this function will remove all the markers in the same read position (read_pos_f)
 * Note that a single read position may have as many markers as up to the number of alignments
 * If BAQ is calculated for read bases a single position may obtain different quality
 * values for different alignments. This function removes markers from the same position
 * if in at least one alignment it obtained a BAQ lower than the threshold.
 *
 * @param markers_p     Pointer to a stList of markers each of which saved as ptMarker struct
 * @param threshold		Minimum quality value for a marker
 *
 */
void filter_lowq_markers(stList** markers_p, int threshold);


/**
 * If there exist a read position which is located in an insertion in at least
 * one alignment this function removes all the related markers in that position
 *
 * @param markers_p             Pointer to a stList of markers each of which saved as ptMarker struct
 * @param alignments            An array of alignments each of which saved as ptAlignment struct
 * @param alignments_len        Length of alignments array
 *
 */
void filter_ins_markers(stList** markers_p, ptAlignment** alignments, int alignments_len);



/**
 * If there exist a read position which is mismatch in all alignments
 * this function will remove those related markers since it is probably
 * a read error.
 *
 * @param markers_p             Pointer to a stList of markers each of which saved as ptMarker struct
 * @param alignments_len        Length of alignments array
 *
 */

void remove_all_mismatch_markers(stList** markers_p, int alignments_len);


/**
 * It is assumed that all markers in the given list are mismatch (is_match == false)
 * A mismatch marker maybe a match in another alignment.
 * This function sorts markers first based on read_pos_f and then alignment index.
 * It adds match markers in a way that after calling this function for each marker
 * position we will have exactly one match or mismatch marker per alignment.
 *
 * @param markers_p       	Pointer to a stList of markers each of which saved as ptMarker struct
 *                      	All of these markers are assumed to be mismatch
 * @param alignments    	An array of alignments each of which saved as ptAlignment struct
 * @param alignments_len	Length of alignments array
 *
 */
void sort_and_fill_markers(stList** markers_p, ptAlignment** alignments, int alignments_len);


/**
 * Quality value (q) is calculated by this formula as a function of probability (p)
 * q(p) = -10 * log(p)
 * This function receives q(p) and returns q(1-p)
 * As edge cases:
 * 1. If q == 0 it reutrns  93
 * 2. If q >= 93 it returns 0
 *
 * @param q		The base quality q(p) (either raw or BAQ)
 * @return rev_q 	The reversed base quality q(1-p)
 */
double reverse_quality(uint8_t q);


//
/**
 * Calculate the alignment scores based on the base qualities
 * of the selected markers (could be raw base quality or BAQ)
 * Each score is in logarithmic scale and represents the probability
 * of all markers being match in the corresponding alignment.
 * We will have two cases:
 *  1. If a marker (base) is match :
 *      The probability of being a real match is equal to 1 - P(error)
 *      Q = -10 * log(P(error)) is saved as base quality so
 *      by calling reverse_quality we can calculate 10 * log(1 - P(error)) (with a negative sign)
 *      is which a quality value of being a real match. It's usually close to 0.0
 *  2. If a marker (base) is mismatch:
 *      The probability of being a real match is equal to P(error) / 3
 *      assuming that it is evenly probable to be a base other than what
 *      is now in the read sequence and one of those bases can be match to the
 *      reference sequence. Therefore the related quality value can be
 *      calculated by 10 * log(P(error)) - 10 * log(3) = -Q - 10 * log(3)
 *      Note that the quality value does not have a negative sign in the formula so it is
 *      always a negative value. The total score of each alignment is the summation of the
 *      quality values of markers. One marker may obtain different quality values for
 *      different alignments both because of different BAQ values and being match or
 *      mismatch in different alignments.An alignment with higher score is more probable to
 *      be to the correct haplotype.
 *
 * @param markers	    stList of markers each of which saved as ptMarker struct
 * @param alignments    An array of alignments each of which saved as ptAlignment struct
 */
void calc_alignment_score(stList* markers, ptAlignment** alignments);


/**
 * Find the confident blocks for the given alignment.
 * Confident blocks are maximal blocks with no insertion longer
 * than the given threshold
 *
 * @param alignment	Alignment record saved as a ptAlignment struct
 * @param threshold	Maximum length of an insertion that may exist in a confident block
 */
stList* find_confident_blocks(ptAlignment* alignment, int threshold);


/**
 * Get the intersection of two arrays of blocks
 *
 * @param blocks1     1st set of blocks
 * @param blocks2     2nd set of blocks
 *
 * @return          stList of intersected blocks
 *
 * @note        Blocks should be sorted by rds_f and no overlap
 * 		        should exist between the blocks of each set
 */
stList* intersect_by_rd_f(stList* blocks1, stList* blocks2);


/**
 * Given a list of alignments for a single read
 * find the confident blocks for each alignment
 * and update the corresponding attribute (alignment->conf_blocks)
 * for each alignment
 *
 * @param alignments		An array of alignments each of which saved as ptAlignment struct
 * @param alignments_len	Length of alignments array
 * @param threshold		    Maximum length of an insertion in each confident block
 */
void set_confident_blocks(ptAlignment** alignments, int alignments_len, int threshold);


/**
 * Given an alignment and the list of markers
 * find the blocks flanking the markers
 * and update the corresponding attribute (alignment->flank_blocks)
 *
 * @param alignments		An array of alignments each of which saved as ptAlignment struct
 * @param alignments_len	Length of alignments array
 * @param margin		    The flanking block around each marker
 * 				            will extend from each side by this margin
 *
 * @note 			        The markers should be sorted by read_pos_f
 */
stList* find_flanking_blocks(ptAlignment* alignment, stList* markers, int margin);


/// Given a list of alignments for a single read
/// find the flanking blocks within each alignment
/// and update the corresponding attribute (alignment->flank_blocks)
/**
 * @param alignments		An array of ptAlignment structure
 * @param alignments_len	Length of alignments array
 * @param margin		    For each marker we will dedicate a flanking block that spans
 * 				            a window of length (2 * margin + 1) with the marker in its center
 */
void set_flanking_blocks(ptAlignment** alignments, int alignments_len, stList* markers, int margin);


/**
 * Given a list of alignments for a single read,
 * find the consensus confident blocks using all previously
 * found confident blocks and flanking blocks. In other words
 * it takes the intersection of all sets of confident and
 * flanking blocks for a single read.
 * Before calling this function there should exist one
 * set of confident blocks and one set of flanking blocks per
 * alignment. They can be set by calling the functions set_flanking_blocks
 * and set_confident_blocks
 *
 * @param alignments		An array of ptAlignment structure
 * @param alignments_len	Length of alignments array
 * @param threshold	        Indels strictly longer than this threshold will cut the confident blocks
 *
 * @return 			        The number of output confident blocks
 *
 * @note 			        The conf_blocks attribute of all the alignments will be updated.
 *                          So all alignments will have the same conf_blocks after
 * 		                    calling this function.
 */
int correct_conf_blocks(ptAlignment** alignments, int alignments_len, int threshold);



/**
 * Given a list of alignments for a single read,
 * Determine if there is any confident block whose length is exceeding the given
 * threshold. This length can be either in the reference or read coordinates
 * This function can be used to shorten the blocks not to be
 * longer than a threshold
 *
 * @param alignments     Array of ptAlignment structures
 * @param alignments_len     Length of alignments array
 * @param threshold     Length threshold
 * @param sam_hdr	SAM file header structure
 */
bool needs_to_find_blocks(ptAlignment** alignments, int alignments_len, int threshold, sam_hdr_t* sam_hdr);


/**
 * Calculate BAQ for all bases in the confident blocks of the given
 * alignment
 *
 * @param fai			Fasta index
 * @param contig_name		Contig name
 * @param alignment 		Alignment saved as a ptAlignment struct
 * @param alignment_idx		Index of the alignment 
 * @param markers		stList of markers
 * @param conf_d		Gap open probability for realignment
 * @param conf_e		Gap extension probability for realignment
 * @param conf_bw		Bandwidth for realignment
 * @param set_q			Base qualities are set to this value prior to realignment
 *
 */
void calc_local_baq(const faidx_t* fai, const char* contig_name, ptAlignment* alignment, int alignment_idx, stList* markers, double conf_d, double conf_e, double conf_bw, int set_q);


void calc_update_baq_all(const faidx_t* fai, 
		        ptAlignment** alignments, int alignments_len, 
			stList* markers, 
			const sam_hdr_t *h,
			double conf_d, double conf_e, double conf_b, int set_q);


void print_markers(stList* markers);

void ptMarker_add_marker_blocks_by_contig(stHash *blocks_per_contig, char* contig, int alignment_idx, stList* markers);

ptBlock* ptMarker_convert_to_block(ptMarker* marker);

#endif /* PT_MARKER_H */
