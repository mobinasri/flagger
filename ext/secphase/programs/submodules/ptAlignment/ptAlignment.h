#ifndef PT_ALIGNMENT_H
#define PT_ALIGNMENT_H

#include "sam.h"
#include "faidx.h"
#include "sonLib.h"
#include <assert.h>
#include <math.h>
#include <float.h>
#include "cigar_it.h"
#include "common.h"

/*! @typedef
 * @abstract Structure for an alignment record and its marker consistency score
 * @field record	Pointer to a bam1_t struct that contains the alignment record
 * @field score		marker consistency score (summation of the base quality values of the markers)
 * @field conf_blocks	confident blocks (blocks around markers. They do not contain any insertions longer than a threshold. 
 * 			They are useful for faster BAQ calculation)
 * @field flank_blocks	flanking blocks (windows of a specific length around markers. 
 * 			They are useful for faster BAQ calculation)
 * @field rfs		start coordinate on reference (0-based inclusive)
 * @field rfe		end coordinate on reference (0-based inclusive)
 * @field rds_f		start coordinate on read's forward strand (0-based inclusive)
 * @field rde_f		end coordinate on read's forward strand (0-based inclusive)
 */
typedef struct {
    bam1_t* record;
    char contig[200];
    double score;
    stList* conf_blocks;
    stList* flank_blocks;
    int rfs;
    int rfe;
    int rds_f;
    int rde_f;
}ptAlignment;

/**
 * Construct a ptAlignment struct
 *
 * @param record 	Alignment record
 * @param sam_hdr   SAM header
 * @return alignment	Alignment saved as a ptAlignment struct
 *
 */
ptAlignment* ptAlignment_construct(bam1_t* record, sam_hdr_t *sam_hdr);


/**
 * Set start and end coordinates of an alignment
 * It will set the attributes rfs, rfe, rds_f, rde_f
 *
 * @param alignment        Alignment saved as a ptAlignment struct
 *
 */
void ptAlignment_init_coordinates(ptAlignment* alignment);


/**
 * Destruct a ptAlignment struct
 *
 * @param alignment        Alignment saved as a ptAlignment struct 
 *
 */
void ptAlignment_destruct(ptAlignment* alignment);

/**
 * Get the number of supplementary alignments
 */
int ptAlignment_supplementary_count(ptAlignment **alignments, int alignments_len);

/**
 * Get the number of primary alignments
 */
int ptAlignment_primary_count(ptAlignment **alignments, int alignments_len);


/**
 * After setting alignment score values by calling calc_alignment_score
 * this function can be called to get the index of the alignment
 * with the highest score
 *
 * @param alignments            An array of alignments each of which saved as ptAlignment struct
 * @param alignments_len        Length of alignments array
 * @param prim_margin		    The score of the best alignment should be higher than the
 * 				                score of the current primary alignment by this value
 * @param min_score		        The minimum score of the best alignment. If the best alignment
 * 				                had a score lower than this value the primary alignment will be
 * 				                returned as the best alignment
 * @param prim_margin_random	If the scores of the best alignment and current primary alignment
 * 				                was closer than this value then one of them will be returned randomly
 *
 */

int get_best_record_index(ptAlignment** alignments, int alignments_len, double prim_margin, double min_score, double prim_margin_random);


void print_contigs(ptAlignment** alignments, int alignments_len);


int get_primary_index(ptAlignment **alignments, int alignments_len);

#endif /* PT_ALIGNMENT_H */
