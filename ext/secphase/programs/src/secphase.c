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


#define ARRAY_SIZE(arr) (sizeof((arr)) / sizeof((arr)[0]))

//#define DEBUG

#ifdef DEBUG
//void DEBUG_PRINT(const char *, ...);
void DEBUG_PRINT(const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
}

#else
static inline void DEBUG_PRINT(const char *fmt, ...) {};
#endif


/*! @typedef
 @abstract Structure for a marker information
 @field alignment_idx	The index of the alignment (could be secondary or primary) from where this marker is taken
 @field read_pos_f	The position of the marker in the read (All positions should be computed w.r.t to
			the forward strand for the inconsistency between negative and positive alignments)
 @field base_q	The base quality of the marker (could be the raw quality or any processed form e.g BAQ)
 @field is_match		True if the marker base is equal to the reference and False otherwise
 */
typedef struct {
	int32_t alignment_idx;
	int32_t read_pos_f;
	int32_t base_idx;
	int32_t base_q;
	bool is_match;
	int32_t ref_pos;
}ptMarker;


/// Create a ptMarker structure
/**
   @param alignment_idx   The index of the alignment (could be secondary or primary) from where this marker is taken
   @param read_pos_f        The position of the marker in the read (All positions should be computed w.r.t to
   			  the forward strand for the inconsistency between negative and positive alignments)
   @param base_q    The base quality of the marker (could be the raw quality or any processed form e.g BAQ)
   @param is_match         True if the marker base is equal to the reference and False otherwise
   @return Pointer to a new marker struct on success; NULL on failure
   The ptMarker struct returned by a successful call should be freed
   via free() when it is no longer needed.
 */
ptMarker* ptMarker_construct(int32_t alignment_idx, int32_t base_idx, int32_t read_pos_f, uint8_t base_q, bool is_match, int32_t ref_pos){
	ptMarker* marker = (ptMarker*) malloc(sizeof(ptMarker));
	marker->alignment_idx = alignment_idx;
        marker->is_match = is_match;
        marker->read_pos_f = read_pos_f;
	marker->base_idx = base_idx;
        marker->base_q = base_q;
	marker->ref_pos = ref_pos;
	return marker;
}

ptMarker* ptMarker_copy(ptMarker* src){
        ptMarker* dest = (ptMarker*) malloc(sizeof(ptMarker));
        dest->alignment_idx = src->alignment_idx;
        dest->is_match = src->is_match;
        dest->read_pos_f = src->read_pos_f;
        dest->base_idx = src->base_idx;
        dest->base_q = src->base_q;
	dest->ref_pos = src->ref_pos;
        return dest;
}

/*! @typedef
 @abstract Structure for an alignment record and its likelihood
 @field record	Pointer to a bam1_t struct that contains the alignment record
 @field p	The likelihood of being the correct alignment
 */
typedef struct {
        bam1_t* record; 
	double qv; // The probability of being the correct alignment in log scale (-10*log10(P))
	stList* conf_blocks; // confident blocks
	stList* flank_blocks; // flanking blocks (blocks that exclude the regions far from markers they are useful for faster BAQ calculation)
	int rfs; // start coordinate on reference (0-based inclusive)
	int rfe; // end coordinate on reference (0-based inclusive)
	int rds_f; // start coordinate on read's forward strand (0-based inclusive)
	int rde_f; // end coordinate on read's forward strand (0-based inclusive)
}ptAlignment;

ptAlignment* ptAlignment_construct(bam1_t* record, double qv){
        ptAlignment* alignment = (ptAlignment*) malloc(sizeof(ptAlignment));
        bam1_t* record_cpy = bam_init1();
        assert(bam_copy1(record_cpy, record) != NULL);
        alignment->record = record_cpy;
        alignment->qv = qv;
        alignment->conf_blocks = NULL;
	alignment->flank_blocks = NULL;
	alignment->rfs = -1;
	alignment->rfe = -1;
	alignment->rds_f = -1;
	alignment->rde_f = -1;
        return alignment;
}

void ptAlignment_destruct(ptAlignment* alignment){
        if (alignment->record){
                bam_destroy1(alignment->record);
                alignment->record = NULL;
        }
        if (alignment->conf_blocks){
                stList_destruct(alignment->conf_blocks);
                alignment->conf_blocks = NULL;
        }
	if (alignment->flank_blocks){
                stList_destruct(alignment->flank_blocks);
                alignment->flank_blocks = NULL;
        }
        free(alignment);
}


typedef struct {
        int rfs; // ref start
	int rfe; // ref end
	int sqs; // seq start
	int sqe; // seq end
	int rds_f;
	int rde_f;
}ptBlock;

ptBlock* ptBlock_construct(int rfs, int rfe, int sqs, int sqe, int rds_f, int rde_f){
	ptBlock* block = malloc(sizeof(ptBlock));
	block->rfs = rfs;
	block->rfe = rfe;
	block->sqs = sqs;
	block->sqe = sqe;
	block->rds_f = rds_f;
	block->rde_f = rde_f;
	return block;
}

int ptBlock_cmp_rds_f(const void *a, const void *b){
        ptBlock* b1 = (ptBlock*) a;
        ptBlock* b2 = (ptBlock*) b;
        return b1->rds_f - b2->rds_f;
}

int ptBlock_cmp_sqs(const void *a, const void *b){
        ptBlock* b1 = (ptBlock*) a;
        ptBlock* b2 = (ptBlock*) b;
        return b1->sqs - b2->sqs;
}


/// Compare two ptMarker structures
/**
   @param a	Pointer to the first ptMarker structure
   @param b	Pointer to the second ptMarker structure
   @return 	Nagative when a should be before b and positive otherwise
   @note	The comparison is first based on their positions in the read and if both given markers
		are located at the same position then we look at the alignment_idx 
 */
int ptMarker_cmp(const void *a, const void *b){
        ptMarker* m1 = (ptMarker*) a;
        ptMarker* m2 = (ptMarker*) b;
	if (m1->read_pos_f == m2->read_pos_f) {
		return m1->alignment_idx - m2->alignment_idx;
	}
	else {
		return m1->read_pos_f - m2->read_pos_f;
	}
}


ptMarker* ptMarker_construct_match(ptAlignment** alignments, int32_t alignment_idx, int32_t read_pos_f){
	bam1_t* record = alignments[alignment_idx]->record;
        uint8_t* q = bam_get_qual(record);
	ptMarker* match;
	uint32_t* cigar = bam_get_cigar(record);
	int lclip=0;
	int rclip=0;
	if (bam_cigar_op(cigar[0]) == BAM_CHARD_CLIP) {
                lclip = bam_cigar_oplen(cigar[0]);
        }
        if (bam_cigar_op(cigar[record->core.n_cigar - 1]) == BAM_CHARD_CLIP) {
                rclip = bam_cigar_oplen(cigar[record->core.n_cigar - 1]);
        }
	// construct the match marker
	if (bam_is_rev(record)){
		match = ptMarker_construct(alignment_idx,
					   record->core.l_qseq + rclip - read_pos_f - 1,
                                           read_pos_f,
                                           q[record->core.l_qseq + rclip - read_pos_f - 1],
                                           true,
					   -1);
        }
        else {
		match = ptMarker_construct(alignment_idx,
				           read_pos_f - lclip,
                                           read_pos_f,
                                           q[read_pos_f - lclip],
                                           true,
					   -1);
        }
	return match;
}


void filter_lowq_markers(stList** markers_p, int threshold){
	stList* markers = *markers_p;
	int markers_len = stList_length(markers);
	if(markers_len == 0) return;
        stList* keep_markers = stList_construct3(0, free);
	ptMarker* pre_marker =NULL;
	ptMarker* marker;
	bool keep_flag = false;
	int idx_s = 0;
	int min_q = 100;
	ptMarker* marker_copy;
	for(int j=0; j < markers_len; j++){
		marker = stList_get(markers, j);
		if (pre_marker && marker->read_pos_f != pre_marker->read_pos_f){
			if (min_q > threshold){
				for(int k=idx_s; k<j; k++){
					// copy all markers with the same read_pos_f and add them
					marker_copy = ptMarker_copy(stList_get(markers, k));
					marker_copy->base_q = min_q;
					stList_append(keep_markers, marker_copy);
				}
			}
			idx_s = j;
			//keep_flag = false;
			min_q = 100;
		}
		if (min_q > marker->base_q) min_q = marker->base_q;
		/*
		if(marker->base_q > threshold) {
			keep_flag |= true;
			DEBUG_PRINT("q = %d ( > %d) marker->read_pos_f = %d\n", marker->base_q, threshold, marker->read_pos_f);
		}*/
		pre_marker = marker;
	}
	if (min_q > threshold){
		for(int k=idx_s; k<markers_len; k++){
			marker_copy = ptMarker_copy(stList_get(markers, k));
                        marker_copy->base_q = min_q;
                        stList_append(keep_markers, marker_copy);
                }
        }
	stList_destruct(markers);
	*markers_p = keep_markers;
}

void filter_ins_markers(stList** markers_p, ptAlignment** alignments, int alignments_len){
	stList* markers = *markers_p;
	int markers_len = stList_length(markers);
	if(markers_len == 0) return;
	bool* markers_keep_flags = (bool*) malloc(markers_len * sizeof(bool));
	for(int i=0; i< markers_len; i++) markers_keep_flags[i] =true;
	bam1_t* b;
	ptMarker* marker;
	ptCigarIt* cigar_it;
	int j; // marker index
        for(int i=0; i < alignments_len; i++){
		b = alignments[i]->record;
		j = bam_is_rev(b) ? stList_length(markers) - 1 : 0;
		marker = stList_get(markers, j);
		cigar_it = ptCigarIt_construct(b, true, true);
		//iterate over cigars
		while(ptCigarIt_next(cigar_it)){
			//check if marker is located within the cigar operation
			while(marker &&
			      (cigar_it->rds_f <= marker->read_pos_f) && 
			      (cigar_it->rde_f >= marker->read_pos_f)){
				//check if the marker is located within an insertion
				if(cigar_it->op == BAM_CINS ||
				   cigar_it->op == BAM_CSOFT_CLIP ||
				   cigar_it->op == BAM_CHARD_CLIP){
					markers_keep_flags[j] = false;
				}
				//update reference positions in the meanwhile!
				if(cigar_it->op == BAM_CEQUAL && marker->alignment_idx == i){
                                	marker->ref_pos = bam_is_rev(b) ? cigar_it->rfs + cigar_it->rde_f - marker->read_pos_f : cigar_it->rfs + marker->read_pos_f - cigar_it->rds_f;
                                }
				j += bam_is_rev(b) ? -1 : 1;
				marker = ((j < markers_len) && (j >= 0)) ? stList_get(markers, j) : NULL;
			}	
		}//end of iteration over cigar ops
		ptCigarIt_destruct(cigar_it);
	}// end of iteration over alignments
	// make a list for saving the remaining markers
	stList* keep_markers = stList_construct3(0, free);
	for(j=0; j < markers_len; j++){
		if(markers_keep_flags[j]){
			marker = stList_get(markers, j);
			stList_append(keep_markers, ptMarker_copy(marker));
		}
	}
	stList_destruct(markers);
        *markers_p = keep_markers;
}


//should be run before matches are added
void remove_all_mismatch_markers(stList** markers_p, int alignments_len){
	stList* markers = *markers_p;
	stList_sort(markers, ptMarker_cmp);
	int markers_len = stList_length(markers);
	if(markers_len == 0) return;
        bool* markers_keep_flags = (bool*) malloc(markers_len * sizeof(bool));
        memset(markers_keep_flags, 1, markers_len);
	ptMarker* marker;
	ptMarker* pre_marker = NULL;
	int occ = 0;
	for(int i=0; i < markers_len; i++){
		marker = stList_get(markers, i);
		if (pre_marker  &&
                    pre_marker->read_pos_f < marker->read_pos_f){
			if(occ == alignments_len){// if all markers are mismatch
				memset(markers_keep_flags + i - alignments_len, 0, alignments_len);
			}
			occ = 0;
		}
		pre_marker = marker;
		occ += 1;
	}
	if(occ == alignments_len){
		memset(markers_keep_flags + markers_len - alignments_len, 0, alignments_len);
        }
	// make a list for saving the remaining markers
        stList* keep_markers = stList_construct3(0, free);
        for(int j=0; j < markers_len; j++){
                if(markers_keep_flags[j]){
                        marker = stList_get(markers, j);
                        stList_append(keep_markers, ptMarker_copy(marker));
                }
        }
	stList_destruct(markers);
	*markers_p = keep_markers;
}


void sort_and_fill_markers(stList** markers_p, ptAlignment** alignments, int alignments_len){
	stList* markers = *markers_p;
	stList_sort(markers, ptMarker_cmp);
	int idx=0;
	ptMarker* match;
	ptMarker* marker;
	ptMarker* pre_marker = NULL;
	int markers_len = stList_length(markers);
	for(int64_t i=0; i < markers_len; i++){
		marker = stList_get(markers, i);
		assert(!marker->is_match);
		if (pre_marker  &&
		    pre_marker->read_pos_f < marker->read_pos_f){
			for (int64_t j=idx; j < alignments_len; j++){
                        	match = ptMarker_construct_match(alignments, j, pre_marker->read_pos_f);
                        	stList_append(markers, match);
                	}
			idx = 0;
		}
		pre_marker = marker;
		for (int64_t j=idx; j < marker->alignment_idx; j++){
			match = ptMarker_construct_match(alignments, j, marker->read_pos_f);
                	stList_append(markers, match);
        	}
		idx = marker->alignment_idx + 1;
	}
	if(markers_len > 0){
		// add the last matches if there still exist to be added
		for (int64_t j=idx; j < alignments_len; j++){
			match = ptMarker_construct_match(alignments, j, marker->read_pos_f);
        		stList_append(markers, match);
        	}
	}
	stList_sort(markers, ptMarker_cmp); // sort again to locate the matches in their correct positions
}

double revert_quality(uint8_t q){
	if (q >= 93) return 0;
	if (q == 0) return 93;
	double p = 1 - pow(10, (double)q / -10);
	double rev_q = -10 * log(p);
	return rev_q;
}

void calc_likelihood(stList* markers, ptAlignment** alignments){
        ptMarker* marker;
        for(int64_t j=0; j < stList_length(markers); j++){
                marker =  stList_get(markers, j);
		if (marker->is_match){
			alignments[marker->alignment_idx]->qv -= revert_quality(marker->base_q);
		}
		else{
			alignments[marker->alignment_idx]->qv -= marker->base_q;
		}
	}
}

int get_best_record(ptAlignment** alignments, int alignments_len, double prim_margin, double min_qv, double prim_margin_random){
	assert(alignments_len > 0);
	if (alignments_len == 1) return 0;
	double max_qv = -DBL_MAX;
	int max_idx = -1;
	double prim_qv = -DBL_MAX;
	int prim_idx = -1;
	for(int i=0; i < alignments_len; i++){
		if ((alignments[i]->record->core.flag & BAM_FSECONDARY) == 0){
			prim_idx = i;
			prim_qv = alignments[i]->qv;
		}
		else if (max_qv < alignments[i]->qv){
			max_idx = i;
			max_qv = alignments[i]->qv;
		}
	}
	int rnd = rand() % 2; // 50% chance for rnd=0 (same for rnd=1)
	if (abs(max_qv - prim_qv) < prim_margin_random){
		return rnd == 0 ? prim_idx : max_idx; // if the scores were closer than prim_margin_random return one randomly
	}
	if (prim_idx == -1 || max_qv < (prim_qv + prim_margin) || max_qv < min_qv) return prim_idx;
	else return max_idx;
}

stList* find_confident_blocks(ptAlignment* alignment, int threshold){
	bam1_t* b = alignment->record;
	stList* conf_blocks = stList_construct3(0, free);
	ptCigarIt* cigar_it = ptCigarIt_construct(b, true, true);
	int conf_sqs = 0; // the seq start of the confident block
	int conf_rfs = b->core.pos; // the ref start of the confident block
	int conf_rd_f = bam_is_rev(b) ? cigar_it->rde_f : cigar_it->rds_f;
	ptBlock* block = NULL;
	while(ptCigarIt_next(cigar_it)){
		switch(cigar_it->op) {
                	case BAM_CINS:
			case BAM_CDEL:
                        	if(cigar_it->len > threshold && 
				   conf_sqs < cigar_it->sqs &&
				   conf_rfs < cigar_it->rfs){
					if (bam_is_rev(b)) {
						block = ptBlock_construct(conf_rfs, cigar_it->rfs - 1,
								          conf_sqs, cigar_it->sqs - 1,
									  cigar_it->rde_f + 1, conf_rd_f);
					}
					else {
						block = ptBlock_construct(conf_rfs, cigar_it->rfs - 1,
								          conf_sqs, cigar_it->sqs - 1,
									  conf_rd_f, cigar_it->rds_f - 1);
					}
					stList_append(conf_blocks, block);
				}
				if(cigar_it->len > threshold){
					conf_sqs = cigar_it->sqe + 1;
                                	conf_rfs = cigar_it->rfe + 1;
					conf_rd_f = bam_is_rev(b) ? cigar_it->rds_f - 1 : cigar_it->rde_f + 1;
				}
                        	break;
                	case BAM_CSOFT_CLIP:
			case BAM_CHARD_CLIP:
                        	if(conf_sqs < cigar_it->sqs &&
				   conf_rfs < cigar_it->rfs){
					if (bam_is_rev(b)) {
                                        	block = ptBlock_construct(conf_rfs, cigar_it->rfs - 1,
                                                                  	  conf_sqs, cigar_it->sqs - 1,
									  cigar_it->rde_f + 1, conf_rd_f);
					}
					else {
						block = ptBlock_construct(conf_rfs, cigar_it->rfs - 1,
                                                                          conf_sqs, cigar_it->sqs - 1,
                                                                          conf_rd_f, cigar_it->rds_f - 1);
					}
					stList_append(conf_blocks, block);
				}
                                conf_sqs = cigar_it->sqe + 1;
                                conf_rfs = cigar_it->rfe + 1;
				conf_rd_f = bam_is_rev(b) ? cigar_it->rds_f - 1 : cigar_it->rde_f + 1;
                        	break;
        	}

	}
	// when the last block is not terminated with SOFT or HARD clipping
	if(conf_sqs <= cigar_it->sqe){
		if (bam_is_rev(b)) {
			block = ptBlock_construct(conf_rfs, cigar_it->rfe,
                                                  conf_sqs, cigar_it->sqe,
						  cigar_it->rds_f, conf_rd_f);
		}
		else {
			block = ptBlock_construct(conf_rfs, cigar_it->rfe,
                                                  conf_sqs, cigar_it->sqe,
						  conf_rd_f, cigar_it->rde_f);
		}
		stList_append(conf_blocks, block);
	}
	ptCigarIt_destruct(cigar_it);
	return conf_blocks;
}

/// Get the intersection of the two given blocks
/**
   @param blocks1     1st set of blocks
   @param blocks2     2nd set of blocks
   @return      A stList of intersected blocks
   @note        blocks should be sorted by rds_f and no overlaps should exist between the blocks within each set
 */
stList* intersect(stList* blocks1, stList* blocks2){
	stList* blocks_intersect = stList_construct3(0, free);
	if (stList_length(blocks1) == 0 || stList_length(blocks2) == 0) return blocks_intersect;
	ptBlock* b1;
	ptBlock* b2;
	ptBlock* b;
	int j = 0; 
	b2 = stList_get(blocks2, j);
	for (int i = 0; i < stList_length(blocks1); i++){
		b1 = stList_get(blocks1, i);
		// Find the first block (b2) in the 2nd set which end point is after 
		// the start point of the current block (b1) in the 1st set
                while (b2 && b2->rde_f < b1->rds_f) {
			j++;
                        b2 = j == stList_length(blocks2) ? NULL : stList_get(blocks2, j);
                }
		// Because of the previous loop, b2's end point is now after b1's
		// start point. If b2's start point is also before the b1's end point
		// it means there exist an overlap so we save that overlap in 
		// blocks_intersect and go to the next block (update b2)
                while (b2 && b2->rds_f < b1->rde_f){
			b = ptBlock_construct(-1, 
					      -1, 
					      -1, 
					      -1, 
					      max(b1->rds_f, b2->rds_f), 
					      min(b1->rde_f, b2->rde_f));
                        stList_append(blocks_intersect, b);
                        if (b2->rde_f <= b1->rde_f) {
                                j++;
                        	b2 = j == stList_length(blocks2) ? NULL : stList_get(blocks2, j);
                        }
                        else break;
                }
        }
	return blocks_intersect;
}


/// Given a list of alignments for a single read 
/// find the confident blocks within each alignment 
/// and update the corresponding attribute (alignment->conf_blocks)
/**
   @param alignments     An array of ptAlignment structure
   @param alignments_len     Length of alignments array
   @param threshold	The indels strictly longer than this manually determined 
   			threshold will cut the confident blocks
 */
void set_confident_blocks(ptAlignment** alignments, int alignments_len, int threshold){
	for(int i=0; i < alignments_len; i++){
		alignments[i]->conf_blocks = find_confident_blocks(alignments[i], threshold);
	}
}

// @note: the markers should be sorted by read_pos_f
stList* find_flanking_blocks(ptAlignment* alignment, stList* markers, int margin){
	bam1_t* b = alignment->record;
        stList* flank_blocks = stList_construct3(0, free);
        ptCigarIt* cigar_it = ptCigarIt_construct(b, true, true);
	// get the first merker for initializing start and end
	int i = 0;
	ptMarker* marker = stList_get(markers, i);
	int start = max(alignment->rds_f, marker->read_pos_f - margin);
	int end = min(alignment->rde_f, marker->read_pos_f + margin);
	i++;
	int curr_start;
	int curr_end;
	ptBlock* block = NULL;
        while(i < stList_length(markers)){
		marker = stList_get(markers, i);
		curr_start = max(alignment->rds_f, marker->read_pos_f - margin);
		curr_end =  min(alignment->rde_f, marker->read_pos_f + margin);
		if(curr_start < end) {
			end = curr_end;
		}
		else {
			// add the previous flanking block
			block = ptBlock_construct(-1, -1,
                                                  -1, -1,
                                                  start, end);
                        stList_append(flank_blocks, block);
			//update start and end for the next flanking block
			start = curr_start;
			end = curr_end;
		}
		i++;
	}
	// add the last flanking block
	block = ptBlock_construct(-1, -1,
                                  -1, -1,
                                  start, end);
        stList_append(flank_blocks, block);
	ptCigarIt_destruct(cigar_it);
	return flank_blocks;
}

/// Given a list of alignments for a single read
/// find the flanking blocks within each alignment
/// and update the corresponding attribute (alignment->flank_blocks)
/**
   @param alignments     An array of ptAlignment structure
   @param alignments_len     Length of alignments array
   @param margin     For each marker we will append a flanking block that spans
   		     a window of length 2 * margin with the marker in its center
 */
void set_flanking_blocks(ptAlignment** alignments, int alignments_len, stList* markers, int margin){
        for(int i=0; i < alignments_len; i++){
                alignments[i]->flank_blocks = find_flanking_blocks(alignments[i], markers, margin);
        }
}


/// Given a list of alignments for a single read, 
/// find the consensus confident blocks using all previously 
/// found confident blocks. In other words it takes the 
/// intersection of all sets of confident blocks (one set per alignment).
/**
   @param alignments     An array of ptAlignment structure
   @param alignments_len     Length of alignments array
   @param threshold     The indels strictly longer than this manually determined
                        threshold will cut the confident blocks
 */
int correct_conf_blocks(ptAlignment** alignments, int alignments_len, int threshold){
	assert(alignments_len > 0);
	stList_sort(alignments[0]->conf_blocks, ptBlock_cmp_rds_f);
        stList* blocks = stList_copy(alignments[0]->conf_blocks, NULL);
	ptBlock* block;
        stList* blocks_new;
	//intersect confident blocks
        for(int i=1; i < alignments_len; i++){
		stList_sort(alignments[i]->conf_blocks, ptBlock_cmp_rds_f);
                blocks_new = intersect(blocks, alignments[i]->conf_blocks);
                stList_destruct(blocks);
                blocks = blocks_new;
        }
	//intersect flanking blocks
	for(int i=0; i < alignments_len; i++){
                stList_sort(alignments[i]->flank_blocks, ptBlock_cmp_rds_f);
                blocks_new = intersect(blocks, alignments[i]->flank_blocks);
                stList_destruct(blocks);
                blocks = blocks_new;
        }
	if(stList_length(blocks) == 0){
		DEBUG_PRINT("\t\t\t### No Consensus Blocks!\n");
		for(int i=0; i < alignments_len; i++){
			stList_destruct(alignments[i]->conf_blocks);
			alignments[i]->conf_blocks = stList_construct3(0, free);
		}
		return 0;
	}
	stList* corrected_conf_blocks;
	// For each alignment project the consensus confident blocks to the ref and seq coordinates
	// The main function of the loop below is to fill the attributes rfs, rfe, sqs, and sqe for each block
	// based on the previously calculated rds_f and rde_f
        for(int i=0; i < alignments_len; i++){
		bam1_t* b = alignments[i]->record;
		int j = bam_is_rev(b) ? stList_length(blocks) - 1 : 0;
		corrected_conf_blocks = stList_construct3(0, free);
        	ptCigarIt* cigar_it = ptCigarIt_construct(alignments[i]->record, true, true);
        	int block_rds_f, block_rde_f, cigar_rds_f, cigar_rde_f, rfs, rfe, sqs, sqe;
		bool del_flag = false;
        	block = stList_get(blocks, j);
		// reverse each block interval to make it work with the projecting algorithm below
		if (bam_is_rev(b)){
			block_rds_f = -1 * block->rde_f;
                        block_rde_f = -1 * block->rds_f;
                }
                else {
                        block_rds_f = block->rds_f;
                	block_rde_f = block->rde_f;
               	}
		//iterate over cigar operations from left to right (w.r.t reference)
        	while(ptCigarIt_next(cigar_it)){
			// reverse each cigar operation interval to make it work with the projecting algorithm below
			if (bam_is_rev(b)){
				cigar_rds_f = -1 * cigar_it->rde_f;
				cigar_rde_f = -1 * cigar_it->rds_f;
			}
			else {
				cigar_rds_f = cigar_it->rds_f;
                                cigar_rde_f = cigar_it->rde_f;
			}
			// match/ mismatch/ insertion
			// for M and I the main algorithm is the same except
			// some minor parts that are corrected by 
			// the conditional statement "cigar_it->op == BAM_CINS ? : "
			if(cigar_it->op == BAM_CMATCH ||
			   cigar_it->op == BAM_CEQUAL ||
			   cigar_it->op == BAM_CDIFF ||
			   cigar_it->op == BAM_CINS) {
				//The loop below iterate the blocks overlap with the current cigar operation
				//There may exist multiple blocks within the current op
				while(block && block_rde_f <= cigar_rde_f){
					// update start locations of the projected coordinates
					if(cigar_rds_f <= block_rds_f && !(del_flag && cigar_rds_f == block_rds_f)){
						rfs = cigar_it->op == BAM_CINS ? cigar_it->rfs : cigar_it->rfs + (block_rds_f - cigar_rds_f);
						sqs = cigar_it->sqs + (block_rds_f - cigar_rds_f);
					}
					// update end locations of the projected coordinates
					rfe = cigar_it->op == BAM_CINS ? cigar_it->rfe : cigar_it->rfs + (block_rde_f - cigar_rds_f);
					sqe = cigar_it->sqs + (block_rde_f - cigar_rds_f);
					// add block
					stList_append(corrected_conf_blocks, ptBlock_construct(rfs, 
									                       rfe, 
									                       sqs, 
											       sqe,
											       block->rds_f,
											       block->rde_f));
					// update block index; j
					if (bam_is_rev(b) && j > 0){
                                		j--;
                                                block = stList_get(blocks, j);
                                                block_rds_f = -1 * block->rde_f;
                                                block_rde_f = -1 * block->rds_f;
                        		}
					else if (!bam_is_rev(b) && j < stList_length(blocks) - 1) {
						j++;
                                		block = stList_get(blocks, j);
                                		block_rds_f = block->rds_f;
                                		block_rde_f = block->rde_f;
                        		}
					else if (j == 0 || j == stList_length(blocks) - 1){
						block = NULL;
					}
				}
				if(block == NULL) break;// no more block remaining so break the cigar iteration
				// !(del_flag && cigar_rds_f == block_rds_f) means that,
				// del_flag and cigar_rds_f == block_rds_f cannot be true at the same time
				// if both were true it means that there is a deletion (shorter than threshold)
				// on the edge of the block (This is a rare event can happen when there is a deletion
				// right before an insertion)
				if(cigar_rds_f <= block_rds_f && 
				   block_rds_f <= cigar_rde_f && 
				   !(del_flag && cigar_rds_f == block_rds_f)) {
					rfs = cigar_it->op == BAM_CINS ? cigar_it->rfs : cigar_it->rfs + (block_rds_f - cigar_rds_f);
                                        sqs = cigar_it->sqs + (block_rds_f - cigar_rds_f);
				}
				del_flag = false;
			}
			//deletion
			else if (cigar_it->op == BAM_CDEL){
				int prev_j = bam_is_rev(b) ? j + 1 : j - 1;
				ptBlock* prev_block =  (0 <= prev_j && prev_j <= stList_length(blocks) - 1) ? stList_get(blocks, prev_j) : NULL;
				// if this is a short deletion where the previous block ended right before it
				// ref end coordinate has to be corrected 
				if (prev_block && cigar_it->len <= threshold){
					if ((bam_is_rev(b) && prev_block->rds_f == cigar_it->rds_f) ||
					    (!bam_is_rev(b) && prev_block->rde_f == cigar_it->rde_f)){
						prev_block->rfe = cigar_it->rfe;
					}
				}
				// if this is a short deletion where the current block starts right after it
				// ref start coordinate should be updated
				if (block && block_rds_f == cigar_rds_f && cigar_it->len <= threshold){
                                        del_flag = true;
                                        rfs = cigar_it->rfs;
					sqs = cigar_it->sqs;
                                }
			}
			// assume that the blocks do not overlap with BAM_CSOFT and BAM_CHARD 
		}
		ptCigarIt_destruct(cigar_it);
		// delete previous blocks
		stList_destruct(alignments[i]->conf_blocks);
		//sort by seq (or ref) start coordinates
		stList_sort(corrected_conf_blocks, ptBlock_cmp_sqs);
		// update confident blocks for each alignment object
		alignments[i]->conf_blocks = corrected_conf_blocks;
	}
	return stList_length(blocks);
}


bool needs_to_find_blocks(ptAlignment** alignments, int alignments_len, int threshold, sam_hdr_t* sam_hdr){
	stList* blocks;
	ptBlock* block;
	bool flag = false;
	int tid; const char* contig_name;
	for(int j=0; j < alignments_len; j++){
		blocks = alignments[j]->conf_blocks;
		tid = alignments[j]->record->core.tid;
                contig_name = sam_hdr_tid2name(sam_hdr, tid);
		if (blocks == NULL) return true;
		if (stList_length(blocks) == 0) return true;
		for(int i=0; i < stList_length(blocks); i++){
                	block = stList_get(blocks, i);
                	DEBUG_PRINT("\t\t\t###     block#%d: seq[%d:%d] ref[%s:%d-%d]\n", i, block->sqs, block->sqe, contig_name, block->rfs, block->rfe);
			if ((block->sqe - block->sqs) > threshold || (block->rfe - block->rfs)>threshold) flag = true;
		}
	}
	return flag;
}

void calc_local_baq(const faidx_t* fai, const char* contig_name, ptAlignment* alignment, int alignment_idx, stList* markers, double conf_d, double conf_e, double conf_bw, int set_q){

	uint8_t* tseq; // translated seq A=>0,C=>1,G=>2,T=>3,other=>4
	uint8_t* tref; // translated ref
	uint8_t* bq;
	uint8_t* block_qual;
	int* state;
	uint8_t* q; // Probability of incorrect alignment from probaln_glocal()
	probaln_par_t conf = { conf_d, conf_e, conf_bw };
	char* ref;
	char* reg;
	bam1_t* b = alignment->record;
	uint8_t* seq = bam_get_seq(b);
	uint8_t* qual = bam_get_qual(b);
	stList* blocks = alignment->conf_blocks;
	int markers_len = stList_length(markers);
	int j = bam_is_rev(b) ? stList_length(markers) - 1 : 0;
	ptMarker* marker = stList_get(markers, j);
	ptCigarIt* cigar_it = ptCigarIt_construct(b, true, true);
	ptBlock* block;
	int ref_len; int seq_len;
	DEBUG_PRINT("\t\t\t## Number of Blocks: %ld\n", stList_length(blocks));
	int block_margin = 10;
	for(int i=0; i < stList_length(blocks); i++){
		block = stList_get(blocks, i);
		//DEBUG_PRINT("\t\t\t###     block#%d: seq[%d:%d] ref[%s:%d-%d]\n", i, block->sqs, block->sqe, contig_name, block->rfs, block->rfe);
		while((cigar_it->sqe < block->sqs) || (cigar_it->rfe < block->rfs)){
			if(ptCigarIt_next(cigar_it) == 0) break;
		}// end of while cigar_it->seq_start should be now less than or equal to block->seq_start
		//assert(cigar_it->sqs <= block->sqs);
		//assert(cigar_it->rfs <= block->rfs);
		//printf("$%d-%d : %d\n", block->sqs, block->sqe, marker->base_idx);
		while(marker &&
                      ((marker->base_idx < block->sqs + block_margin) ||
		        (marker->alignment_idx != alignment_idx))){
			j += bam_is_rev(alignment->record) ? -1 : 1;
                        marker = ((j < markers_len) && (j >= 0)) ? stList_get(markers, j) : NULL;
		}// end of while the marker is now the first one located within the current block (or after)
		/*if(marker == NULL) printf("NULL\n");
		else printf("%d-%d : %d\n", block->sqs, block->sqe, marker->base_idx);*/
		if(marker &&
                   (block->sqe >= marker->base_idx + block_margin) &&
		   (block->sqs <= marker->base_idx - block_margin)){
			seq_len = block->sqe - block->sqs + 1;
			ref_len = block->rfe - block->rfs + 1;
			DEBUG_PRINT("\t\t\t###     block#%d: seq[%d:%d] ref[%s:%d-%d]\n", i, block->sqs, block->sqe, contig_name, block->rfs, block->rfe);
			DEBUG_PRINT("\t\t\t###     seq_len:%d\tref_len:%d\n", seq_len, ref_len);
			if (seq_len < 5 || ref_len < 5) continue; // skip small blocks
			//allocate read sequence
			tseq = (uint8_t*) malloc(seq_len);

			for(int k=0; k < seq_len; k++)
				tseq[k] = seq_nt16_int[bam_seqi(seq, block->sqs + k)];
			//allocate reference sequenc
			tref = (uint8_t*) malloc(ref_len);
			reg = malloc(200);
			memset(reg,'\0',200);
			sprintf(reg, "{%s}:%d-%d", contig_name, block->rfs + 1 , block->rfe + 1);
			int len;
			ref = fai_fetch(fai, reg, &len);
			assert(len == (block->rfe - block->rfs + 1));
			for(int k=0; k < ref_len; k++)
                                tref[k] = seq_nt16_int[seq_nt16_table[(unsigned char)ref[k]]];
			//allocate the quality of this confident block
			block_qual = (uint8_t*) malloc(seq_len);
                        for(int k=0; k < seq_len; k++)
				block_qual[k] = set_q;
			//allocate neccessary arrays for HMM BAQ
			state = (int*) malloc((block->sqe - block->sqs + 1) * sizeof(int));
			q = (uint8_t*) malloc(block->sqe - block->sqs + 1);
			//DEBUG_PRINT("Starting local BAQ :))\n");
			conf.bw = abs(ref_len - seq_len) + conf_bw;
			if (probaln_glocal(tref, ref_len, 
					   tseq, seq_len,
					   block_qual, &conf, state, q) == INT_MIN) {
            			fprintf(stderr, "probaln_glocal ERROR\n");
				fprintf(stderr, "%s:%d-%d", contig_name, block->rfs + 1 , block->rfe + 1);
			}
			/*printf("##%d\n",alignment_idx);
			for(int k=0; k < seq_len; k++)
                                printf("%d\t%d\t%c\n", block->sqs + k, q[k], q[k] < 20 ? '*' : ' ');*/
			//DEBUG_PRINT("local BAQ Finished:))\n");
			// the state and q are now updated if there is any marker located within the block
			//apply BAQ to the quality array
			bq = (uint8_t*) malloc(seq_len);
			memcpy(bq, block_qual, seq_len);
			int x;
			int y;
			while(cigar_it->sqs <= block->sqe || cigar_it->rfs <= block->rfe){
				x = cigar_it->rfs - block->rfs;
				x = x < 0 ? 0 : x;
				y = cigar_it->sqs - block->sqs;
				y = y < 0 ? 0 : y;
				//DEBUG_PRINT("rf:%d\t%d\n", cigar_it->rfe, block->rfe);
				//DEBUG_PRINT("sq:%d\t%d\n", cigar_it->sqe, block->sqe);
				if (cigar_it->op == BAM_CMATCH || 
			    	    cigar_it->op == BAM_CEQUAL || 
			            cigar_it->op == BAM_CDIFF) {
					int len = min(cigar_it->len, min(cigar_it->sqe, block->sqe) - max(cigar_it->sqs, block->sqs) + 1);
					//DEBUG_PRINT("\t\t\t#len:%d\tstart:%d\tend:%d\n", len, y, y +len -1);
					for (int t = y; t < (y + len); t++) {
						assert(t < seq_len);
                        			if (((state[t]&3) != 0) || (state[t]>>2 != x + (t - y))) bq[t] = 0;
						else bq[t] = qual[block->sqs + t] < q[t] ? qual[block->sqs + t] : q[t];
                    			}
				}
				if(cigar_it->sqe <= block->sqe || cigar_it->rfe <= block->rfe){
					if(ptCigarIt_next(cigar_it) == 0) break;
				}
				else break;
			}
			for(int k=0; k < seq_len; k++) qual[block->sqs + k] = bq[k] < 94 ? bq[k] : 93;
			free(tseq);
			free(tref);
			free(state);
			free(q);
			free(bq);
			free(block_qual);
			free(ref);
			free(reg);
		}
		else { // for the markers close to the borders of the confident blocks
			while (marker &&
                      	       marker->base_idx <= block->sqe){
				if (marker->alignment_idx == alignment_idx) qual[marker->base_idx] = 0;
                        	j += bam_is_rev(alignment->record) ? -1 : 1;
                        	marker = ((j < markers_len) && (j >= 0)) ? stList_get(markers, j) : NULL;
                	}
		}
	}
	ptCigarIt_destruct(cigar_it);
}

void calc_update_baq_all(const faidx_t* fai, 
		        ptAlignment** alignments, int alignments_len, 
			stList* markers, 
			const sam_hdr_t *h,
			double conf_d, double conf_e, double conf_b, int set_q){
	const char* contig_name;
	int tid;
	for(int i=0; i < alignments_len; i++){
		tid = alignments[i]->record->core.tid;
		contig_name = sam_hdr_tid2name(h, tid);
		calc_local_baq(fai, contig_name, alignments[i], i, markers, conf_d, conf_e, conf_b, set_q);
	}
	ptMarker* marker;
	uint8_t* q;
	//DEBUG_PRINT("\t\t## start updating marker base qualities\n");
	for(int i=0; i < stList_length(markers); i++){
		marker = stList_get(markers, i);
		q = bam_get_qual(alignments[marker->alignment_idx]->record);
		marker->base_q = q[marker->base_idx];
	}
}

bool contain_supp(ptAlignment** alignments, int alignments_len){
	for(int i=0;i < alignments_len; i++){
		if (alignments[i]->record->core.flag & BAM_FSUPPLEMENTARY) return 1;
	}
	return 0;
}

void print_markers(stList* markers){
	if(stList_length(markers) == 0) DEBUG_PRINT("NO MARKERS!\n");
	DEBUG_PRINT("#alignment_idx\tseq_pos\tread_pos_f\tq\tis_match\tref_pos\n");
	for(int t =0 ; t < stList_length(markers); t++) {
        	ptMarker* m = stList_get(markers, t);
                DEBUG_PRINT("%d\t%d\t%d\t%d\t%d\t%d\n", m->alignment_idx, m->base_idx, m->read_pos_f, m->base_q, m->is_match, m->ref_pos);
        }
}

void print_contigs(ptAlignment** alignments, int alignments_len, sam_hdr_t* h){
        DEBUG_PRINT("\nalignments:\n");
	DEBUG_PRINT("#alignment_idx\tcontig_name\t\n");
        for(int t =0 ; t < alignments_len; t++) {
                bam1_t* b = alignments[t]->record;
                DEBUG_PRINT("%d\t%s\n", t, sam_hdr_tid2name(h, b->core.tid));
        }
}

static struct option long_options[] =
{
    {"inputBam", required_argument, NULL, 'i'},
    {"inputFasta", required_argument, NULL, 'f'},
    {"baq", no_argument, NULL, 'q'},
    {"gapOpen", required_argument, NULL, 'd'},
    {"gapExt", required_argument, NULL, 'e'},
    {"bandwidth", required_argument, NULL, 'b'},
    {"consensus", no_argument, NULL, 'c'},
    {"indelThreshold", required_argument, NULL, 't'},
    {"initQ", required_argument, NULL, 's'},
    {"minQ", required_argument, NULL, 'm'},
    {"primMarginScore", required_argument, NULL, 'p'},
    {"primMarginRandom", required_argument, NULL, 'r'},
    {"minScore", required_argument, NULL, 'n'},
    {"hifi", no_argument, NULL, 'x'},
    {"ont", no_argument, NULL, 'y'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char *argv[]){
	int c;
	bool baq_flag = false;
	bool consensus = false;
	int threshold = 10;
	int min_q = 20;
	int min_score = -50;
	double prim_margin_score = 50;
	double prim_margin_random = 50;
	int set_q = 40;
	double conf_d=1e-4;
	double conf_e=0.1;
	double conf_b=20;
	char* inputPath;
	char* fastaPath;
   	char *program;
	bool preset_ont = false;
	bool preset_hifi = false;
   	(program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
   	while (~(c=getopt_long(argc, argv, "i:p:f:qd:e:b:n:r:m:ct:s:xyh", long_options, NULL))) {
		switch (c) {
                        case 'i':
                                inputPath = optarg;
                                break;
			case 'f':
				fastaPath = optarg;
				break;
			case 'x':
				preset_hifi = true;
				baq_flag = true;
			       	consensus = true;
				threshold = 10; // indel size threshold
			       	conf_d = 1e-4;
			       	conf_e = 0.1;
				conf_b = 20;
			       	min_q = 10;
			        set_q = 40;
			       	prim_margin_score = 50;
				prim_margin_random = 50;
			        min_score = -50;
				break;
			case 'y':
				preset_ont = true;
                                baq_flag = true;
                                consensus = true;
                                threshold = 20; // indel size threshold
                                conf_d = 1e-3;
                                conf_e = 0.1;
                                conf_b = 20;
                                min_q = 10;
                                set_q = 20;
                                prim_margin_score = 10;
				prim_margin_random = 10;
                                min_score = -50;
				break;
			case 'q':
				baq_flag = true;
				break;
			case 'd':
                                conf_d = atof(optarg);
                                break;
			case 'e':
                                conf_e = atof(optarg);
                                break;
			case 'b':
                                conf_b = atof(optarg);
                                break;
			case 'c':
                                consensus = true;
                                break;
			case 't':
				threshold = atoi(optarg);
				break;
			case 's':
                                set_q = atoi(optarg);
                                break;
			case 'm':
                                min_q = atoi(optarg);
                                break;
			case 'p':
				prim_margin_score = atof(optarg);
				break;
			case 'r':
                                prim_margin_random = atof(optarg);
                                break;
			case 'n':
				min_score = atoi(optarg);
				break;
                        default:
                                if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
			help:
                                fprintf(stderr, "\nUsage: %s  -i <INPUT_BAM> -f <FASTA> \n", program);
                                fprintf(stderr, "Options:\n");
                                fprintf(stderr, "         --inputBam, -i         Input bam file\n");
				fprintf(stderr, "         --inputFasta, -f         Input fasta file\n");
				fprintf(stderr, "         --hifi, -x         hifi preset params [-q -c -t10 -d 1e-4 -e 0.1 -b20 -m10 -s40 -p50 -r50 -n -50] (Only one of --hifi or --ont should be enabled)\n");
				fprintf(stderr, "         --ont, -y        ont present params [-q -c -t20 -d 1e-3 -e 0.1 -b20 -m10 -s20 -p50 -r50 -n -50] (Only one of --hifi or --ont should be enabled) \n");
				fprintf(stderr, "         --baq, -q         Calculate BAQ [Disabled by default]\n");
				fprintf(stderr, "         --gapOpen, -d         Gap prob [Default: 1e-4, (for ONT use 1e-2)]\n");
				fprintf(stderr, "         --gapExt, -e         Gap extension [Default: 0.1]\n");
				fprintf(stderr, "         --bandwidth, -b         DP bandwidth [Default: 20]\n");
				fprintf(stderr, "         --consensus, -c         Use consensus confident blocks [Disabled by default]\n");
				fprintf(stderr, "         --indelThreshold, -t         Indel size threshold for confident blocks [Default: 10 (for ONT use 20)]\n");
				fprintf(stderr, "         --initQ, -s         Before calculating BAQ set all base qualities to this number [Default: 40 (for ONT use 20)]\n");
				fprintf(stderr, "         --minQ, -m         Minimum base quality (or BAQ if -q is set) to be considered as a marker  [Default: 20 (for ONT use 10)]\n");
				fprintf(stderr, "         --primMarginScore, -p         Minimum margin between the consistensy score of primary and secondary alignment to select the secondary alignment [Default: 50]\n");
				fprintf(stderr, "         --primMarginRandom, -r         Maximum margin between the consistensy score of primary and secondary alignment to select one randomly [Default: 50]\n");
				
				fprintf(stderr, "         --minScore, -n         Minimum consistency score of the selected secondary alignment [Default: -50]\n");
                                return 1;
		}
	}

	if(preset_ont && preset_hifi){
		fprintf(stderr, "Presets --hifi and --ont cannot be enabled at the same time. Select one of them!\n");
		exit(EXIT_FAILURE);
	}

	faidx_t* fai = fai_load(fastaPath);
	samFile* fp = sam_open(inputPath, "r");
	sam_hdr_t* sam_hdr = sam_hdr_read(fp);
	bam1_t* b = bam_init1();
	char read_name[100];
	char read_name_new[100];
	memset(read_name, '\0', 100);
	memset(read_name_new, '\0', 100);
        uint8_t* quality;
	ptMarker* marker;
	stList* markers = stList_construct3(0, free);
	int alignments_len=0;
	int32_t read_pos;
	ptAlignment* alignments[10];
	int bytes_read;
	const char* contig_name;
	bool conf_blocks_length;
	while(true) {
		bytes_read = sam_read1(fp, sam_hdr, b);
		if (bytes_read > - 1){
			strcpy(read_name_new, bam_get_qname(b));
			if (read_name[0] == '\0') {
				strcpy(read_name, read_name_new);
			}
		}
		// If read name has changed or file is finished
		if ((strcmp(read_name_new, read_name) != 0) || (bytes_read <= -1)){
			// If we have at least one marker and at least two alignments
			// then we can decide which one is the best alignment
			if ((stList_length(markers) > 0) && (alignments_len > 1) && (alignments_len <= 10) && !contain_supp(alignments, alignments_len)){

DEBUG_PRINT("\n@@ READ NAME: %s\n\t$ Number of alignments: %d\t Read l_qseq: %d\n", read_name, alignments_len, alignments[0]->record->core.l_qseq);
				print_contigs(alignments, alignments_len, sam_hdr);
				//DEBUG_PRINT("Initial markers:\n");
				//print_markers(markers);
				remove_all_mismatch_markers(&markers, alignments_len);
				//DEBUG_PRINT("No all-mismatch markers:\n");
                                //print_markers(markers);
				sort_and_fill_markers(&markers, alignments, alignments_len);
				//DEBUG_PRINT("Match-added markers:\n");
                                //print_markers(markers);
				filter_ins_markers(&markers, alignments, alignments_len);
				DEBUG_PRINT("Insertion-removed markers\n");
				print_markers(markers);
				if(markers && stList_length(markers) > 0){
					int flank_margin = 625;
					DEBUG_PRINT("\t# Set confident & flanking blocks\n");
					while(needs_to_find_blocks(alignments, alignments_len, 5000, sam_hdr)){
						flank_margin *= 0.8;
						DEBUG_PRINT("\t# Margin = %d\n", flank_margin);
						set_confident_blocks(alignments, alignments_len, threshold);
						set_flanking_blocks(alignments, alignments_len, markers, flank_margin);
						if (consensus) conf_blocks_length = correct_conf_blocks(alignments, alignments_len, threshold);//TODO: remove if
						if (conf_blocks_length == 0) break;
					}
					if (conf_blocks_length > 0 || consensus == false){
						if(baq_flag){
							DEBUG_PRINT("\t# There are more than 0 markers and confident blocks are not empty \n\t# So calc baq and update quality\n");
							calc_update_baq_all(fai, 
									    alignments, alignments_len,
									    markers, sam_hdr,
								            conf_d, conf_e, conf_b, set_q);
						}
						filter_lowq_markers(&markers, min_q);
						DEBUG_PRINT("High-quality markers\n");
						print_markers(markers);
						DEBUG_PRINT("#Calculate likelihoods\n");
						calc_likelihood(markers, alignments);
					}
				}
			}
			if (alignments_len > 0){ // maybe the previous alignment was unmapped
				// get the best alignment
				int best_idx = get_best_record(alignments, alignments_len, prim_margin_score, min_score, prim_margin_random);
				bam1_t* best = 0 <= best_idx ? alignments[best_idx]->record : NULL;
				// write all alignments without any change if they are either chimeric or best alignment is primary
				if(best &&
				   !contain_supp(alignments, alignments_len) && 
				   (best->core.flag & BAM_FSECONDARY)){
					printf("$\t%s\n", read_name);
					for(int i=0; i<alignments_len; i++){
						// change primary to secondary
						if ((alignments[i]->record->core.flag & BAM_FSECONDARY) == 0){
							printf("*\t");
							alignments[i]->record->core.flag |= BAM_FSECONDARY;
						}
						//change secondary to primary for the best alignment
						else if (i == best_idx){
							printf("@\t");
							alignments[i]->record->core.flag &= ~BAM_FSECONDARY;
						}
						else printf("!\t");
						contig_name = sam_hdr_tid2name(sam_hdr, alignments[i]->record->core.tid);
						printf("%.2f\t%s\t%ld\t%ld\n", alignments[i]->qv, contig_name, alignments[i]->record->core.pos, alignments[i]->rfe);
					}
					printf("\n");
				}
				stList_destruct(markers);
				markers = stList_construct3(0, free);
				for(int i = 0; i < alignments_len; i++){
					ptAlignment_destruct(alignments[i]);
                			alignments[i] = NULL;
        			}
				// initialize for new alignments
				alignments_len = 0;
			}
			strcpy(read_name, read_name_new);
		}
		if (bytes_read <= -1) break; // file is finished so break
		if(b->core.flag & BAM_FUNMAP) continue; // unmapped
		if(alignments_len > 9) continue;
		alignments[alignments_len] = ptAlignment_construct(b, 0.0);
		ptCigarIt* cigar_it = ptCigarIt_construct(b, true, true);
		quality = bam_get_qual(b);
		while(ptCigarIt_next(cigar_it)){
			if (cigar_it->op == BAM_CDIFF) {
				for(int j=0; j < cigar_it->len; j++){
					if(quality[cigar_it->sqs + j] < min_q) continue;
                                        if (bam_is_rev(b)) {
						marker = ptMarker_construct(alignments_len,
                                                                            cigar_it->sqs + j,
                                                                            cigar_it->rde_f - j,
                                                                            quality[cigar_it->sqs + j],
                                                                            false,
									    cigar_it->rfs + j);
					}
					else{
						marker = ptMarker_construct(alignments_len,
                                                                            cigar_it->sqs + j,
                                                                            cigar_it->rds_f + j,
                                                                            quality[cigar_it->sqs + j],
                                                                            false,
									    cigar_it->rfs + j);
					}
					stList_append(markers, marker);
				}
			}
			//set the start coordinates of the alignment
			if (alignments[alignments_len]->rfs == -1 &&
			    (cigar_it->op == BAM_CMATCH ||
                             cigar_it->op == BAM_CEQUAL ||
                             cigar_it->op == BAM_CDIFF)) {
				alignments[alignments_len]->rfs = cigar_it->rfs;
				if(bam_is_rev(b)){
					alignments[alignments_len]->rde_f = cigar_it->rde_f;
				}
				else{ 
					alignments[alignments_len]->rds_f = cigar_it->rds_f;
				}
			}
			//set the end coordinates of the alignment
                        //if the alignment ends with hard or soft clipping
			//alignments[alignments_len]->rfs != -1 is to make sure we have reached
			//the end of the alignment
			if (alignments[alignments_len]->rfe == -1 &&
			    alignments[alignments_len]->rfs != -1 &&
                            (cigar_it->op == BAM_CHARD_CLIP ||
                             cigar_it->op == BAM_CSOFT_CLIP )) {
                                alignments[alignments_len]->rfe = cigar_it->rfe;
				if(bam_is_rev(b)){
					alignments[alignments_len]->rds_f = cigar_it->rde_f + 1;
				}else{
					alignments[alignments_len]->rde_f = cigar_it->rds_f - 1;
				}
                        }
		}
		//set the start coordinates of the alignment
                //if the alignment ends with mis/matches
                if (alignments[alignments_len]->rfe == -1 &&
		    (cigar_it->op == BAM_CMATCH ||
                     cigar_it->op == BAM_CEQUAL ||
                     cigar_it->op == BAM_CDIFF)) {
			alignments[alignments_len]->rfe = cigar_it->rfe;
			if(bam_is_rev(b)){
				alignments[alignments_len]->rds_f = cigar_it->rds_f;
			}
			else{
                		alignments[alignments_len]->rde_f = cigar_it->rde_f;
			}
		}
		alignments_len += 1;
		ptCigarIt_destruct(cigar_it);
	}
	// free memory
	stList_destruct(markers);
	fai_destroy(fai);
	sam_hdr_destroy(sam_hdr);
	sam_close(fp);
	bam_destroy1(b);
}

//main();
