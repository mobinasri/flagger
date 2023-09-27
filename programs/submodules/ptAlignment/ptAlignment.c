#include "ptAlignment.h"

/*
 * @abstract Structure for saving a block
 * @field rfs           Start coordinate on reference (0-based inclusive)
 * @field rfe           End coordinate on reference (0-based inclusive)
 * @field sqs           Start coordinate on sequence (0-based inclusive)
 * @field sqe           End coordinate on sequence (0-based inclusive)
 * @field rds_f         Start coordinate on read's forward strand (0-based inclusive)
 * @field rde_f         End coordinate on read's forward strand (0-based inclusive)
 * @field data		Pointer to some data attached to the block
 * @field free_data	Function for freeing the memory allocated to data
 */
struct ptBlock_ {
    int rfs; // ref start
    int rfe; // ref end
    int sqs; // seq start
    int sqe; // seq end
    int rds_f;
    int rde_f;
    void *data;

    void (*destruct_data)(void *);

    void *(*copy_data)(void *);

    void (*extend_data)(void *, void *);
};

ptAlignment *ptAlignment_construct(bam1_t *record, sam_hdr_t *sam_hdr) {
    ptAlignment *alignment = (ptAlignment *) malloc(sizeof(ptAlignment));
    alignment->record = bam_init1();
    assert(bam_copy1(alignment->record, record) != NULL);
    strcpy(alignment->contig, sam_hdr_tid2name(sam_hdr, record->core.tid));
    alignment->score = 0.0;
    alignment->conf_blocks = NULL;
    alignment->flank_blocks = NULL;
    alignment->mapq = record->core.qual;
    ptAlignment_init_coordinates(alignment);
    return alignment;
}

void ptAlignment_init_coordinates(ptAlignment *alignment) {
    // initialize to -1
    alignment->rfs = -1;
    alignment->rfe = -1;
    alignment->rde_f = -1;
    alignment->rds_f = -1;

    // initialize to 0
    alignment->r_clip = 0;
    alignment->l_clip = 0;

    bam1_t *b = alignment->record;
    ptCigarIt *cigar_it = ptCigarIt_construct(b, true, true);
    uint8_t *quality = bam_get_qual(b);
    while (ptCigarIt_next(cigar_it)) {
        //set the start coordinates of the alignment
        if (alignment->rfs == -1 &&
            (cigar_it->op == BAM_CMATCH ||
             cigar_it->op == BAM_CEQUAL ||
             cigar_it->op == BAM_CDIFF)) {
            alignment->rfs = cigar_it->rfs;
            if (bam_is_rev(b)) {
                alignment->rde_f = cigar_it->rde_f;
            } else {
                alignment->rds_f = cigar_it->rds_f;
            }
        }
        //set the size of the left clipping (left w.r.t. ref)
        if (alignment->rfs == -1 &&
            (cigar_it->op == BAM_CHARD_CLIP ||
             cigar_it->op == BAM_CSOFT_CLIP)) {
            alignment->l_clip += cigar_it->len;
        }

        //set the size of the right clipping (right w.r.t. ref)
        if (alignment->rfs != -1 &&
            (cigar_it->op == BAM_CHARD_CLIP ||
             cigar_it->op == BAM_CSOFT_CLIP)) {
            alignment->r_clip += cigar_it->len;
        }


        //set the end coordinates of the alignment
        //if the alignment ends with hard or soft clipping
        //alignment->rfs != -1 is to make sure we have reached
        //the end of the alignment
        if (alignment->rfe == -1 &&
            alignment->rfs != -1 &&
            (cigar_it->op == BAM_CHARD_CLIP ||
             cigar_it->op == BAM_CSOFT_CLIP)) {
            alignment->rfe = cigar_it->rfe;
            if (bam_is_rev(b)) {
                alignment->rds_f = cigar_it->rde_f + 1;
            } else {
                alignment->rde_f = cigar_it->rds_f - 1;
            }
        }
    }
    //set the end coordinates of the alignment
    //if the alignment ends with mis/matches
    if (alignment->rfe == -1 &&
        (cigar_it->op == BAM_CMATCH ||
         cigar_it->op == BAM_CEQUAL ||
         cigar_it->op == BAM_CDIFF)) {
        alignment->rfe = cigar_it->rfe;
        if (bam_is_rev(b)) {
            alignment->rds_f = cigar_it->rds_f;
        } else {
            alignment->rde_f = cigar_it->rde_f;
        }
    }
    ptCigarIt_destruct(cigar_it);
}

void ptAlignment_destruct(ptAlignment *alignment) {
    if (alignment->record) {
        bam_destroy1(alignment->record);
        alignment->record = NULL;
    }
    if (alignment->conf_blocks) {
        stList_destruct(alignment->conf_blocks);
        alignment->conf_blocks = NULL;
    }
    if (alignment->flank_blocks) {
        stList_destruct(alignment->flank_blocks);
        alignment->flank_blocks = NULL;
    }
    free(alignment);
}

int ptAlignment_supplementary_count(ptAlignment **alignments, int alignments_len) {
    int supp_count = 0;
    for (int i = 0; i < alignments_len; i++) {
        if (alignments[i]->record->core.flag & BAM_FSUPPLEMENTARY) supp_count += 1;
    }
    return supp_count;
}

int ptAlignment_primary_count(ptAlignment **alignments, int alignments_len) {
    int primary_count = 0;
    for (int i = 0; i < alignments_len; i++) {
        if ((alignments[i]->record->core.flag & BAM_FSECONDARY) == 0) primary_count += 1;
    }
    return primary_count;
}

void print_contigs(ptAlignment **alignments, int alignments_len) {
    DEBUG_PRINT("\nalignments:\n");
    DEBUG_PRINT("#alignment_idx\tcontig_name\t\n");
    for (int t = 0; t < alignments_len; t++) {
        DEBUG_PRINT("%d\t%s\n", t, alignments[t]->contig);
    }
}

int get_best_record_index(ptAlignment **alignments, int alignments_len, double prim_margin, double min_score,
                          double prim_margin_random) {
    assert(alignments_len > 0);
    if (alignments_len == 1) return 0;
    double max_score = -DBL_MAX;
    int max_idx = -1;
    double prim_score = -DBL_MAX;
    int prim_idx = -1;
    for (int i = 0; i < alignments_len; i++) {
        if ((alignments[i]->record->core.flag & BAM_FSECONDARY) == 0) {
            prim_idx = i;
            prim_score = alignments[i]->score;
        } else if (max_score < alignments[i]->score) {
            max_idx = i;
            max_score = alignments[i]->score;
        }
    }
    // make a list of all secondary alignments with max score
    stList *sec_alignment_indices_with_max_score = stList_construct3(0, free);
    for (int i = 0; i < alignments_len; i++) {
        if (((alignments[i]->record->core.flag & BAM_FSECONDARY) != 0) && (max_score <= alignments[i]->score)) {
            int *index = malloc(sizeof(int));
            *index = i;
            stList_append(sec_alignment_indices_with_max_score, index);
        }
    }
    // if there are multiple secondary alignments with the maximum score,
    // select one of them randomly
    if (stList_length(sec_alignment_indices_with_max_score) > 1) {
        int *max_idx_ptr = stList_get(sec_alignment_indices_with_max_score,
                                      rand() % stList_length(sec_alignment_indices_with_max_score));
        max_idx = *max_idx_ptr;
    }
    stList_destruct(sec_alignment_indices_with_max_score);
    int rnd = rand() % 2; // 50% chance for rnd=0 (same for rnd=1)
    if (abs(max_score - prim_score) < prim_margin_random) {
        return rnd == 0 ? prim_idx : max_idx; // if the scores were closer than prim_margin_random return one randomly
    }
    if (prim_idx == -1 || max_score <= (prim_score + prim_margin) || max_score < min_score) return prim_idx;
    else return max_idx;
}

int get_primary_index(ptAlignment **alignments, int alignments_len) {
    assert(alignments_len > 0);
    for (int i = 0; i < alignments_len; i++) {
        if ((alignments[i]->record->core.flag & BAM_FSECONDARY) == 0) {
            return i;
        }
    }
    // if no primary exists return -1
    return -1;
}


