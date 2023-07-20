#include "ptMarker.h"
#include "ptBlock.h"
#include "cigar_it.h"


ptMarker *ptMarker_construct(int32_t alignment_idx, int32_t base_idx, int32_t read_pos_f, uint8_t base_q, bool is_match,
                             int32_t ref_pos) {
    ptMarker *marker = (ptMarker *) malloc(sizeof(ptMarker));
    marker->alignment_idx = alignment_idx;
    marker->is_match = is_match;
    marker->read_pos_f = read_pos_f;
    marker->base_idx = base_idx;
    marker->base_q = base_q;
    marker->ref_pos = ref_pos;
    return marker;
}


ptMarker *ptMarker_copy(ptMarker *src) {
    ptMarker *dest = (ptMarker *) malloc(sizeof(ptMarker));
    dest->alignment_idx = src->alignment_idx;
    dest->is_match = src->is_match;
    dest->read_pos_f = src->read_pos_f;
    dest->base_idx = src->base_idx;
    dest->base_q = src->base_q;
    dest->ref_pos = src->ref_pos;
    return dest;
}


int ptMarker_cmp(const void *a, const void *b) {
    ptMarker *m1 = (ptMarker *) a;
    ptMarker *m2 = (ptMarker *) b;
    if (m1->read_pos_f == m2->read_pos_f) {
        return m1->alignment_idx - m2->alignment_idx;
    } else {
        return m1->read_pos_f - m2->read_pos_f;
    }
}


stList *ptMarker_get_initial_markers(ptAlignment **alignments, int alignments_len, int min_q) {
    stList *markers = stList_construct3(0, free);
    ptMarker *marker = NULL;
    for (int i = 0; i < alignments_len; i++) {
        bam1_t* b = alignments[i]->record;
        ptCigarIt *cigar_it = ptCigarIt_construct(b, true, true);
        uint8_t *quality = bam_get_qual(b);
        while (ptCigarIt_next(cigar_it)) {
            if (cigar_it->op == BAM_CDIFF) {
                for (int j = 0; j < cigar_it->len; j++) {
                    if (quality[cigar_it->sqs + j] < min_q) continue;
                    if (bam_is_rev(b)) {
                        marker = ptMarker_construct(i,
                                                    cigar_it->sqs + j,
                                                    cigar_it->rde_f - j,
                                                    quality[cigar_it->sqs + j],
                                                    false,
                                                    cigar_it->rfs + j);
                    } else {
                        marker = ptMarker_construct(i,
                                                    cigar_it->sqs + j,
                                                    cigar_it->rds_f + j,
                                                    quality[cigar_it->sqs + j],
                                                    false,
                                                    cigar_it->rfs + j);
                    }
                    stList_append(markers, marker);
                }
            }
        }
        ptCigarIt_destruct(cigar_it);
    }
    return markers;
}

ptMarker *ptMarker_construct_match(ptAlignment **alignments, int32_t alignment_idx, int32_t read_pos_f) {
    bam1_t *record = alignments[alignment_idx]->record;
    uint8_t *q = bam_get_qual(record);
    ptMarker *match;
    uint32_t *cigar = bam_get_cigar(record);
    int lclip = 0;
    int rclip = 0;
    if (bam_cigar_op(cigar[0]) == BAM_CHARD_CLIP) {
        lclip = bam_cigar_oplen(cigar[0]);
    }
    if (bam_cigar_op(cigar[record->core.n_cigar - 1]) == BAM_CHARD_CLIP) {
        rclip = bam_cigar_oplen(cigar[record->core.n_cigar - 1]);
    }
    // construct the match marker
    if (bam_is_rev(record)) {
        match = ptMarker_construct(alignment_idx,
                                   record->core.l_qseq + rclip - read_pos_f - 1,
                                   read_pos_f,
                                   q[record->core.l_qseq + rclip - read_pos_f - 1],
                                   true,
                                   -1);
    } else {
        match = ptMarker_construct(alignment_idx,
                                   read_pos_f - lclip,
                                   read_pos_f,
                                   q[read_pos_f - lclip],
                                   true,
                                   -1);
    }
    return match;
}


void filter_lowq_markers(stList **markers_p, int threshold) {
    stList *markers = *markers_p;
    int markers_len = stList_length(markers);
    if (markers_len == 0) return;
    stList *keep_markers = stList_construct3(0, free);
    ptMarker *pre_marker = NULL;
    ptMarker *marker;
    bool keep_flag = false;
    int idx_s = 0;
    int min_q = 100;
    ptMarker *marker_copy;
    for (int j = 0; j < markers_len; j++) {
        marker = stList_get(markers, j);
        if (pre_marker && marker->read_pos_f != pre_marker->read_pos_f) {
            if (min_q > threshold) {
                for (int k = idx_s; k < j; k++) {
                    // copy all markers with the same read_pos_f and add them
                    marker_copy = ptMarker_copy(stList_get(markers, k));
                    marker_copy->base_q = min_q;
                    stList_append(keep_markers, marker_copy);
                }
            }
            idx_s = j;
            // reset min_q for the next read position
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
    if (min_q > threshold) {
        for (int k = idx_s; k < markers_len; k++) {
            marker_copy = ptMarker_copy(stList_get(markers, k));
            marker_copy->base_q = min_q;
            stList_append(keep_markers, marker_copy);
        }
    }
    stList_destruct(markers);
    *markers_p = keep_markers;
}


void filter_ins_markers(stList **markers_p, ptAlignment **alignments, int alignments_len) {
    stList *markers = *markers_p;
    int markers_len = stList_length(markers);
    if (markers_len == 0) return;
    bool *markers_keep_flags = (bool *) malloc(markers_len * sizeof(bool));
    for (int i = 0; i < markers_len; i++) markers_keep_flags[i] = true;
    bam1_t *b;
    ptMarker *marker;
    ptCigarIt *cigar_it;
    int j; // marker index
    for (int i = 0; i < alignments_len; i++) {
        b = alignments[i]->record;
        j = bam_is_rev(b) ? stList_length(markers) - 1 : 0;
        marker = stList_get(markers, j);
        cigar_it = ptCigarIt_construct(b, true, true);
        //iterate over cigars
        while (ptCigarIt_next(cigar_it)) {
            //check if marker is located within the cigar operation
            while (marker &&
                   (cigar_it->rds_f <= marker->read_pos_f) &&
                   (cigar_it->rde_f >= marker->read_pos_f)) {
                //check if the marker is located within an insertion
                if (cigar_it->op == BAM_CINS ||
                    cigar_it->op == BAM_CSOFT_CLIP ||
                    cigar_it->op == BAM_CHARD_CLIP) {
                    markers_keep_flags[j] = false;
                }
                //update reference positions in the meanwhile!
                if (cigar_it->op == BAM_CEQUAL && marker->alignment_idx == i) {
                    marker->ref_pos = bam_is_rev(b) ? cigar_it->rfs + cigar_it->rde_f - marker->read_pos_f :
                                      cigar_it->rfs + marker->read_pos_f - cigar_it->rds_f;
                }
                j += bam_is_rev(b) ? -1 : 1;
                marker = ((j < markers_len) && (j >= 0)) ? stList_get(markers, j) : NULL;
            }
        }//end of iteration over cigar ops
        ptCigarIt_destruct(cigar_it);
    }// end of iteration over alignments
    // make a list for saving the remaining markers
    stList *keep_markers = stList_construct3(0, free);
    for (j = 0; j < markers_len; j++) {
        if (markers_keep_flags[j]) {
            marker = stList_get(markers, j);
            stList_append(keep_markers, ptMarker_copy(marker));
        }
    }
    stList_destruct(markers);
    *markers_p = keep_markers;
    // free memory
    free(markers_keep_flags);
}


void remove_all_mismatch_markers(stList **markers_p, int alignments_len) {
    stList *markers = *markers_p;
    stList_sort(markers, ptMarker_cmp);
    int markers_len = stList_length(markers);
    if (markers_len == 0) return;
    bool *markers_keep_flags = (bool *) malloc(markers_len * sizeof(bool));
    memset(markers_keep_flags, 1, markers_len);
    ptMarker *marker;
    ptMarker *pre_marker = NULL;
    int occ = 0;
    for (int i = 0; i < markers_len; i++) {
        marker = stList_get(markers, i);
        if (pre_marker &&
            pre_marker->read_pos_f < marker->read_pos_f) {
            if (occ == alignments_len) {// if all markers are mismatch
                memset(markers_keep_flags + i - alignments_len, 0, alignments_len);
            }
            occ = 0;
        }
        pre_marker = marker;
        occ += 1;
    }
    if (occ == alignments_len) {
        memset(markers_keep_flags + markers_len - alignments_len, 0, alignments_len);
    }
    // make a list for saving the remaining markers
    stList *keep_markers = stList_construct3(0, free);
    for (int j = 0; j < markers_len; j++) {
        if (markers_keep_flags[j]) {
            marker = stList_get(markers, j);
            stList_append(keep_markers, ptMarker_copy(marker));
        }
    }
    // free previous list of markers
    stList_destruct(markers);
    // update the list
    *markers_p = keep_markers;
    // free
    free(markers_keep_flags);
}


void sort_and_fill_markers(stList **markers_p, ptAlignment **alignments, int alignments_len) {
    stList *markers = *markers_p;
    stList_sort(markers, ptMarker_cmp);
    int idx = 0;
    ptMarker *match;
    ptMarker *marker;
    ptMarker *pre_marker = NULL;

    // markers_len is the number of initially given markers (all mismatch)
    // it does not get updated while adding match markers
    int markers_len = stList_length(markers);
    for (int64_t i = 0; i < markers_len; i++) {
        marker = stList_get(markers, i);
        // none of the initially given markers should be match
        assert(!marker->is_match);
        // if the previous location of read is passed
        // then we start adding matched markers for the
        // remaining alignments (if it existed any)
        if (pre_marker &&
            pre_marker->read_pos_f < marker->read_pos_f) {
            for (int64_t j = idx; j < alignments_len; j++) {
                match = ptMarker_construct_match(alignments, j, pre_marker->read_pos_f);
                stList_append(markers, match);
            }
            idx = 0;
        }
        pre_marker = marker;
        // add match markers for the alignments whose indices start
        // from idx and end one before marker->alignment_idx
        for (int64_t j = idx; j < marker->alignment_idx; j++) {
            match = ptMarker_construct_match(alignments, j, marker->read_pos_f);
            stList_append(markers, match);
        }
        // update the first index of the next alignments that may lack this marker
        idx = marker->alignment_idx + 1;
    }
    if (markers_len > 0) {
        // add the last matches if there still exist to be added
        for (int64_t j = idx; j < alignments_len; j++) {
            match = ptMarker_construct_match(alignments, j, marker->read_pos_f);
            stList_append(markers, match);
        }
    }
    stList_sort(markers, ptMarker_cmp); // sort again to locate the matches in their correct positions
}


double reverse_quality(uint8_t q) {
    if (q >= 93) return 0;
    if (q == 0) return 93;
    double p = 1 - pow(10, (double) q / -10);
    double rev_q = -10 * log(p);
    return rev_q;
}


void calc_alignment_score(stList *markers, ptAlignment **alignments) {
    ptMarker *marker;
    for (int64_t j = 0; j < stList_length(markers); j++) {
        marker = stList_get(markers, j);
        if (marker->is_match) {

            // p(all matches are real match) = (1-e_1) * (1-e_2) * ... (1-e_n1)
            // q(all matches are real match) = -rev(q_1) + -rev(q_2) + ... + -rev(q_n1)
            alignments[marker->alignment_idx]->score += -1 * reverse_quality(marker->base_q);
        } else {
            // p(all mismtaches are real match) = (e1/3) * (e2/3) * ... (en/3)
            // q(all mismatches are real match) = (-q1 - 10log(3)) + (-q2 - 10log(3)) + ... + (-qn - 10log(3))
            alignments[marker->alignment_idx]->score += -1 * marker->base_q - 10 * log(3);
        }
        // p(all markers are real match for this alignment) = p(all matches are real match) * p(all mismtaches are real match)
        // p(all markers are real match for this alignment) = p(this alignment is correct)
        // final q = q(all matches are real match) + q(all mismatches are real match)
    }
}


stList *find_confident_blocks(ptAlignment *alignment, int threshold) {
    bam1_t *b = alignment->record;
    stList *conf_blocks = stList_construct3(0, ptBlock_destruct);
    ptCigarIt *cigar_it = ptCigarIt_construct(b, true, true);
    int conf_sqs = 0; // the seq start of the confident block
    int conf_rfs = b->core.pos; // the ref start of the confident block
    int conf_rd_f = bam_is_rev(b) ? cigar_it->rde_f : cigar_it->rds_f;
    ptBlock *block = NULL;
    while (ptCigarIt_next(cigar_it)) {
        switch (cigar_it->op) {
            case BAM_CINS:
            case BAM_CDEL:
                if (cigar_it->len > threshold && conf_sqs < cigar_it->sqs && conf_rfs < cigar_it->rfs) {
                    if (bam_is_rev(b)) {
                        block = ptBlock_construct(conf_rfs, cigar_it->rfs - 1,
                                                  conf_sqs, cigar_it->sqs - 1,
                                                  cigar_it->rde_f + 1, conf_rd_f);
                    } else {
                        block = ptBlock_construct(conf_rfs, cigar_it->rfs - 1,
                                                  conf_sqs, cigar_it->sqs - 1,
                                                  conf_rd_f, cigar_it->rds_f - 1);
                    }
                    stList_append(conf_blocks, block);
                }
                if (cigar_it->len > threshold) {
                    conf_sqs = cigar_it->sqe + 1;
                    conf_rfs = cigar_it->rfe + 1;
                    conf_rd_f = bam_is_rev(b) ? cigar_it->rds_f - 1 : cigar_it->rde_f + 1;
                }
                break;
            case BAM_CSOFT_CLIP:
            case BAM_CHARD_CLIP:
                if (conf_sqs < cigar_it->sqs &&
                    conf_rfs < cigar_it->rfs) {
                    if (bam_is_rev(b)) {
                        block = ptBlock_construct(conf_rfs, cigar_it->rfs - 1,
                                                  conf_sqs, cigar_it->sqs - 1,
                                                  cigar_it->rde_f + 1, conf_rd_f);
                    } else {
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
    if (conf_sqs <= cigar_it->sqe) {
        if (bam_is_rev(b)) {
            block = ptBlock_construct(conf_rfs, cigar_it->rfe,
                                      conf_sqs, cigar_it->sqe,
                                      cigar_it->rds_f, conf_rd_f);
        } else {
            block = ptBlock_construct(conf_rfs, cigar_it->rfe,
                                      conf_sqs, cigar_it->sqe,
                                      conf_rd_f, cigar_it->rde_f);
        }
        stList_append(conf_blocks, block);
    }
    ptCigarIt_destruct(cigar_it);
    return conf_blocks;
}


stList *intersect_by_rd_f(stList *blocks1, stList *blocks2) {
    stList *blocks_intersect = stList_construct3(0, ptBlock_destruct);
    if (stList_length(blocks1) == 0 || stList_length(blocks2) == 0) return blocks_intersect;
    ptBlock *b1;
    ptBlock *b2;
    ptBlock *b;
    int j = 0;
    b2 = stList_get(blocks2, j);
    for (int i = 0; i < stList_length(blocks1); i++) {
        b1 = stList_get(blocks1, i);
        // Find the first block (b2) in the 2nd set whose end point is after
        // the start point of the current block (b1) in the 1st set
        while (b2 && b2->rde_f < b1->rds_f) {
            j++;
            b2 = j == stList_length(blocks2) ? NULL : stList_get(blocks2, j);
        }
        // Because of the previous loop, b2's end point is now after b1's
        // start point. If b2's start point is also before the b1's end point
        // it means there exist an overlap so we save that overlap in
        // blocks_intersect and go to the next block (update b2)
        while (b2 && b2->rds_f < b1->rde_f) {
            b = ptBlock_construct(-1,
                                  -1,
                                  -1,
                                  -1,
                                  max(b1->rds_f, b2->rds_f),
                                  min(b1->rde_f, b2->rde_f));
            stList_append(blocks_intersect, b);
            // Go to next b2 if current b2's end point is
            // before b1's end point
            if (b2->rde_f <= b1->rde_f) {
                j++;
                b2 = j == stList_length(blocks2) ? NULL : stList_get(blocks2, j);
            } else break;
        }
    }
    return blocks_intersect;
}


void set_confident_blocks(ptAlignment **alignments, int alignments_len, int threshold) {
    for (int i = 0; i < alignments_len; i++) {

        alignments[i]->conf_blocks = find_confident_blocks(alignments[i], threshold);
    }
}


stList *find_flanking_blocks(ptAlignment *alignment, stList *markers, int margin) {
    bam1_t *b = alignment->record;
    stList *flank_blocks = stList_construct3(0, ptBlock_destruct);
    ptCigarIt *cigar_it = ptCigarIt_construct(b, true, true);
    // get the first merker for initializing start and end
    int i = 0;
    ptMarker *marker = stList_get(markers, i);
    int start = max(alignment->rds_f, marker->read_pos_f - margin);
    int end = min(alignment->rde_f, marker->read_pos_f + margin);
    i++;
    int curr_start;
    int curr_end;
    ptBlock *block = NULL;
    while (i < stList_length(markers)) {
        marker = stList_get(markers, i);
        curr_start = max(alignment->rds_f, marker->read_pos_f - margin);
        curr_end = min(alignment->rde_f, marker->read_pos_f + margin);
        if (curr_start < end) {
            end = curr_end;
        } else {
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


void set_flanking_blocks(ptAlignment **alignments, int alignments_len, stList *markers, int margin) {
    for (int i = 0; i < alignments_len; i++) {
        if(alignments[i]->flank_blocks) stList_destruct(alignments[i]->flank_blocks);
        alignments[i]->flank_blocks = find_flanking_blocks(alignments[i], markers, margin);
    }
}


int correct_conf_blocks(ptAlignment **alignments, int alignments_len, int threshold) {
    assert(alignments_len > 0);
    stList_sort(alignments[0]->conf_blocks, ptBlock_cmp_rds_f);
    stList *blocks = stList_copy(alignments[0]->conf_blocks, NULL);
    ptBlock *block;
    stList *blocks_new;
    //intersect confident blocks
    for (int i = 1; i < alignments_len; i++) {
        stList_sort(alignments[i]->conf_blocks, ptBlock_cmp_rds_f);
        blocks_new = intersect_by_rd_f(blocks, alignments[i]->conf_blocks);
        stList_destruct(blocks);
        blocks = blocks_new;
    }
    //intersect flanking blocks
    for (int i = 0; i < alignments_len; i++) {
        stList_sort(alignments[i]->flank_blocks, ptBlock_cmp_rds_f);
        blocks_new = intersect_by_rd_f(blocks, alignments[i]->flank_blocks);
        stList_destruct(blocks);
        blocks = blocks_new;
    }
    if (stList_length(blocks) == 0) {
        DEBUG_PRINT("\t\t\t### No Consensus Blocks!\n");
        for (int i = 0; i < alignments_len; i++) {
            stList_destruct(alignments[i]->conf_blocks);
            alignments[i]->conf_blocks = stList_construct3(0, ptBlock_destruct);
        }
        stList_destruct(blocks);
        return 0;
    }
    stList *corrected_conf_blocks;
    // For each alignment project the consensus confident blocks to the ref and seq coordinates
    // The main function of the loop below is to fill the attributes rfs, rfe, sqs, and sqe for each block
    // based on the previously calculated rds_f and rde_f
    for (int i = 0; i < alignments_len; i++) {
        bam1_t *b = alignments[i]->record;
        int j = bam_is_rev(b) ? stList_length(blocks) - 1 : 0;
        corrected_conf_blocks = stList_construct3(0, ptBlock_destruct);
        ptCigarIt *cigar_it = ptCigarIt_construct(alignments[i]->record, true, true);
        int block_rds_f, block_rde_f, cigar_rds_f, cigar_rde_f, rfs, rfe, sqs, sqe;
        bool del_flag = false;
        block = stList_get(blocks, j);
        // reverse each block interval to make it work with the projecting algorithm below
        if (bam_is_rev(b)) {
            block_rds_f = -1 * block->rde_f;
            block_rde_f = -1 * block->rds_f;
        } else {
            block_rds_f = block->rds_f;
            block_rde_f = block->rde_f;
        }
        //iterate over cigar operations from left to right (w.r.t reference)
        // There is no need to check for BAM_CHARD_CLIP and BAM_CSOFT_CLIP since
        // we know that consensus blocks will have no overlap with those operations
        // Since these blocks are the outcomes of calling find_confident_blocks() and
        // then they are intersected with intersect_by_rd_f()
        while (ptCigarIt_next(cigar_it)) {
            // reverse each cigar operation interval to make it work with the projecting algorithm below
            if (bam_is_rev(b)) {
                cigar_rds_f = -1 * cigar_it->rde_f;
                cigar_rde_f = -1 * cigar_it->rds_f;
            } else {
                cigar_rds_f = cigar_it->rds_f;
                cigar_rde_f = cigar_it->rde_f;
            }
            // match/ mismatch/ insertion
            // for M and I the main algorithm is the same except
            // some minor parts that are corrected by
            // the conditional statement "cigar_it->op == BAM_CINS ? : "
            if (cigar_it->op == BAM_CMATCH ||
                cigar_it->op == BAM_CEQUAL ||
                cigar_it->op == BAM_CDIFF ||
                cigar_it->op == BAM_CINS) {
                //The loop below iterate the blocks overlap with the current cigar operation
                //There may exist multiple blocks within the current op
                while (block && block_rde_f <= cigar_rde_f) {
                    // update start locations of the projected coordinates
                    if (cigar_rds_f <= block_rds_f && !(del_flag && cigar_rds_f == block_rds_f)) {
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
                    if (bam_is_rev(b) && j > 0) {
                        j--;
                        block = stList_get(blocks, j);
                        block_rds_f = -1 * block->rde_f;
                        block_rde_f = -1 * block->rds_f;
                    } else if (!bam_is_rev(b) && j < stList_length(blocks) - 1) {
                        j++;
                        block = stList_get(blocks, j);
                        block_rds_f = block->rds_f;
                        block_rde_f = block->rde_f;
                    } else if (j == 0 || j == stList_length(blocks) - 1) {
                        block = NULL;
                    }
                }
                if (block == NULL) break;// no more block remaining so break the cigar iteration
                // !(del_flag && cigar_rds_f == block_rds_f) means that,
                // del_flag and cigar_rds_f == block_rds_f cannot be true at the same time
                // if both were true it means that there is a deletion (shorter than threshold)
                // on the edge of the block (This is a rare event can happen when there is a deletion
                // right before an insertion)
                if (cigar_rds_f <= block_rds_f &&
                    block_rds_f <= cigar_rde_f &&
                    !(del_flag && cigar_rds_f == block_rds_f)) {
                    rfs = cigar_it->op == BAM_CINS ? cigar_it->rfs : cigar_it->rfs + (block_rds_f - cigar_rds_f);
                    sqs = cigar_it->sqs + (block_rds_f - cigar_rds_f);
                }
                del_flag = false;
            }
                //deletion
            else if (cigar_it->op == BAM_CDEL) {
                int prev_j = bam_is_rev(b) ? j + 1 : j - 1;
                ptBlock *prev_block = (0 <= prev_j && prev_j <= stList_length(blocks) - 1) ? stList_get(blocks, prev_j)
                                                                                           : NULL;
                // if this is a short deletion where the previous block ended right before it
                // ref end coordinate has to be corrected
                if (prev_block && cigar_it->len <= threshold) {
                    if ((bam_is_rev(b) && prev_block->rds_f == cigar_it->rds_f) ||
                        (!bam_is_rev(b) && prev_block->rde_f == cigar_it->rde_f)) {
                        prev_block->rfe = cigar_it->rfe;
                    }
                }
                // if this is a short deletion where the current block starts right after it
                // ref start coordinate should be updated
                if (block && block_rds_f == cigar_rds_f && cigar_it->len <= threshold) {
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
    int len_blocks = stList_length(blocks);
    stList_destruct(blocks);
    return len_blocks;
}

bool needs_to_find_blocks(ptAlignment **alignments, int alignments_len, int threshold, sam_hdr_t *sam_hdr) {
    stList *blocks;
    ptBlock *block;
    bool flag = false;
    char *contig_name;
    for (int j = 0; j < alignments_len; j++) {
        blocks = alignments[j]->conf_blocks;
        contig_name = alignments[j]->contig;
        if (blocks == NULL) return true;
        if (stList_length(blocks) == 0) return true;
        for (int i = 0; i < stList_length(blocks); i++) {
            block = stList_get(blocks, i);
            DEBUG_PRINT("\t\t\t###     block#%d: seq[%d:%d] ref[%s:%d-%d]\n", i, block->sqs, block->sqe, contig_name,
                        block->rfs, block->rfe);
            if ((block->sqe - block->sqs) > threshold || (block->rfe - block->rfs) > threshold) flag = true;
        }
    }
    return flag;
}


void
calc_local_baq(const faidx_t *fai, const char *contig_name, ptAlignment *alignment, int alignment_idx, stList *markers,
               double conf_d, double conf_e, double conf_bw, int set_q) {

    uint8_t *tseq; // translated seq A=>0,C=>1,G=>2,T=>3,other=>4
    uint8_t *tref; // translated ref
    uint8_t *bq;
    uint8_t *block_qual;
    int *state;
    uint8_t *q; // Probability of incorrect alignment from probaln_glocal()
    probaln_par_t conf = {conf_d, conf_e, conf_bw};
    char *ref;
    char *reg;
    bam1_t *b = alignment->record;
    uint8_t *seq = bam_get_seq(b);
    uint8_t *qual = bam_get_qual(b);
    stList *blocks = alignment->conf_blocks;
    int markers_len = stList_length(markers);
    int j = bam_is_rev(b) ? stList_length(markers) - 1 : 0; // marker index
    ptMarker *marker = stList_get(markers, j);
    ptCigarIt *cigar_it = ptCigarIt_construct(b, true, true);
    ptBlock *block;
    int ref_len;
    int seq_len;
    DEBUG_PRINT("\t\t\t## Number of Blocks: %ld\n", stList_length(blocks));
    // block margin determines how close a marker can be to the ends of a confident block
    // if it was closer than this margin the BAQ value will be set to 0
    int block_margin = 10;
    // Start iterating over confident blocks
    for (int i = 0; i < stList_length(blocks); i++) {
        block = stList_get(blocks, i);
        // Find the first cigar operation that has overlap with the block
        while ((cigar_it->sqe < block->sqs) || (cigar_it->rfe < block->rfs)) {
            if (ptCigarIt_next(cigar_it) == 0) break;
        }// end of while cigar_it->sqs should be now less than or equal to block->sqs
        //assert(cigar_it->sqs <= block->sqs);
        //assert(cigar_it->rfs <= block->rfs);
        // Move forward untill the marker is after start of the block + some margin
        // this loop ensures that the alignment index matches the marker's
        while (marker &&
               ((marker->base_idx < block->sqs + block_margin) ||
                (marker->alignment_idx != alignment_idx))) {
            // set quality to zero if it's close to the edge of the block by the block margin
            // The alignment index should match
            if ((marker->alignment_idx == alignment_idx) &&
                (block->sqs <= marker->base_idx)) {
                qual[marker->base_idx] = 0;
            }
            j += bam_is_rev(alignment->record) ? -1 : 1;
            marker = ((j < markers_len) && (j >= 0)) ? stList_get(markers, j) : NULL;
        }// end of while. The marker is now the first one located after the start of the block + some margin

        // If there exists a marker and that is located within the block excluding left and right margins
        // then get realign ref seq and read seq and realign them against each other
        // and calculate base alignment qualities for all bases of the block
        if (marker &&
            (marker->base_idx <= block->sqe - block_margin) &&
            (block->sqs + block_margin <= marker->base_idx)) {
            seq_len = block->sqe - block->sqs + 1;
            ref_len = block->rfe - block->rfs + 1;
            //allocate and fetch read sequence
            tseq = (uint8_t *) malloc(seq_len);
            for (int k = 0; k < seq_len; k++)
                tseq[k] = seq_nt16_int[bam_seqi(seq, block->sqs + k)];

            //allocate and fetch reference sequence
            tref = (uint8_t *) malloc(ref_len);
            reg = malloc(200);
            memset(reg, '\0', 200);
            sprintf(reg, "{%s}:%d-%d", contig_name, block->rfs + 1, block->rfe + 1);
            int len;
            ref = fai_fetch(fai, reg, &len);
            assert(len == (block->rfe - block->rfs + 1));
            for (int k = 0; k < ref_len; k++)
                tref[k] = seq_nt16_int[seq_nt16_table[(unsigned char) ref[k]]];

            //allocate an array for the base qualities of this confident block
            block_qual = (uint8_t *) malloc(seq_len);
            for (int k = 0; k < seq_len; k++)
                block_qual[k] = set_q;

            //allocate neccessary arrays for HMM BAQ
            state = (int *) malloc((block->sqe - block->sqs + 1) * sizeof(int));
            q = (uint8_t *) malloc(block->sqe - block->sqs + 1);
            conf.bw = abs(ref_len - seq_len) + conf_bw;
            if (probaln_glocal(tref, ref_len,
                               tseq, seq_len,
                               block_qual, &conf, state, q) == INT_MIN) {
                fprintf(stderr, "probaln_glocal ERROR\n");
                fprintf(stderr, "%s:%d-%d", contig_name, block->rfs + 1, block->rfe + 1);
            }
            // the state and q are now updated if there is any marker located within the block
            //apply BAQ to the quality array
            bq = (uint8_t *) malloc(seq_len);
            memcpy(bq, block_qual, seq_len);
            int x;
            int y;
            while (cigar_it->sqs <= block->sqe || cigar_it->rfs <= block->rfe) {
                x = cigar_it->rfs - block->rfs;
                x = x < 0 ? 0 : x;
                y = cigar_it->sqs - block->sqs;
                y = y < 0 ? 0 : y;
                if (cigar_it->op == BAM_CMATCH ||
                    cigar_it->op == BAM_CEQUAL ||
                    cigar_it->op == BAM_CDIFF) {
                    int len = min(cigar_it->len, min(cigar_it->sqe, block->sqe) - max(cigar_it->sqs, block->sqs) + 1);
                    for (int t = y; t < (y + len); t++) {
                        assert(t < seq_len);
                        if (((state[t] & 3) != 0) || (state[t] >> 2 != x + (t - y))) bq[t] = 0;
                        else bq[t] = qual[block->sqs + t] < q[t] ? qual[block->sqs + t] : q[t];
                    }
                }
                if (cigar_it->sqe <= block->sqe || cigar_it->rfe <= block->rfe) {
                    if (ptCigarIt_next(cigar_it) == 0) break;
                } else break;
            }
            for (int k = block_margin; k < seq_len - block_margin; k++) qual[block->sqs + k] = bq[k] < 94 ? bq[k] : 93;
            free(tseq);
            free(tref);
            free(state);
            free(q);
            free(bq);
            free(block_qual);
            free(ref);
            free(reg);
        }
        // for the markers located within the ending margin of the confident block
        while (marker &&
               ((marker->base_idx <= block->sqe &&
                 marker->alignment_idx == alignment_idx) || (marker->alignment_idx != alignment_idx))) {
            if (block->sqe - block_margin <= marker->base_idx &&
                marker->alignment_idx == alignment_idx) {
                qual[marker->base_idx] = 0;
            }
            j += bam_is_rev(alignment->record) ? -1 : 1;
            marker = ((j < markers_len) && (j >= 0)) ? stList_get(markers, j) : NULL;
        }
    }
    ptCigarIt_destruct(cigar_it);
}

void calc_update_baq_all(const faidx_t *fai,
                         ptAlignment **alignments, int alignments_len,
                         stList *markers,
                         const sam_hdr_t *h,
                         double conf_d, double conf_e, double conf_b, int set_q) {
    const char *contig_name;
    int tid;
    for (int i = 0; i < alignments_len; i++) {
        tid = alignments[i]->record->core.tid;
        contig_name = sam_hdr_tid2name(h, tid);
        calc_local_baq(fai, contig_name, alignments[i], i, markers, conf_d, conf_e, conf_b, set_q);
    }
    ptMarker *marker;
    uint8_t *q;
    //DEBUG_PRINT("\t\t## start updating marker base qualities\n");
    for (int i = 0; i < stList_length(markers); i++) {
        marker = stList_get(markers, i);
        q = bam_get_qual(alignments[marker->alignment_idx]->record);
        marker->base_q = q[marker->base_idx];
    }
}


void print_markers(stList *markers) {
    if (stList_length(markers) == 0) DEBUG_PRINT("NO MARKERS!\n");
    DEBUG_PRINT("#alignment_idx\tseq_pos\tread_pos_f\tq\tis_match\tref_pos\n");
    for (int t = 0; t < stList_length(markers); t++) {
        ptMarker *m = stList_get(markers, t);
        DEBUG_PRINT("%d\t%d\t%d\t%d\t%d\t%d\n", m->alignment_idx, m->base_idx, m->read_pos_f, m->base_q, m->is_match,
                    m->ref_pos);
    }
}

ptBlock* ptMarker_convert_to_block(ptMarker* marker){
    ptBlock* block = ptBlock_construct(marker->ref_pos, marker->ref_pos,
                                       marker->base_idx, marker->base_idx,
                                       marker->read_pos_f, marker->read_pos_f);
    return block;
}

void ptMarker_add_marker_blocks_by_contig(stHash *blocks_per_contig, char* contig, int alignment_idx, stList* markers){
    stList *blocks = stHash_search(blocks_per_contig, contig);
    if (blocks == NULL) {
        blocks = stList_construct3(0, ptBlock_destruct);
        // contig name should be copied prior to inserting as a key
        stHash_insert(blocks_per_contig, copyString(contig), blocks);
    }
    for (int i =0; i < stList_length(markers); i++){
        ptMarker* marker = stList_get(markers, i);
        // add marker block only if its alignment index matches the given index
        if(marker->alignment_idx == alignment_idx) {
            stList_append(blocks, ptMarker_convert_to_block(marker));
        }
    }
}

