#include "ptBlock.h"


ptBlock *ptBlock_construct(int rfs, int rfe, int sqs, int sqe, int rds_f, int rde_f) {
    ptBlock *block = malloc(sizeof(ptBlock));
    block->rfs = rfs;
    block->rfe = rfe;
    block->sqs = sqs;
    block->sqe = sqe;
    block->rds_f = rds_f;
    block->rde_f = rde_f;
    block->data = NULL;
    block->destruct_data = free;
    block->copy_data = NULL;
    block->extend_data = NULL;
    return block;
}

ptBlock *ptBlock_construct_with_count(int rfs, int rfe, int sqs, int sqe, int rds_f, int rde_f, int count) {
    ptBlock *block = malloc(sizeof(ptBlock));
    block->rfs = rfs;
    block->rfe = rfe;
    block->sqs = sqs;
    block->sqe = sqe;
    block->rds_f = rds_f;
    block->rde_f = rde_f;
    block->data = NULL;
    block->destruct_data = free;
    block->copy_data = NULL;
    block->extend_data = NULL;
    int* data = malloc(sizeof(int));
    *data = count;
    ptBlock_set_data(block, data, destruct_count_data, copy_count_data, extend_count_data);
    return block;
}

int ptBlock_get_count(ptBlock* block){
    int* count = (int*) block->data;
    return *count;
}

void ptBlock_set_data(ptBlock *block, void *data, void (*destruct_data)(void *), void *(*copy_data)(void *),
                      void (*extend_data)(void *, void *)) {
    if (block->data != NULL) {
        block->destruct_data(block->data);
    }
    block->data = data;
    block->destruct_data = destruct_data;
    block->copy_data = copy_data;
    block->extend_data = extend_data;
}

void ptBlock_extend_data(ptBlock *block, void *data) {
    if (block->extend_data != NULL) {
        block->extend_data(block->data, data);
    }
}

void ptBlock_destruct_data(ptBlock *block) {
    if (block->destruct_data != NULL) {
        block->destruct_data(block->data);
    }
}

void *ptBlock_copy_data(ptBlock *block) {
    if (block->copy_data != NULL) {
        return block->copy_data(block->data);
    } else {
        return NULL;
    }
}


void extend_count_data(void *dest_, void *src_){
    int* dest = dest_;
    int* src = src_;
    *dest += *src;
}


void destruct_count_data(void* src){
    free(src);
}


void *copy_count_data(void* src_){
    int* src = src_;
    int* dest = malloc(sizeof(int));
    *dest = *src;
    return dest;
}

ptBlock *ptBlock_copy(ptBlock *block) {
    ptBlock *block_copy = ptBlock_construct(block->rfs,
                                            block->rfe,
                                            block->sqs,
                                            block->sqe,
                                            block->rds_f,
                                            block->rde_f);
    void *data_copy = ptBlock_copy_data(block);
    ptBlock_set_data(block_copy, data_copy, block->destruct_data, block->copy_data, block->extend_data);
    return block_copy;
}


void ptBlock_destruct(ptBlock *block) {
    ptBlock_destruct_data(block);
    free(block);
}


/// Get Functions ////
int ptBlock_get_rfs(ptBlock *block) {
    return block->rfs;
}


int ptBlock_get_rfe(ptBlock *block) {
    return block->rfe;
}

int ptBlock_get_sqs(ptBlock *block) {
    return block->sqs;
}

int ptBlock_get_sqe(ptBlock *block) {
    return block->sqe;
}

int ptBlock_get_rds_f(ptBlock *block) {
    return block->rds_f;
}

int ptBlock_get_rde_f(ptBlock *block) {
    return block->rde_f;
}


//// Set Functions ////

void ptBlock_set_rfs(ptBlock *block, int rfs) {
    block->rfs = rfs;
}

void ptBlock_set_rfe(ptBlock *block, int rfe) {
    block->rfe = rfe;
}

void ptBlock_set_sqs(ptBlock *block, int sqs) {
    block->sqs = sqs;
}

void ptBlock_set_sqe(ptBlock *block, int sqe) {
    block->sqe = sqe;
}

void ptBlock_set_rds_f(ptBlock *block, int rds_f) {
    block->rds_f = rds_f;
}

void ptBlock_set_rde_f(ptBlock *block, int rde_f) {
    block->rde_f = rde_f;
}

//// Comparing Functions ////
int ptBlock_cmp_rfs(const void *a, const void *b) {
    ptBlock *b1 = (ptBlock *) a;
    ptBlock *b2 = (ptBlock *) b;
    return b1->rfs - b2->rfs;
}

int ptBlock_cmp_rds_f(const void *a, const void *b) {
    ptBlock *b1 = (ptBlock *) a;
    ptBlock *b2 = (ptBlock *) b;
    return b1->rds_f - b2->rds_f;
}

int ptBlock_cmp_sqs(const void *a, const void *b) {
    ptBlock *b1 = (ptBlock *) a;
    ptBlock *b2 = (ptBlock *) b;
    return b1->sqs - b2->sqs;
}


bool ptBlock_is_equal(ptBlock* block_1, ptBlock* block_2){

    if((block_1 == NULL) && (block_2 == NULL)){ // if both are NULL
        return true;
    }else if((block_1 == NULL) || (block_2 == NULL)){ // if only one block is NULL
        return false;
    }

    bool is_equal = true;
    is_equal &= (block_1->rfs == block_2->rfs);
    is_equal &= (block_1->rfe == block_2->rfe);
    is_equal &= (block_1->sqs == block_2->sqs);
    is_equal &= (block_1->sqe == block_2->sqe);
    is_equal &= (block_1->rds_f == block_2->rds_f);
    is_equal &= (block_1->rde_f == block_2->rde_f);
    return is_equal;
}


bool ptBlock_is_equal_stList(stList* blocks_1, stList* blocks_2){
    if((blocks_1 == NULL) && (blocks_2 == NULL)){ // if both are NULL
        return true;
    }else if((blocks_1 == NULL) || (blocks_2 == NULL)){ // if only one block list is NULL
        return false;
    }
    if(stList_length(blocks_1) != stList_length(blocks_2)) return false;
    for(int i=0; i < stList_length(blocks_1); i++){
        ptBlock* block_1 = stList_get(blocks_1, i);
        ptBlock* block_2 = stList_get(blocks_2, i);
        if(ptBlock_is_equal(block_1, block_2) == false) return false;
    }
    return true;
}


bool ptBlock_is_equal_stHash(stHash* blocks_per_contig_1, stHash* blocks_per_contig_2){
    bool is_equal = true;
    char* contig_name;
    stList* blocks_1;
    stList* blocks_2;

    // check if all blocks in table 1 is in table 2
    stHashIterator *it_1 = stHash_getIterator(blocks_per_contig_1);
    while ((contig_name = stHash_getNext(it_1)) != NULL) {
        blocks_1 = stHash_search(blocks_per_contig_1, contig_name);
        blocks_2 = stHash_search(blocks_per_contig_2, contig_name);
        is_equal &= ptBlock_is_equal_stList(blocks_1, blocks_2);
    }
    // check if all blocks in table 2 is in table 1
    stHashIterator *it_2 = stHash_getIterator(blocks_per_contig_2);
    while ((contig_name = stHash_getNext(it_2)) != NULL) {
        blocks_1 = stHash_search(blocks_per_contig_1, contig_name);
        blocks_2 = stHash_search(blocks_per_contig_2, contig_name);
        is_equal &= ptBlock_is_equal_stList(blocks_1, blocks_2);
    }
    stHash_destructIterator(it_1);
    stHash_destructIterator(it_2);
    return is_equal;
}


// Functions for block iterator


ptBlockItrPerContig *ptBlockItrPerContig_construct(stHash *blocks_per_contig){
    ptBlockItrPerContig *block_iter = malloc(sizeof(ptBlockItrPerContig));
    block_iter->blocks_per_contig = blocks_per_contig;
    block_iter->ctg_list = stList_construct3(0, free);
    stHashIterator *it = stHash_getIterator(blocks_per_contig);
    char* contig_name;
    while ((contig_name = stHash_getNext(it)) != NULL) {
        stList_append(block_iter->ctg_list, copyString(contig_name));
    }
    stList_sort(block_iter->ctg_list, (int (*)(const void *, const void *))strcmp);
    block_iter->ctg_index = 0;
    block_iter->block_index = 0;
    return block_iter;
}


ptBlock* ptBlockItrPerContig_next(ptBlockItrPerContig *block_iter, char* ctg_name){
    ptBlock* block;
    if(stList_length(block_iter->ctg_list) - 1 < block_iter->ctg_index){
        // if all blocks are iterated already
        // set contig_name to null and also return null as the next block
        ctg_name[0] = '\0';
        block = NULL;
    }else{
        char* ctg_name_new = stList_get(block_iter->ctg_list, block_iter->ctg_index);
        // update contig name
        strcpy(ctg_name, ctg_name_new);
        stList* blocks = stHash_search(block_iter->blocks_per_contig, ctg_name);
        if(stList_length(blocks) - 1 < block_iter->block_index) {
            // if the blocks for this contig is finished
            // go to the next contig and reset the block index
            // call "ptBlockItrPerContig_next" recursively
            block_iter->ctg_index += 1;
            block_iter->block_index = 0;
            block = ptBlockItrPerContig_next(block_iter, ctg_name);
        }else{
            // if there is still block for the current contig
            // get the block and increase the block index
            block = stList_get(blocks, block_iter->block_index);
            block_iter->block_index += 1;
        }
    }
    return block;
}

void ptBlockItrPerContig_destruct(ptBlockItrPerContig *blockItr){
    blockItr->blocks_per_contig = NULL;
    stList_destruct(blockItr->ctg_list);
    free(blockItr);
}

void ptBlock_add_block_to_stList_table(stHash* blocks_per_contig, ptBlock* block, char* ctg_name){
    stList *blocks = stHash_search(blocks_per_contig, ctg_name);
    if (blocks == NULL) {
        blocks = stList_construct3(0, ptBlock_destruct);
        stHash_insert(blocks_per_contig, copyString(ctg_name), blocks);
    }
    stList_append(blocks, block);
}

stList* ptBlock_split_into_batches(stHash *blocks_per_contig, int split_number){
    stList* batches = stList_construct3(0, stHash_destruct);
    int total_size = ptBlock_get_total_length_by_rf(blocks_per_contig);
    // ceil is used here so the last batch might be smaller than other batches
    int batch_size = ceil((double)total_size / split_number);
    int batch_filled_size = 0;
    // make a block iterator
    ptBlockItrPerContig *block_iter = ptBlockItrPerContig_construct(blocks_per_contig);
    char ctg_name[200];
    // get the first block 
    ptBlock *block = ptBlockItrPerContig_next(block_iter, ctg_name);
    int start = block->rfs;
    for (int i =0 ;  i < split_number; i++) {
	// start filling the new batch
        batch_filled_size = 0;
        stHash *batch_blocks_per_contig = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free,
                                                            (void (*)(void *)) stList_destruct);
        while ((start + batch_size - batch_filled_size - 1) >= block->rfe) {
            // create block and add to the batch table
            ptBlock *block_to_add = ptBlock_construct(start, block->rfe,
                                                      -1, -1,
                                                      -1, -1);
            ptBlock_add_block_to_stList_table(batch_blocks_per_contig, block_to_add, ctg_name);
            batch_filled_size += block->rfe - start + 1;
            // go to next block
            block = ptBlockItrPerContig_next(block_iter, ctg_name);
            if (block == NULL) break;
	    // a new block is taken from iterator so start should be equal to block->rfs
            start = block->rfs;
        }
        if ((block != NULL) && (batch_filled_size < batch_size)) {
	    // only a part of the current block will be in the current batch
            ptBlock *block_to_add = ptBlock_construct(start, start + batch_size - batch_filled_size - 1,
                                                      -1, -1,
                                                      -1, -1);
            ptBlock_add_block_to_stList_table(batch_blocks_per_contig, block_to_add, ctg_name);
	    // update start location of this block for the next batch
            start += batch_size - batch_filled_size;
        }
        // add the new batch
        stList_append(batches, batch_blocks_per_contig);
        if (block == NULL) break;
    }
    ptBlockItrPerContig_destruct(block_iter);
    return batches;
}


stHash *ptBlock_parse_bed(char *bed_path) {
    FILE *fp = fopen(bed_path, "r");
    size_t read;
    size_t len;
    char *line = NULL;
    char *contig_name;
    int start;
    int end;
    char *token;
    ptBlock *block = NULL;
    stList *blocks = NULL;
    stHash *blocks_per_contig = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL,
                                                  (void (*)(void *)) stList_destruct);
    while ((read = getline(&line, &len, fp)) != -1) {
        // replace '\n' by '\0'
        if (line[strlen(line) - 1] == '\n') {
            line[strlen(line) - 1] = '\0';
        }
        token = strtok(line, "\t");
        contig_name = copyString(token);
        token = strtok(NULL, "\t");
        start = atoi(token); // 0-based
        token = strtok(NULL, "\t");
        end = atoi(token) - 1; // 0-based
        block = ptBlock_construct(start, end, -1, -1, -1, -1);
        blocks = stHash_search(blocks_per_contig, contig_name);
        // if contig does not exist as a key in the table
        // add contig along with an empty list as its value
        if (blocks == NULL) {
            blocks = stList_construct3(0, ptBlock_destruct);
            stHash_insert(blocks_per_contig, contig_name, blocks);
        }
        stList_append(blocks, block);
    }
    // sort blocks per contig
    stHashIterator *it = stHash_getIterator(blocks_per_contig);
    while ((contig_name = stHash_getNext(it)) != NULL) {
        blocks = stHash_search(blocks_per_contig, contig_name);
        stList_sort(blocks, ptBlock_cmp_rfs);
    }
    stHash_destructIterator(it);
    return blocks_per_contig;
}

void ptBlock_print_blocks_stHash(stHash* blocks_per_contig, bool print_count, FILE* fp){
    char* ctg_name;
    stHashIterator *it = stHash_getIterator(blocks_per_contig);
    while ((ctg_name = stHash_getNext(it)) != NULL) {
        stList* blocks = stHash_search(blocks_per_contig, ctg_name);
        for(int i=0; i < stList_length(blocks); i++){
            ptBlock* block = stList_get(blocks, i);
            if(print_count){
                int *count_ptr = block->data;
                fprintf(fp, "%s\t%d\t%d\t%d\n",ctg_name, block->rfs, block->rfe + 1, *count_ptr);
            }else {
                fprintf(fp, "%s\t%d\t%d\n", ctg_name, block->rfs, block->rfe + 1);
            }
        }
    }
    stHash_destructIterator(it);
}



void ptBlock_sort_stHash_by_rfs(stHash *blocks_per_contig) {
    char *contig_name;
    stList *blocks;
    stHashIterator *it = stHash_getIterator(blocks_per_contig);
    while ((contig_name = stHash_getNext(it)) != NULL) {
        blocks = stHash_search(blocks_per_contig, contig_name);
        stList_sort(blocks, ptBlock_cmp_rfs);
    }
}

stList *ptBlock_merge_blocks(stList *blocks,
                             int (*get_start)(ptBlock *), int (*get_end)(ptBlock *),
                             void (*set_end)(ptBlock *, int)) {
    stList *blocks_merged = stList_construct3(0, ptBlock_destruct);
    if (stList_length(blocks) == 0) return blocks_merged;
    ptBlock *b = NULL;
    ptBlock *b_merged = NULL;
    for (int i = 0; i < stList_length(blocks); i++) {
        b = stList_get(blocks, i);
        //printf("%d\t%d\n", b->rds_f, b->rde_f);
        if (i == 0) { // Initiate b_merged for the first block
            b_merged = ptBlock_copy(b);
            continue;
        }
        if (get_end(b_merged) < get_start(b)) {// no overlap with previous merged block
            //save the merged block
            stList_append(blocks_merged, b_merged);
            // Initiate a new merged block
            b_merged = ptBlock_copy(b);
        } else { //there is overlap
            set_end(b_merged, max(get_end(b_merged), get_end(b))); // extend end pos of the merged block
            if (b->data != NULL) {
                // add the data of the new block to the merged block
                ptBlock_extend_data(b_merged, b->data);
            }
        }
    }

    // Add the last merged block
    if (b_merged) {
        stList_append(blocks_merged, b_merged);
    }

    return blocks_merged;
}

stList *ptBlock_merge_blocks_v2(stList *blocks,
                             int (*get_start)(ptBlock *), int (*get_end)(ptBlock *),
                             void (*set_start)(ptBlock *, int), void (*set_end)(ptBlock *, int)) {
    stList *blocks_merged_finalized = stList_construct3(0, ptBlock_destruct);
    stList *blocks_merged_ongoing = stList_construct3(0, ptBlock_destruct);
    stList *blocks_merged_temp;
    if (stList_length(blocks) == 0) return blocks_merged_finalized;
    ptBlock *b2 = NULL;
    ptBlock *b1 = NULL;
    ptBlock *b_merged = NULL;
    int s1;
    int e1;
    int s2;
    int e2;
    for (int i = 0; i < stList_length(blocks); i++) {
        b2 = stList_get(blocks, i);
        //printf("%d\t%d\n", b->rds_f, b->rde_f);
        if (stList_length(blocks_merged_ongoing) == 0) { // Initiate b_merged for the first block
            b_merged = ptBlock_copy(b2);
            stList_append(blocks_merged_ongoing, b_merged);
            continue;
        }
        e2 = get_end(b2);
        s2 = get_start(b2);
        blocks_merged_temp = blocks_merged_ongoing;
        blocks_merged_ongoing = stList_construct3(0, ptBlock_destruct);
        for (int j = 0; j < stList_length(blocks_merged_temp); j++) {
            b1 = stList_get(blocks_merged_temp, j);
            e1 = get_end(b1);
            s1 = get_start(b1);
            /*
             * finalized:
             *
             *  s1       e1
             * [**********]
             *               [----------]
             *                s2       e2
             */
            if(e1 < s2) {
                b_merged = ptBlock_copy(b1);
                stList_append(blocks_merged_finalized, b_merged);
            }
            else if (s1 <= s2){
                /*
                 * finalized:
                 *   s1       e1
                 *  [***-------]
                 *     [----------]
                 *      s2       e2
                 */
                if (s1 < s2){ // && s2 <= e1
                    b_merged = ptBlock_copy(b1);
                    set_end(b_merged, s2 - 1);
                    stList_append(blocks_merged_finalized, b_merged);
                }
                /*
                 *  ongoing:
                 *
                 *    s1       e1                      s1        e1
                 *   [---*******]           OR        [----*****--]
                 *      [*******--]                       [*****]
                 *       s2      e2                        s2   e2
                 *
                 */
                b_merged = ptBlock_copy(b1);
                set_start(b_merged, s2);
                set_end(b_merged, min(e1, e2));
                if (b2->data != NULL) {
                    // add the data of the new block to the merged block
                    ptBlock_extend_data(b_merged, b2->data);
                }
                stList_append(blocks_merged_ongoing, b_merged);

                /*
                 * finalized:
                 *     s1       e1
                 *    [-----*****]
                 *      [---]
                 *       s2 e2
                 */
                if (e2 < e1){
                    b_merged = ptBlock_copy(b1);
                    set_start(b_merged, e2 + 1);
                    stList_append(blocks_merged_ongoing, b_merged);
                }
            }
            /*
             * ongoing:
             *
             *            s1       e1
             *           [**********]
             *       [----**********---]
             *        s2              e2
             */
            else if(e1 <= e2){ // && s2 < s1
                b_merged = ptBlock_copy(b1);
                if (b2->data != NULL) {
                    // add the data of the new block to the merged block
                    ptBlock_extend_data(b_merged, b2->data);
                }
                stList_append(blocks_merged_ongoing, b_merged);
            }
            else { // e2 < e1 && s2 < s1
                /*
                 * ongoing:
                 *
                 *            s1       e1
                 *           [******----]
                 *       [----******]
                 *        s2       e2
                 */
                if(s1 <= e2) {
                    b_merged = ptBlock_copy(b1);
                    set_end(b_merged, e2);
                    if (b2->data != NULL) {
                        // add the data of the new block to the merged block
                        ptBlock_extend_data(b_merged, b2->data);
                    }
                    stList_append(blocks_merged_ongoing, b_merged);
                }
                /*
                * ongoing:
                *
                *            s1        e1
                *           [------****]
                *       [----------]
                *        s2       e2
                */

                b_merged = ptBlock_copy(b1);
                set_start(b_merged, max(e2 + 1, s1));
                stList_append(blocks_merged_ongoing, b_merged);
             }
        }
        // add the last non-overlapping block
        if (max(e1 + 1, s2) <= e2){
            b_merged = ptBlock_copy(b2);
            set_start(b_merged, max(e1 + 1, s2));
            stList_append(blocks_merged_ongoing, b_merged);
        }
        // destruct the temporary blocks
        stList_destruct(blocks_merged_temp);
    }

    // Add the remaining blocks
    for (int j = 0; j < stList_length(blocks_merged_ongoing); j++) {
        b1 = stList_get(blocks_merged_ongoing, j);
        b_merged = ptBlock_copy(b1);
        stList_append(blocks_merged_finalized, b_merged);
    }

    stList_destruct(blocks_merged_ongoing);

    return blocks_merged_finalized;
}

stHash *ptBlock_merge_blocks_per_contig(stHash *blocks_per_contig,
                                        int (*get_start)(ptBlock *), int (*get_end)(ptBlock *),
                                        void (*set_end)(ptBlock *, int)) {
    char *contig_name;
    stList *blocks;
    stList *merged_blocks;
    stHash *merged_blocks_per_contig = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL,
                                                         (void (*)(void *)) stList_destruct);
    stHashIterator *it = stHash_getIterator(blocks_per_contig);
    while ((contig_name = stHash_getNext(it)) != NULL) {
        // get blocks
        blocks = stHash_search(blocks_per_contig, contig_name);
        // merge blocks
        merged_blocks = ptBlock_merge_blocks(blocks, get_start, get_end, set_end);
        // add merged blocks to the new table
        stHash_insert(merged_blocks_per_contig, contig_name, merged_blocks);
    }
    return merged_blocks_per_contig;
}


stHash *ptBlock_merge_blocks_per_contig_by_rf(stHash *blocks_per_contig) {
    return ptBlock_merge_blocks_per_contig(blocks_per_contig, ptBlock_get_rfs,
                                              ptBlock_get_rfe,
                                              ptBlock_set_rfe);
}

stHash *ptBlock_merge_blocks_per_contig_by_rd_f(stHash *blocks_per_contig) {
    return ptBlock_merge_blocks_per_contig(blocks_per_contig, ptBlock_get_rds_f,
                                              ptBlock_get_rde_f,
                                              ptBlock_set_rde_f);
}

stHash *ptBlock_merge_blocks_per_contig_by_sq(stHash *blocks_per_contig) {
    return ptBlock_merge_blocks_per_contig(blocks_per_contig, ptBlock_get_sqs,
                                           ptBlock_get_sqe,
                                           ptBlock_set_sqe);
}

stHash *ptBlock_merge_blocks_per_contig_v2(stHash *blocks_per_contig,
                                        int (*get_start)(ptBlock *), int (*get_end)(ptBlock *),
                                           void (*set_start)(ptBlock *, int), void (*set_end)(ptBlock *, int)) {
    char *contig_name;
    stList *blocks;
    stList *merged_blocks;
    stHash *merged_blocks_per_contig = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL,
                                                         (void (*)(void *)) stList_destruct);
    stHashIterator *it = stHash_getIterator(blocks_per_contig);
    while ((contig_name = stHash_getNext(it)) != NULL) {
        // get blocks
        blocks = stHash_search(blocks_per_contig, contig_name);
        // merge blocks
        merged_blocks = ptBlock_merge_blocks_v2(blocks, get_start, get_end, set_start, set_end);
        // add merged blocks to the new table
        stHash_insert(merged_blocks_per_contig, contig_name, merged_blocks);
    }
    return merged_blocks_per_contig;
}

stHash *ptBlock_merge_blocks_per_contig_by_rf_v2(stHash *blocks_per_contig) {
    return ptBlock_merge_blocks_per_contig_v2(blocks_per_contig, ptBlock_get_rfs,
                                           ptBlock_get_rfe,
                                           ptBlock_set_rfs,
                                           ptBlock_set_rfe);
}

stHash *ptBlock_merge_blocks_per_contig_by_rd_f_v2(stHash *blocks_per_contig) {
    return ptBlock_merge_blocks_per_contig_v2(blocks_per_contig, ptBlock_get_rds_f,
                                           ptBlock_get_rde_f,
                                           ptBlock_set_rds_f,
                                           ptBlock_set_rde_f);
}

stHash *ptBlock_merge_blocks_per_contig_by_sq_v2(stHash *blocks_per_contig) {
    return ptBlock_merge_blocks_per_contig_v2(blocks_per_contig, ptBlock_get_sqs,
                                           ptBlock_get_sqe,
                                           ptBlock_set_sqs,
                                           ptBlock_set_sqe);
}

int ptBlock_get_total_number(stHash *blocks_per_contig) {
    int n = 0;
    char *contig_name;
    stList *blocks;
    stHashIterator *it = stHash_getIterator(blocks_per_contig);
    while ((contig_name = stHash_getNext(it)) != NULL) {
        blocks = stHash_search(blocks_per_contig, contig_name);
        n += stList_length(blocks);
    }
    return n;
}

int ptBlock_get_total_length(stHash *blocks_per_contig, int (*get_start)(ptBlock *), int (*get_end)(ptBlock *)) {
    int total_len = 0;
    char *contig_name;
    stList *blocks;
    ptBlock *block;
    stHashIterator *it = stHash_getIterator(blocks_per_contig);
    while ((contig_name = stHash_getNext(it)) != NULL) {
        blocks = stHash_search(blocks_per_contig, contig_name);
        for (int i = 0; i < stList_length(blocks); i++) {
            block = stList_get(blocks, i);
            total_len += get_end(block) - get_start(block) + 1;
        }
    }
    return total_len;
}

int ptBlock_get_total_length_by_rf(stHash *blocks_per_contig) {
    return ptBlock_get_total_length(blocks_per_contig, ptBlock_get_rfs, ptBlock_get_rfe);
}

int ptBlock_get_total_length_by_rd_f(stHash *blocks_per_contig) {
    return ptBlock_get_total_length(blocks_per_contig, ptBlock_get_rds_f, ptBlock_get_rde_f);
}

int ptBlock_get_total_length_by_sq(stHash *blocks_per_contig) {
    return ptBlock_get_total_length(blocks_per_contig, ptBlock_get_sqs, ptBlock_get_sqe);
}


void ptBlock_add_alignment(stHash *blocks_per_contig, ptAlignment *alignment, bool init_count_data) {
    ptBlock *block = ptBlock_construct(alignment->rfs,
                                       alignment->rfe,
                                       -1, -1,
                                       -1, -1);
    if (init_count_data){
        int* count_data = malloc(sizeof(int));
        *count_data = 1;
        ptBlock_set_data(block, count_data,
                         destruct_count_data,
                         copy_count_data,
                         extend_count_data);
    }
    stList *blocks = stHash_search(blocks_per_contig, alignment->contig);
    if (blocks == NULL) {
        blocks = stList_construct3(0, ptBlock_destruct);
        // contig name should be copied prior to inserting as a key
        stHash_insert(blocks_per_contig, copyString(alignment->contig), blocks);
    }
    stList_append(blocks, block);
}

void ptBlock_save_in_bed(stHash *blocks_per_contig, char* bed_path, bool print_count_data){
    // get contigs as a list and sort them
    stList * contigs = stHash_getKeys(blocks_per_contig);
    stList_sort(contigs, (int (*)(const void *, const void *))strcmp);

    // open bed to write
    FILE *fp = fopen(bed_path, "w");
    if (fp == NULL) {
        fprintf(stderr, "[%s] Error: Failed to open file %s.\n", get_timestamp(), bed_path);
    }
    // iterate over sorted contig names
    for(int i = 0; i < stList_length(contigs); i++){
        char* contig = stList_get(contigs, i);
        stList *blocks = stHash_search(blocks_per_contig, contig);
        // iterate over blocks in this contig
        for(int j = 0; j < stList_length(blocks); j++) {
            ptBlock* block = stList_get(blocks, j);
            // it may happen when the whole projected variant block is within an insertion in another haplotype
            // TODO: experiment if it can help to output these coordinates
            if(block->rfe < block->rfs) continue;
            if(print_count_data) {
                fprintf(fp, "%s\t%d\t%d\t%d\n", contig, block->rfs, block->rfe + 1, *((int*)block->data)); // end should be 1-based in BED
            }else{
                fprintf(fp, "%s\t%d\t%d\n", contig, block->rfs, block->rfe + 1); // end should be 1-based in BED
            }
        }

    }
    fclose(fp);
}

void ptBlock_add_blocks_by_contig(stHash *blocks_per_contig, char* contig, stList *blocks_to_add) {
    stList *blocks = stHash_search(blocks_per_contig, contig);
    if (blocks == NULL) {
        blocks = stList_construct3(0, ptBlock_destruct);
        // contig name should be copied prior to inserting as a key
        stHash_insert(blocks_per_contig, copyString(contig), blocks);
    }
    for (int i =0; i < stList_length(blocks_to_add); i++){
        stList_append(blocks, ptBlock_copy(stList_get(blocks_to_add,i)));
    }
}

stList* ptBlock_copy_stList(stList* blocks) {
    stList *copy_blocks = stList_construct3(0, ptBlock_destruct);
    for (int i = 0; i < stList_length(blocks); i++) {
        ptBlock* block = stList_get(blocks, i);
        stList_append(copy_blocks, ptBlock_copy(block));
    }
    return copy_blocks;
}

