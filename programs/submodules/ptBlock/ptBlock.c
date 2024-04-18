#include "ptBlock.h"
#include "stdlib.h"
#include "stdio.h"
#include "track_reader.h"
#include <zlib.h>


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

char *get_string_count_data(void* src_){
    int* src = src_;
    char* str = malloc(10);
    sprintf(str, "%d", *src);
    return str;
}

CoverageInfo *CoverageInfo_construct(int32_t annotation_flag,
                            u_int16_t coverage,
                            u_int16_t coverage_high_mapq,
                            u_int16_t coverage_high_clip){
    CoverageInfo *cov_info = malloc(sizeof(CoverageInfo));
    cov_info->annotation_flag = annotation_flag;
    cov_info->coverage = coverage;
    cov_info->coverage_high_mapq = coverage_high_mapq;
    cov_info->coverage_high_clip = coverage_high_clip;
    return cov_info;
}

CoverageInfo *CoverageInfo_copy(CoverageInfo *src){
    return CoverageInfo_construct(src->annotation_flag,
                                  src->coverage,
                                  src->coverage_high_mapq,
                                  src->coverage_high_clip);
}


void CoverageInfo_reset(CoverageInfo *coverageInfo){
    coverageInfo->annotation_flag = 0;
    coverageInfo->coverage = 0;
    coverageInfo->coverage_high_mapq = 0;
    coverageInfo->coverage_high_clip = 0;
}


CoverageInfo **CoverageInfo_construct1DArray(int len){
    CoverageInfo **array = (CoverageInfo **) malloc(len * sizeof(CoverageInfo *));
    for(int i=0; i < len; i++){
        array[i] = CoverageInfo_construct(0,0,0,0);
    }
    return array;
}


CoverageInfo **CoverageInfo_copy1DArray(CoverageInfo **src, int len){
    CoverageInfo **dest = (CoverageInfo **) malloc(len * sizeof(CoverageInfo *));
    for(int i=0; i < len; i++){
        dest[i] = CoverageInfo_copy(src[i]);
    }
    return dest;
}

void CoverageInfo_destruct1DArray(CoverageInfo **coverageInfo1DArray, int len){
    for(int i=0; i < len; i++){
        CoverageInfo_destruct(coverageInfo1DArray[i]);
    }
    free(coverageInfo1DArray);
}

void CoverageInfo_destruct(CoverageInfo *coverageInfo){
    free(coverageInfo);
}

int32_t CoverageInfo_getRegionBitRepresentation(int regionIndex){
	return 1 << regionIndex;
}

int CoverageInfo_getRegionIndex(CoverageInfo *coverageInfo){
	return getFirstIndexWithNonZeroBitFromRight(coverageInfo->annotation_flag);
}


u_int16_t CoverageInfo_getCoverage(CoverageInfo *coverageInfo){
	return coverageInfo->coverage;
}
u_int16_t CoverageInfo_getCoverageHighMapq(CoverageInfo *coverageInfo){
	return coverageInfo->coverage_high_mapq;
}
u_int16_t CoverageInfo_getCoverageHighClip(CoverageInfo *coverageInfo){
	return coverageInfo->coverage_high_clip;
}

CoverageInfo *CoverageInfo_construct_from_alignment(ptAlignment *alignment, int min_mapq, double min_clipping_ratio){
    u_int16_t coverage_high_mapq = min_mapq <= alignment->mapq ? 1 : 0;
    int max_clip = max(alignment->r_clip, alignment->l_clip);
    int alignment_len = alignment->rfe - alignment->rfs + 1;
    u_int16_t coverage_high_clip = min_clipping_ratio <= ((double) max_clip / alignment_len) ? 1 : 0;
    return CoverageInfo_construct(0, 1, coverage_high_mapq, coverage_high_clip);
}

void extend_cov_info_data(void *dest_, void *src_){
    CoverageInfo * dest = dest_;
    CoverageInfo * src = src_;
    dest->annotation_flag |= src->annotation_flag;
    dest->coverage += src->coverage;
    dest->coverage_high_mapq += src->coverage_high_mapq;
    dest->coverage_high_clip += src->coverage_high_clip;
}


void destruct_cov_info_data(void* src){
    free(src);
}


void *copy_cov_info_data(void* src_){
    CoverageInfo * src = src_;
    CoverageInfo * dest = malloc(sizeof(CoverageInfo));
    dest->annotation_flag = src->annotation_flag;
    dest->coverage = src->coverage;
    dest->coverage_high_mapq = src->coverage_high_mapq;
    dest->coverage_high_clip = src->coverage_high_clip;
    return dest;
}

char *get_string_cov_info_data_format_1(void* src_){
    CoverageInfo * src = src_;
    char *str = malloc(150);
    sprintf(str,
            "annot=%d, cov=%d, cov_mapq=%d, cov_clip=%d",
            CoverageInfo_getRegionIndex(src),
            src->coverage,
            src->coverage_high_mapq,
            src->coverage_high_clip);
    return str;
}

char *get_string_cov_info_data_format_2(void* src_){
    CoverageInfo * src = src_;
    char *str = malloc(150);
    sprintf(str,
            "%d\t%d\t%d\t%d",
            src->coverage,
            src->coverage_high_mapq,
            src->coverage_high_clip,
            CoverageInfo_getRegionIndex(src));
    return str;
}

char *get_string_cov_info_data_format_only_total(void* src_){
    CoverageInfo * src = src_;
    char *str = malloc(50);
    sprintf(str,
            "%d",
            src->coverage);
    return str;
}

char *get_string_cov_info_data_format_only_high_mapq(void* src_){
    CoverageInfo * src = src_;
    char *str = malloc(50);
    sprintf(str,
            "%d",
            src->coverage_high_mapq);
    return str;
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
    int64_t total_size = ptBlock_get_total_length_by_rf(blocks_per_contig);
    // ceil is used here so the last batch might be smaller than other batches
    int64_t batch_size = ceil((double)total_size / split_number);
    int64_t batch_filled_size = 0;
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

void ptBlock_print_blocks_stHash_in_bed(stHash* blocks_per_contig,
                                        char * (*get_string_function)(void *),
                                        void* file_ptr,
                                        bool is_compressed){
    char* ctg_name;
    char line[1000];
    stList *sorted_contig_list = ptBlock_get_sorted_contig_list(blocks_per_contig);
    for(int ctg_i=0; ctg_i < stList_length(sorted_contig_list); ctg_i++){
        ctg_name = stList_get(sorted_contig_list, ctg_i);
	stList* blocks = stHash_search(blocks_per_contig, ctg_name);
    	for(int i=0; i < stList_length(blocks); i++){
            ptBlock* block = stList_get(blocks, i);
            if(get_string_function != NULL){
                sprintf(line,
                        "%s\t%d\t%d\t%s\n",
                        ctg_name,
                        block->rfs,
                        block->rfe + 1,
                        get_string_function((void *) block->data));
            }else {
                sprintf(line,
                        "%s\t%d\t%d\n",
                        ctg_name,
                        block->rfs,
                        block->rfe + 1);
            }
            if(is_compressed){
                gzFile* gzFile_ptr = file_ptr;
                gzprintf(*gzFile_ptr, line);
            }else{
                FILE* fp = file_ptr;
                fprintf(fp,line);
            }
        }
    }
    stList_destruct(sorted_contig_list);
}


stList *ptBlock_get_sorted_contig_list(stHash* blocks_per_contig){
	char* ctg_name;
	stHashIterator *it = stHash_getIterator(blocks_per_contig);
	stList* contig_list = stList_construct3(0, free);
	while ((ctg_name = stHash_getNext(it)) != NULL) {
		stList_append(contig_list, copyString(ctg_name));
	}
	// sort contig list
	stList_sort(contig_list, strcmp);
	return contig_list;
}

void ptBlock_print_blocks_stHash_in_cov(stHash* blocks_per_contig,
                                        char * (*get_string_function)(void *),
                                        void* file_ptr,
                                        bool is_compressed,
                                        stHash* ctg_to_len){

    char* ctg_name;
    char line[1000];
    stList *sorted_contig_list = ptBlock_get_sorted_contig_list(blocks_per_contig);
    for(int ctg_i=0; ctg_i < stList_length(sorted_contig_list); ctg_i++){
	ctg_name = stList_get(sorted_contig_list, ctg_i);
        // get the contig length and print it beside the contig header
        int *ctg_len_ptr = stHash_search(ctg_to_len, ctg_name);
        if(ctg_len_ptr == NULL) {
            fprintf(stderr, "[%s] Error: contig %s is not present in the bam/sam header or fai file\n", get_timestamp(), ctg_name);
            exit(EXIT_FAILURE);
        }
        sprintf(line,
                ">%s %d\n",
                ctg_name,
                *ctg_len_ptr);
        // print contig header
        if(is_compressed){
            gzFile* gzFile_ptr = file_ptr;
            gzprintf(*gzFile_ptr, line);
        }else{
            FILE* fp = file_ptr;
            fprintf(fp,line);
        }
        // iterate over the blocks in this contig and write them
        stList* blocks = stHash_search(blocks_per_contig, ctg_name);
        for(int i=0; i < stList_length(blocks); i++){
            ptBlock* block = stList_get(blocks, i);
            if(get_string_function != NULL){
                sprintf(line,
                        "%d\t%d\t%s\n",
                        block->rfs + 1, // start is 1-based in cov format
                        block->rfe + 1,
                        get_string_function((void *) block->data));
            }else { // warning: not having coverage is not meaningful when we want to write in cov format
                sprintf(line,
                        "%d\t%d\n",
                        block->rfs + 1, // start is 1-based in cov format
                        block->rfe + 1);
            }
            if(is_compressed){
                gzFile* gzFile_ptr = file_ptr;
                gzprintf(*gzFile_ptr, line);
            }else{
                FILE* fp = file_ptr;
                fprintf(fp,line);
            }
        }
    }
    stList_destruct(sorted_contig_list);
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
        stHash_insert(merged_blocks_per_contig, copyString(contig_name), merged_blocks);
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

int64_t ptBlock_get_total_number(stHash *blocks_per_contig) {
    int64_t n = 0;
    char *contig_name;
    stList *blocks;
    stHashIterator *it = stHash_getIterator(blocks_per_contig);
    while ((contig_name = stHash_getNext(it)) != NULL) {
        blocks = stHash_search(blocks_per_contig, contig_name);
        n += stList_length(blocks);
    }
    return n;
}

int64_t ptBlock_get_total_length(stHash *blocks_per_contig, int (*get_start)(ptBlock *), int (*get_end)(ptBlock *)) {
    int64_t total_len = 0;
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

int64_t ptBlock_get_total_length_by_rf(stHash *blocks_per_contig) {
    return ptBlock_get_total_length(blocks_per_contig, ptBlock_get_rfs, ptBlock_get_rfe);
}

int64_t ptBlock_get_total_length_by_rd_f(stHash *blocks_per_contig) {
    return ptBlock_get_total_length(blocks_per_contig, ptBlock_get_rds_f, ptBlock_get_rde_f);
}

int64_t ptBlock_get_total_length_by_sq(stHash *blocks_per_contig) {
    return ptBlock_get_total_length(blocks_per_contig, ptBlock_get_sqs, ptBlock_get_sqe);
}

void ptBlock_add_alignment_as_CoverageInfo(stHash *blocks_per_contig,
                                           ptAlignment *alignment,
                                           int min_mapq,
                                           double min_clipping_ratio) {
    ptBlock *block = ptBlock_construct(alignment->rfs,
                                       alignment->rfe,
                                       -1, -1,
                                       -1, -1);
    CoverageInfo * cov_info_data = CoverageInfo_construct_from_alignment(alignment,min_mapq, min_clipping_ratio);
    ptBlock_set_data(block, cov_info_data,
                     destruct_cov_info_data,
                     copy_cov_info_data,
                     extend_cov_info_data);
    stList *blocks = stHash_search(blocks_per_contig, alignment->contig);
    if (blocks == NULL) {
        blocks = stList_construct3(0, ptBlock_destruct);
        // contig name should be copied prior to inserting as a key
        stHash_insert(blocks_per_contig, copyString(alignment->contig), blocks);
    }
    stList_append(blocks, block);
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

void ptBlock_extend_block_tables(stHash *blocks_per_contig_dest, stHash *blocks_per_contig_src) {
    stHashIterator *it = stHash_getIterator(blocks_per_contig_src);
    char* contig_name;
    while ((contig_name = stHash_getNext(it)) != NULL) {
        stList* blocks_to_add = stHash_search(blocks_per_contig_src, contig_name);
        ptBlock_add_blocks_by_contig(blocks_per_contig_dest, contig_name, ptBlock_copy_stList(blocks_to_add));
    }
}

void ptBlock_add_data_to_all_blocks_stHash(stHash *blocks_per_contig,
                                           void* data,
                                           void (*destruct_data)(void *),
                                           void *(*copy_data)(void *),
                                           void (*extend_data)(void *, void *)){
    ptBlockItrPerContig * block_iter = ptBlockItrPerContig_construct(blocks_per_contig);
    char ctg_name[200];
    ptBlock* block;
    while ((block = ptBlockItrPerContig_next(block_iter, ctg_name)) != NULL) {
        ptBlock_set_data(block,
                         copy_data(data),
                         destruct_data,
                         copy_data,
                         extend_data);
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


stHash* ptBlock_get_contig_length_stHash_from_bam(char* bam_path){
    stHash *ctg_to_len = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    samFile *fp = sam_open(bam_path, "r");
    sam_hdr_t *sam_hdr = sam_hdr_read(fp);
    for(int i = 0; i < sam_hdr->n_targets; i++){
        char* ctg_name = sam_hdr->target_name[i];
        int* ctg_len_ptr = malloc(sizeof(int));
        *ctg_len_ptr = sam_hdr->target_len[i];
        stHash_insert(ctg_to_len, copyString(ctg_name), ctg_len_ptr);
    }
    sam_hdr_destroy(sam_hdr);
    sam_close(fp);
    return ctg_to_len;
}

stHash* ptBlock_get_contig_length_stHash_from_fai(char* fai_path){
    stHash *ctg_to_len = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    char *line = malloc(1000);
    FILE *fp = fopen(fai_path, "r");
    if (fp == NULL){
	    fprintf(stderr, "Error: Unable to open %s\n", fai_path);
    }
    int len;
    int read;
    char *token;
    char *ctg_name;
    while((read = getline(&line, &len, fp)) > 0){
	    //fprintf(stderr, "%s %d\n",line, read);
	    if (0 < read && line[read-1] == '\n') line[read-1] = '\0';
	    Splitter* splitter = Splitter_construct(line, '\t');
	    // get contig name
	    token = Splitter_getToken(splitter);
	    ctg_name = copyString(token);
	    //fprintf(stderr, "%s %s %s\n",ctg_name, token, line);
	    // get contig len
	    token = Splitter_getToken(splitter);
            int * ctg_len_ptr = malloc(sizeof(int));
	    *ctg_len_ptr = atoi(token);
	    // add to table
	    stHash_insert(ctg_to_len, ctg_name, ctg_len_ptr);
	    Splitter_destruct(splitter);
    }
    fclose(fp);
    free(line);
    return ctg_to_len;
}

int ptBlock_get_max_contig_length(stHash *ctg_to_len){
	char* ctg_name;
	int max_len = 0;
        stHashIterator *it = stHash_getIterator(ctg_to_len);
        while ((ctg_name = stHash_getNext(it)) != NULL) {
                int * ctg_len_ptr = stHash_search(ctg_to_len, ctg_name);
		if (max_len < *ctg_len_ptr){
			max_len = *ctg_len_ptr;
		}
        }
	stHash_destructIterator(it);
	return max_len;
}

// for one thread of parsing alignments
void _update_coverage_blocks_with_alignments(void * arg_){
    // get the arguments
    work_arg_t *arg = arg_;
    ArgumentsCovExt *argsCovExt = arg->data;
    stHash *coverage_blocks_per_contig = argsCovExt->coverage_blocks_per_contig;
    stHash *ref_blocks_per_contig_to_parse = argsCovExt->ref_blocks_per_contig_to_parse;
    char *bam_path = argsCovExt->bam_path;
    pthread_mutex_t *mutexPtr = argsCovExt->mutexPtr;
    int min_mapq = argsCovExt->min_mapq;
    double min_clipping_ratio = argsCovExt->min_clipping_ratio;

    //open bam file
    samFile *fp = sam_open(bam_path, "r");
    sam_hdr_t *sam_hdr = sam_hdr_read(fp);
    bam1_t *b = bam_init1();
    // load the bam index
    hts_idx_t * sam_idx = sam_index_load(fp, bam_path);
    char* ctg_name;

    int count_parsed_reads = 0;
    // iterate over all contigs for this batch
    stHashIterator *it = stHash_getIterator(ref_blocks_per_contig_to_parse);
    while ((ctg_name = stHash_getNext(it)) != NULL) {
        stList* ref_blocks_to_parse = stHash_search(ref_blocks_per_contig_to_parse, ctg_name);
        // get the contig id
        int tid = sam_hdr_name2tid(sam_hdr, ctg_name);
        // iterate over all blocks in this contig
        for(int i = 0; i < stList_length(ref_blocks_to_parse); i++){
            ptBlock * block = stList_get(ref_blocks_to_parse, i);
            fprintf(stderr, "[%s] Started parsing reads in the block: %s\t%d\t%d\n", get_timestamp(), ctg_name, block->rfs, block->rfe+1);
            // iterate over all alignments overlapping this block
            hts_itr_t * sam_itr = sam_itr_queryi(sam_idx, tid, block->rfs, block->rfe);
            int bytes_read;
            while (sam_itr != NULL) {
                bytes_read = sam_itr_next(fp, sam_itr, b);
                if (bytes_read <= -1)break;
                if (b->core.flag & BAM_FUNMAP) continue; // skip unmapped
                if ((b->core.flag & BAM_FSECONDARY) > 0) continue; // skip secondary alignments
                if (b->core.pos < block->rfs) continue; // make sure the alignment starts after block->rfs
                ptAlignment * alignment = ptAlignment_construct(b, sam_hdr);
                // lock the mutex, add the coverage block and unlock the mutex
                pthread_mutex_lock(mutexPtr);
                bool init_count_data = true;
                ptBlock_add_alignment_as_CoverageInfo(coverage_blocks_per_contig,
                                                      alignment,
                                                      min_mapq,
                                                      min_clipping_ratio);
		count_parsed_reads += 1;
		// log after parsing every 100k reads
		if(count_parsed_reads % 100000 == 0){
			fprintf(stderr, "[%s][%s:%d-%d] Parsed %d reads.\n", get_timestamp(), ctg_name, block->rfs, block->rfe+1, count_parsed_reads);
		}
                pthread_mutex_unlock(mutexPtr);
		ptAlignment_destruct(alignment);
            }
            if (sam_itr != NULL) hts_itr_destroy(sam_itr);
        }
    }
    free(argsCovExt);
    stHash_destructIterator(it);
    hts_idx_destroy(sam_idx);
    sam_hdr_destroy(sam_hdr);
    sam_close(fp);
    bam_destroy1(b);
}

stHash* ptBlock_get_whole_genome_blocks_per_contig(char* bam_path){
    stHash *blocks_per_contig = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free,
                                                  (void (*)(void *)) stList_destruct);
    samFile *fp = sam_open(bam_path, "r");
    sam_hdr_t *sam_hdr = sam_hdr_read(fp);
    int64_t total_len = 0;
    for(int i = 0; i < sam_hdr->n_targets; i++){
        char* ctg_name = sam_hdr->target_name[i];
        int ctg_len = sam_hdr->target_len[i];
        ptBlock * block = ptBlock_construct(0, ctg_len-1,
                                            -1, -1,
                                            -1, -1);
        stList* block_list = stList_construct3(0, ptBlock_destruct);
        stList_append(block_list, block);
        stHash_insert(blocks_per_contig, copyString(ctg_name), block_list);
        total_len += ctg_len;
    }
    fprintf("[%s] Size of the whole genome = %ld (n=%d)\n", get_timestamp(), total_len, sam_hdr->n_targets);
    sam_hdr_destroy(sam_hdr);
    sam_close(fp);
    return blocks_per_contig;
}

stHash* ptBlock_multi_threaded_coverage_extraction(char* bam_path,
                                                   int threads,
                                                   int min_mapq,
                                                   double min_clipping_ratio){
    stHash* whole_genome_blocks_per_contig = ptBlock_get_whole_genome_blocks_per_contig(bam_path);
    stList* block_batches = ptBlock_split_into_batches(whole_genome_blocks_per_contig, threads);
    stHash *coverage_blocks_per_contig = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL,
                                                           (void (*)(void *)) stList_destruct);
    // create a thread pool
    tpool_t *tm = tpool_create(threads);
    pthread_mutex_t *mutexPtr = malloc(sizeof(pthread_mutex_t));
    pthread_mutex_init(mutexPtr, NULL);
    for(int i=0; i < stList_length(block_batches); i++) {
        stHash* batch = stList_get(block_batches, i);
        // make the args struct to pass to the function that
        // is going to be run in each thread
        ArgumentsCovExt *argsCovExt = malloc(sizeof(ArgumentsCovExt));
        argsCovExt->coverage_blocks_per_contig = coverage_blocks_per_contig;
        argsCovExt->ref_blocks_per_contig_to_parse = batch;
        argsCovExt->bam_path = bam_path;
        argsCovExt->mutexPtr = mutexPtr;
        argsCovExt->min_mapq = min_mapq;
        argsCovExt->min_clipping_ratio = min_clipping_ratio;
        work_arg_t *arg = malloc(sizeof(work_arg_t));
        arg->data = (void*) argsCovExt;
        // Add a new job to the thread pool
        tpool_add_work(tm,
                       _update_coverage_blocks_with_alignments,
                       (void*) arg);
        fprintf(stderr, "[%s] Created thread for parsing batch %d (length=%ld)\n",get_timestamp(), i, ptBlock_get_total_length_by_rf(batch));
    }
    tpool_wait(tm);
    tpool_destroy(tm);

    fprintf(stderr, "[%s] All batches are parsed.\n", get_timestamp());

    pthread_mutex_destroy(mutexPtr);


    stHash_destruct(whole_genome_blocks_per_contig);
    stList_destruct(block_batches);

    fprintf(stderr, "[%s] Started sorting and merging blocks.\n", get_timestamp());
    //sort
    ptBlock_sort_stHash_by_rfs(coverage_blocks_per_contig);
    //merge
    stHash * coverage_blocks_per_contig_merged = ptBlock_merge_blocks_per_contig_by_rf_v2(coverage_blocks_per_contig);
    fprintf(stderr, "[%s] Merging coverage blocks is done.\n", get_timestamp());

    // free unmerged blocks
    stHash_destruct(coverage_blocks_per_contig);

    return coverage_blocks_per_contig_merged;
}

int get_annotation_index(stList* annotation_names, char* annotation_name){
    for(int i=0; i < stList_length(annotation_names); i++){
        if(strcmp(stList_get(annotation_names, i), annotation_name) == 0){
            return i;
        }
    }
    return -1;
}

stList *parse_annotation_names_and_save_in_stList(char* json_path){
    int buffer_size = 0;
    char* json_buffer = read_whole_file(json_path, &buffer_size, "r");
    fwrite(json_buffer,1,buffer_size,stderr);
    cJSON *annotation_json = cJSON_ParseWithLength(json_buffer, buffer_size);
    if (annotation_json == NULL)
    {
        const char *error_ptr = cJSON_GetErrorPtr();
        if (error_ptr != NULL)
        {
            fprintf(stderr, "Error before: %s\n", error_ptr);
        }
        return NULL;
    }

    int annotation_count = cJSON_GetArraySize(annotation_json);
    stList* annotation_names = stList_construct3(0, free);

    // iterate over key-values in json
    // each key is an index
    // each value is a path to a bed file
    cJSON *element = NULL;
    cJSON_ArrayForEach(element, annotation_json){
        char* annotation_name = element->string;
        stList_append(annotation_names, copyString(annotation_name));
    }
    cJSON_Delete(annotation_json);
    return annotation_names;
}

stList *parse_annotation_paths_and_save_in_stList(char* json_path){
    int buffer_size = 0;
    char* json_buffer = read_whole_file(json_path, &buffer_size, "r");
    fwrite(json_buffer,1,buffer_size,stderr);
    cJSON *annotation_json = cJSON_ParseWithLength(json_buffer, buffer_size);
    if (annotation_json == NULL)
    {
        const char *error_ptr = cJSON_GetErrorPtr();
        if (error_ptr != NULL)
        {
            fprintf(stderr, "Error before: %s\n", error_ptr);
        }
        return NULL;
    }

    stList* annotation_paths = stList_construct3(0, free);

    // iterate over key-values in json
    // each key is an index
    // each value is a path to a bed file
    cJSON *element = NULL;
    cJSON_ArrayForEach(element, annotation_json){
        char* annotation_path = cJSON_GetStringValue(element);
        stList_append(annotation_paths, copyString(annotation_path));
    } 
    cJSON_Delete(annotation_json);
    return annotation_paths;
}


stList* parse_all_annotations_and_save_in_stList(char* json_path){
    int buffer_size = 0;
    char* json_buffer = read_whole_file(json_path, &buffer_size, "r");
    fwrite(json_buffer,1,buffer_size,stderr);
    cJSON *annotation_json = cJSON_ParseWithLength(json_buffer, buffer_size);
    if (annotation_json == NULL)
    {
        const char *error_ptr = cJSON_GetErrorPtr();
        if (error_ptr != NULL)
        {
            fprintf(stderr, "Error before: %s\n", error_ptr);
        }
        return NULL;
    }

    int annotation_count = cJSON_GetArraySize(annotation_json);
    stList* block_table_list = stList_construct3(0, stHash_destruct);

    // iterate over key-values in json
    // each key is an index
    // each value is a path to a bed file
    cJSON *element = NULL;
    cJSON_ArrayForEach(element, annotation_json)
    {
        if (cJSON_IsString(element)){
            char* bed_path = cJSON_GetStringValue(element);
            stHash *annotation_block_table = ptBlock_parse_bed(bed_path);
            fprintf(stderr, "[%s] Parsed  annotation %s:%s\n", get_timestamp(), element->string, cJSON_GetStringValue(element));
            stList_append(block_table_list, annotation_block_table);
        }
    }
    cJSON_Delete(annotation_json);
    fprintf(stderr, "[%s] Number of parsed annotations = %d\n", get_timestamp(), stList_length(block_table_list));
    //for(int i=0;i < stList_length(block_table_list);i++){
    //    ptBlock_print_blocks_stHash_in_bed(stList_get(block_table_list, i), NULL, stderr, false);
    //}
    return block_table_list;

}


stHash *ptBlock_parse_coverage_info_blocks(char *filePath){
	TrackReader *trackReader = TrackReader_construct(filePath, NULL, true); //0-based coors = true
	stHash *coverage_blocks_per_contig = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL,
                                                           (void (*)(void *)) stList_destruct);
	stList *blocks = NULL;
	while(0 < TrackReader_next(trackReader)){
		// create a ptBlock based on the parsed track
		ptBlock *block = ptBlock_construct(trackReader->s, trackReader->e,
				                   -1, -1,
						   -1, -1);
		// annotation_flag is the first attribute
		int annotation_flag = CoverageInfo_getRegionBitRepresentation(atoi(trackReader->attrbs[3]));
		CoverageInfo * cov_info_data = CoverageInfo_construct(annotation_flag,
				                                      atoi(trackReader->attrbs[0]),
								      atoi(trackReader->attrbs[1]),
				                                      atoi(trackReader->attrbs[2]));
		// add coverageInfo data to block
		ptBlock_set_data(block, cov_info_data, 
				destruct_cov_info_data,
				copy_cov_info_data,
				extend_cov_info_data);

		// add block to the block stHash table
		blocks = stHash_search(coverage_blocks_per_contig, trackReader->ctg);
		// add new contig key to the table if it does not exist
		if (blocks == NULL){
			blocks = stList_construct3(0,ptBlock_destruct);
			stHash_insert(coverage_blocks_per_contig, copyString(trackReader->ctg), blocks);
		}
		stList_append(blocks, block);
	}
	TrackReader_destruct(trackReader);
	ptBlock_sort_stHash_by_rfs(coverage_blocks_per_contig);
	return coverage_blocks_per_contig;
}

// annotation block tables need to have a CoverageInfo object as their "data"
// this CoverageInfo will have zero values for the coverage related attributes
// the annotation flag is set based on the order of the annotation blocks in the
// given stList as the input
// The data augmentation is happening in place
void add_coverage_info_to_all_annotation_block_tables(stList *block_table_list){
    for(int i=0; i < stList_length(block_table_list); i++){
        stHash * block_per_contig = stList_get(block_table_list, i);
        // Each annotation has an associated flag represented by a bit-vector
        // the size of the bit-vector is 32, so it can be saved in an int32_t variable
        // for example for i=0 -> flag = 1 and for i=6 -> flag= 64
        int32_t annotation_flag = CoverageInfo_getRegionBitRepresentation(i);
        CoverageInfo * cov_info = CoverageInfo_construct(annotation_flag, 0, 0, 0);
        // The cov_info object created above will be copied and added as "data" to all blocks
        // for the current annotation. This process is happening in place
        ptBlock_add_data_to_all_blocks_stHash(block_per_contig,
                                              cov_info,
                                              destruct_cov_info_data,
                                              copy_cov_info_data,
                                              extend_cov_info_data);
    }
}

stHash *ptBlock_multi_threaded_coverage_extraction_with_zero_coverage_and_annotation(char* bam_path,
                                                                                     char* json_path,
                                                                                     int threads,
                                                                                     int min_mapq,
                                                                                     double min_clipping_ratio) {
    // parse annotation bed files
    stList *annotation_block_table_list = parse_all_annotations_and_save_in_stList(json_path);

    // add coverage info objects as data to all annotation blocks
    // each coverage info will contain only the related annotation flag with 0 coverage
    add_coverage_info_to_all_annotation_block_tables(annotation_block_table_list);

    // parse alignments and make a block table that contains the necessary coverage values per block
    // the coverage values will be related to the total alignments, alignments with high mapq and
    // each block in the output table is a maximal contiguous block with no change in the depth of coverage
    stHash *coverage_block_table = ptBlock_multi_threaded_coverage_extraction(bam_path,
                                                                              threads,
                                                                              min_mapq,
                                                                              min_clipping_ratio);
    // print len/number stats for the coverage block table
    fprintf(stderr, "[%s] Created block table with coverage data : tot_len=%ld, number=%ld\n", get_timestamp(),
            ptBlock_get_total_length_by_rf(coverage_block_table),
            ptBlock_get_total_number(coverage_block_table));

    // cover the whole genome with blocks that have zero coverage
    // this is useful to save the blocks with no coverage
    stHash *whole_genome_block_table = ptBlock_get_whole_genome_blocks_per_contig(bam_path);
    // print len/number stats for the whole genome block table
    fprintf(stderr, "[%s] Created block table for whole genome  : tot_len=%ld, number=%ld\n", get_timestamp(),
            ptBlock_get_total_length_by_rf(whole_genome_block_table),
            ptBlock_get_total_number(whole_genome_block_table));

    CoverageInfo *cov_info = CoverageInfo_construct(0, 0, 0, 0);
    // The cov_info object created above will be copied and added as "data" to all whole genome blocks
    // This process is happening in place.
    ptBlock_add_data_to_all_blocks_stHash(whole_genome_block_table,
                                          cov_info,
                                          destruct_cov_info_data,
                                          copy_cov_info_data,
                                          extend_cov_info_data);

    // print len/number stats for the whole genome block table
    fprintf(stderr, "[%s] Created block table for whole genome  : tot_len=%ld, number=%ld\n", get_timestamp(),
            ptBlock_get_total_length_by_rf(whole_genome_block_table),
            ptBlock_get_total_number(whole_genome_block_table));


    // add annotation and coverage blocks to the whole genome blocks in place
    ptBlock_extend_block_tables(whole_genome_block_table, coverage_block_table);
    fprintf(stderr, "[%s] Added 0-coverage whole genome blocks to coverage block tables : tot_len=%ld, number=%ld\n",
            get_timestamp(),
            ptBlock_get_total_length_by_rf(whole_genome_block_table),
            ptBlock_get_total_number(whole_genome_block_table));
    for (int i = 0; i < stList_length(annotation_block_table_list); i++) {
        ptBlock_extend_block_tables(whole_genome_block_table, stList_get(annotation_block_table_list, i));
    }
    fprintf(stderr, "[%s] Added annotation blocks to coverage block tables: tot_len=%ld, number=%ld\n", get_timestamp(),
            ptBlock_get_total_length_by_rf(whole_genome_block_table),
            ptBlock_get_total_number(whole_genome_block_table));

    fprintf(stderr, "[%s] Started sorting and merging blocks\n", get_timestamp());
    //sort
    ptBlock_sort_stHash_by_rfs(whole_genome_block_table);
    fprintf(stderr, "[%s] Merged blocks : tot_len=%ld, number=%ld\n", get_timestamp(),
            ptBlock_get_total_length_by_rf(whole_genome_block_table),
            ptBlock_get_total_number(whole_genome_block_table));

    //merge and create the final block table
    stHash *final_block_table = ptBlock_merge_blocks_per_contig_by_rf_v2(whole_genome_block_table);

    fprintf(stderr, "[%s] Created final block table : tot_len=%ld, number=%ld\n", get_timestamp(),
            ptBlock_get_total_length_by_rf(final_block_table),
            ptBlock_get_total_number(final_block_table));


    // free unmerged blocks
    stHash_destruct(coverage_block_table);
    stHash_destruct(whole_genome_block_table);
    stList_destruct(annotation_block_table_list);

    return final_block_table;
}
