//
// Created by mobin on 7/16/23.
//

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
#include <string.h>
#include "ptBlock.h"
#include "ptVariant.h"
#include "ptAlignment.h"
#include "ptMarker.h"
#include "tpool.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

bool test_sortingBlocks(){
    stHash* blocks_per_contig = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL,
                                                    (void (*)(void *)) stList_destruct);

    char ctg1_name[10] = "ctg1";
    char ctg2_name[10] = "ctg2";

    // create list of unsorted blocks

    stList* ctg1_blocks = stList_construct3(0, ptBlock_destruct);
    stList* ctg2_blocks = stList_construct3(0, ptBlock_destruct);

    stList_append(ctg1_blocks, ptBlock_construct(15, 50, -1, -1, -1, -1));
    stList_append(ctg1_blocks, ptBlock_construct(5, 15, -1, -1, -1, -1));
    stList_append(ctg1_blocks, ptBlock_construct(10, 20, -1, -1, -1, -1));
    stList_append(ctg1_blocks, ptBlock_construct(0, 10, -1, -1, -1, -1));
    stList_append(ctg1_blocks, ptBlock_construct(15, 60, -1, -1, -1, -1));

    stList_append(ctg2_blocks, ptBlock_construct(50, 60, -1, -1, -1, -1));
    stList_append(ctg2_blocks, ptBlock_construct(5, 15, -1, -1, -1, -1));
    stList_append(ctg2_blocks, ptBlock_construct(10, 10, -1, -1, -1, -1));
    stList_append(ctg2_blocks, ptBlock_construct(8, 8, -1, -1, -1, -1));
    stList_append(ctg2_blocks, ptBlock_construct(0, 10, -1, -1, -1, -1));


    // add blocks to the table
    stHash_insert(blocks_per_contig, copyString(ctg1_name), ctg1_blocks);
    stHash_insert(blocks_per_contig, copyString(ctg2_name), ctg2_blocks);

    // sort blocks by rfs
    ptBlock_sort_stHash_by_rfs(blocks_per_contig); // sort in place

    // get sorted blocks
    stList* ctg1_sorted_blocks = stHash_search(blocks_per_contig, ctg1_name);
    stList* ctg2_sorted_blocks = stHash_search(blocks_per_contig, ctg2_name);

    int ctg1_sorted_start[5] = {0, 5, 10, 15 ,15};
    int ctg2_sorted_start[5] = {0, 5, 8, 10 ,50};

    bool allPassed = true;
    for(int i =0; i < 5; i ++){
        ptBlock* ctg1_block = stList_get(ctg1_sorted_blocks, i);
        ptBlock* ctg2_block = stList_get(ctg2_sorted_blocks, i);

        allPassed &= ctg1_sorted_start[i] == ctg1_block->rfs;
        allPassed &= ctg2_sorted_start[i] == ctg2_block->rfs;
    }

    stHash_destruct(blocks_per_contig);

    return allPassed;
}

bool test_mergingBlocksWithoutCount(){
    stHash* blocks_per_contig = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL,
                                                  (void (*)(void *)) stList_destruct);

    char ctg1_name[10] = "ctg1";
    char ctg2_name[10] = "ctg2";

    // create list of unsorted blocks

    stList* ctg1_blocks = stList_construct3(0, ptBlock_destruct);
    stList* ctg2_blocks = stList_construct3(0, ptBlock_destruct);

    stList_append(ctg1_blocks, ptBlock_construct(30, 50, -1, -1, -1, -1));
    stList_append(ctg1_blocks, ptBlock_construct(5, 6, -1, -1, -1, -1));
    stList_append(ctg1_blocks, ptBlock_construct(10, 20, -1, -1, -1, -1));
    stList_append(ctg1_blocks, ptBlock_construct(0, 10, -1, -1, -1, -1));
    stList_append(ctg1_blocks, ptBlock_construct(50, 60, -1, -1, -1, -1));

    stList_append(ctg2_blocks, ptBlock_construct(50, 60, -1, -1, -1, -1));
    stList_append(ctg2_blocks, ptBlock_construct(5, 15, -1, -1, -1, -1));
    stList_append(ctg2_blocks, ptBlock_construct(10, 10, -1, -1, -1, -1));
    stList_append(ctg2_blocks, ptBlock_construct(8, 8, -1, -1, -1, -1));
    stList_append(ctg2_blocks, ptBlock_construct(0, 10, -1, -1, -1, -1));


    // add blocks to the table
    stHash_insert(blocks_per_contig, copyString(ctg1_name), ctg1_blocks);
    stHash_insert(blocks_per_contig, copyString(ctg2_name), ctg2_blocks);

    // sort blocks by rfs
    ptBlock_sort_stHash_by_rfs(blocks_per_contig); // sort in place

    // merge blocks by rfs
    stHash* merged_blocks_per_contig = ptBlock_merge_blocks_per_contig_by_rf(blocks_per_contig);

    // get merged blocks
    stList* ctg1_merged_blocks = stHash_search(merged_blocks_per_contig, ctg1_name);
    stList* ctg2_merged_blocks = stHash_search(merged_blocks_per_contig, ctg2_name);


    // truth
    int ctg1_merged_blocks_start[2] = {0, 30};
    int ctg1_merged_blocks_end[2] = {20, 60};

    int ctg2_merged_blocks_start[2] = {0, 50};
    int ctg2_merged_blocks_end[2] = {15, 60};


    // check number of merged blocks
    if (stList_length(ctg1_merged_blocks) != 2 || stList_length(ctg2_merged_blocks) != 2){
        return false;
    }

    bool allPassed = true;

    for(int i =0; i < 2; i ++){
        ptBlock* ctg1_merged_block = stList_get(ctg1_merged_blocks, i);
        ptBlock* ctg2_merged_block = stList_get(ctg2_merged_blocks, i);

        allPassed &= ctg1_merged_blocks_start[i] == ctg1_merged_block->rfs;
        allPassed &= ctg1_merged_blocks_end[i] == ctg1_merged_block->rfe;

        allPassed &= ctg2_merged_blocks_start[i] == ctg2_merged_block->rfs;
        allPassed &= ctg2_merged_blocks_end[i] == ctg2_merged_block->rfe;
    }

    stHash_destruct(blocks_per_contig);
    stHash_destruct(merged_blocks_per_contig);

    return allPassed;
}


bool test_mergingBlocksWithCount_v2(){
    stHash* blocks_per_contig = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL,
                                                  (void (*)(void *)) stList_destruct);

    char ctg1_name[10] = "ctg1";
    char ctg2_name[10] = "ctg2";

    // create list of unsorted blocks

    stList* ctg1_blocks = stList_construct3(0, ptBlock_destruct);
    stList* ctg2_blocks = stList_construct3(0, ptBlock_destruct);

    stList_append(ctg1_blocks, ptBlock_construct_with_count(30, 50, -1, -1, -1, -1, 1));
    stList_append(ctg1_blocks, ptBlock_construct_with_count(5, 6, -1, -1, -1, -1, 1));
    stList_append(ctg1_blocks, ptBlock_construct_with_count(10, 20, -1, -1, -1, -1, 1));
    stList_append(ctg1_blocks, ptBlock_construct_with_count(0, 10, -1, -1, -1, -1, 1));
    stList_append(ctg1_blocks, ptBlock_construct_with_count(50, 60, -1, -1, -1, -1, 1));

    stList_append(ctg2_blocks, ptBlock_construct_with_count(50, 60, -1, -1, -1, -1, 1));
    stList_append(ctg2_blocks, ptBlock_construct_with_count(5, 15, -1, -1, -1, -1, 1));
    stList_append(ctg2_blocks, ptBlock_construct_with_count(10, 10, -1, -1, -1, -1, 1));
    stList_append(ctg2_blocks, ptBlock_construct_with_count(8, 8, -1, -1, -1, -1, 1));
    stList_append(ctg2_blocks, ptBlock_construct_with_count(0, 10, -1, -1, -1, -1, 1));


    // add blocks to the table
    stHash_insert(blocks_per_contig, copyString(ctg1_name), ctg1_blocks);
    stHash_insert(blocks_per_contig, copyString(ctg2_name), ctg2_blocks);

    // sort blocks by rfs
    ptBlock_sort_stHash_by_rfs(blocks_per_contig); // sort in place

    // merge blocks by rfs
    stHash* merged_blocks_per_contig = ptBlock_merge_blocks_per_contig_by_rf_v2(blocks_per_contig);

    // get merged blocks
    stList* ctg1_merged_blocks = stHash_search(merged_blocks_per_contig, ctg1_name);
    stList* ctg2_merged_blocks = stHash_search(merged_blocks_per_contig, ctg2_name);


    // truth
    int ctg1_merged_blocks_start[8] = {0,5,7,10,11,30,50,51};
    int ctg1_merged_blocks_end[8] = {4,6,9,10,20,49,50,60};
    int ctg1_merged_blocks_count[8] = {1,2,1,2,1,1,2,1};

    int ctg2_merged_blocks_start[7] = {0,5,8,9,10,11, 50};
    int ctg2_merged_blocks_end[7] = {4,7,8,9,10,15, 60};
    int ctg2_merged_blocks_count[7] = {1, 2, 3, 2, 3, 1, 1};

    // check number of merged blocks
    if (stList_length(ctg1_merged_blocks) != 8 || stList_length(ctg2_merged_blocks) != 7){
        return false;
    }

    bool allPassed = true;
    for(int i =0; i < 8; i ++){
        ptBlock* ctg1_merged_block = stList_get(ctg1_merged_blocks, i);

        allPassed &= ctg1_merged_blocks_start[i] == ctg1_merged_block->rfs;
        allPassed &= ctg1_merged_blocks_end[i] == ctg1_merged_block->rfe;
        allPassed &= ctg1_merged_blocks_count[i] == ptBlock_get_count(ctg1_merged_block);
    }

    for(int i =0; i < 7; i ++){
        ptBlock* ctg2_merged_block = stList_get(ctg2_merged_blocks, i);

        allPassed &= ctg2_merged_blocks_start[i] == ctg2_merged_block->rfs;
        allPassed &= ctg2_merged_blocks_end[i] == ctg2_merged_block->rfe;
        allPassed &= ctg2_merged_blocks_count[i] == ptBlock_get_count(ctg2_merged_block);
    }

    stHash_destruct(blocks_per_contig);
    stHash_destruct(merged_blocks_per_contig);

    return allPassed;
}




int main(int argc, char *argv[]) {
    fprintf(stdout, "Start testing ....\n");
    fprintf(stdout, "Test sorting blocks:");
    fprintf(stdout, test_sortingBlocks() ? "\x1B[32m PASSED \x1B[0m\n" : "\x1B[31m FAILED \x1B[0m\n");
    fprintf(stdout, "Test merging blocks without count:");
    fprintf(stdout, test_mergingBlocksWithoutCount() ? "\x1B[32m PASSED \x1B[0m\n" : "\x1B[31m FAILED \x1B[0m\n");
    fprintf(stdout, "Test merging blocks with count (v2):");
    fprintf(stdout, test_mergingBlocksWithCount_v2() ? "\x1B[32m PASSED \x1B[0m\n" : "\x1B[31m FAILED \x1B[0m\n");
}
