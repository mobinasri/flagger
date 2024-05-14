#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include "sonLib.h"
#include "ptBlock.h"
#include "ptAlignment.h"
#include "common.h"
#include "hmm.h"

bool Test_parse_bed(char *bed_path) {
    stHash *blocks_per_contig_truth = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL,
                                                        (void (*)(void *)) stList_destruct);

    //ctg1
    stList *ctg1_list = stList_construct3(0, ptBlock_destruct);
    stList_append(ctg1_list, ptBlock_construct(10, 49,
                                               -1, -1,
                                               -1, -1));
    stList_append(ctg1_list, ptBlock_construct(100, 149,
                                               -1, -1,
                                               -1, -1));
    stHash_insert(blocks_per_contig_truth, "ctg1", ctg1_list);

    //ctg2
    stList *ctg2_list = stList_construct3(0, ptBlock_destruct);
    stList_append(ctg2_list, ptBlock_construct(0, 9,
                                               -1, -1,
                                               -1, -1));
    stList_append(ctg2_list, ptBlock_construct(10, 14,
                                               -1, -1,
                                               -1, -1));
    stHash_insert(blocks_per_contig_truth, "ctg2", ctg2_list);

    //ctg3
    stList *ctg3_list = stList_construct3(0, ptBlock_destruct);
    stList_append(ctg3_list, ptBlock_construct(0, 99,
                                               -1, -1,
                                               -1, -1));
    stHash_insert(blocks_per_contig_truth, "ctg3", ctg3_list);

    // parse blocks from bed file
    stHash *blocks_per_contig_out = ptBlock_parse_bed(bed_path);
    bool is_equal = ptBlock_is_equal_stHash(blocks_per_contig_truth, blocks_per_contig_out);

    stHash_destruct(blocks_per_contig_truth);
    stHash_destruct(blocks_per_contig_out);

    return is_equal;
}

bool Test_split_blocks(char *bed_path) {
    int split_number = 4;
    stList *ctg_list;
    stHash *blocks_per_contig = ptBlock_parse_bed(bed_path);
    stList *batches_out = ptBlock_split_into_batches(blocks_per_contig, split_number);

    // make truth batches
    stList *batches_truth = stList_construct3(0, stHash_destruct);

    //batch_0
    stHash *batch_0 = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL,
                                        (void (*)(void *)) stList_destruct);
    ctg_list = stList_construct3(0, ptBlock_destruct);
    stList_append(ctg_list, ptBlock_construct(10, 49,
                                              -1, -1,
                                              -1, -1));
    stList_append(ctg_list, ptBlock_construct(100, 111,
                                              -1, -1,
                                              -1, -1));
    stHash_insert(batch_0, "ctg1", ctg_list);

    //batch_1
    stHash *batch_1 = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL,
                                        (void (*)(void *)) stList_destruct);
    ctg_list = stList_construct3(0, ptBlock_destruct);
    stList_append(ctg_list, ptBlock_construct(112, 149,
                                              -1, -1,
                                              -1, -1));
    stHash_insert(batch_1, "ctg1", ctg_list);

    ctg_list = stList_construct3(0, ptBlock_destruct);
    stList_append(ctg_list, ptBlock_construct(0, 9,
                                              -1, -1,
                                              -1, -1));
    stList_append(ctg_list, ptBlock_construct(10, 13,
                                              -1, -1,
                                              -1, -1));
    stHash_insert(batch_1, "ctg2", ctg_list);

    //batch_2
    stHash *batch_2 = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL,
                                        (void (*)(void *)) stList_destruct);
    ctg_list = stList_construct3(0, ptBlock_destruct);
    stList_append(ctg_list, ptBlock_construct(14, 14,
                                              -1, -1,
                                              -1, -1));
    stHash_insert(batch_2, "ctg2", ctg_list);

    ctg_list = stList_construct3(0, ptBlock_destruct);
    stList_append(ctg_list, ptBlock_construct(0, 50,
                                              -1, -1,
                                              -1, -1));
    stHash_insert(batch_2, "ctg3", ctg_list);

    //batch_3
    stHash *batch_3 = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL,
                                        (void (*)(void *)) stList_destruct);
    ctg_list = stList_construct3(0, ptBlock_destruct);
    stList_append(ctg_list, ptBlock_construct(51, 99,
                                              -1, -1,
                                              -1, -1));
    stHash_insert(batch_3, "ctg3", ctg_list);


    //add all batches

    stList_append(batches_truth, batch_0);
    stList_append(batches_truth, batch_1);
    stList_append(batches_truth, batch_2);
    stList_append(batches_truth, batch_3);

    if (stList_length(batches_out) != stList_length(batches_truth)) {
        return false;
    }

    bool test_passed = true;
    for (int i = 0; i < stList_length(batches_truth); i++) {
        test_passed &= ptBlock_is_equal_stHash(stList_get(batches_out, i), stList_get(batches_truth, i));
    }
    stList_destruct(batches_truth);
    stList_destruct(batches_out);

    return test_passed;
}


bool test_CoverageInfo_getAnnotationFlagFromArray() {
    int x[3];
    x[0] = 1;
    x[1] = 3;
    x[2] = 9;
    uint64_t flag = CoverageInfo_getAnnotationFlagFromArray(x, 3);
    return flag == 0x105;
}


bool test_CoverageInfo_setAndGetRegionIndex() {
    CoverageInfo *covInfo = CoverageInfo_construct(0, 0, 0, 0);

    bool correct = true;

    CoverageInfo_setRegionIndex(covInfo, 3);
    correct &= (CoverageInfo_getRegionIndex(covInfo) == 3);

    CoverageInfo_setRegionIndex(covInfo, 6);
    correct &= (CoverageInfo_getRegionIndex(covInfo) == 6);

    CoverageInfo_destruct(covInfo);
    return correct;
}


bool test_CoverageInfo_getAnnotationIndices() {
    CoverageInfo *covInfo = CoverageInfo_construct(0, 0, 0, 0);
    bool correct = true;
    int x[3];
    x[0] = 1;
    x[1] = 3;
    x[2] = 9;
    covInfo->annotation_flag = CoverageInfo_getAnnotationFlagFromArray(x, 3);

    int len = -1;
    int *y = CoverageInfo_getAnnotationIndices(covInfo, &len);
    if (len != 3) return false;
    correct &= (y[0] == 1) && (y[1] == 3) && (y[2] == 9);
    free(y);

    // another test for 0 annotation
    x[0] = 0;
    covInfo->annotation_flag = CoverageInfo_getAnnotationFlagFromArray(x, 1);

    y = CoverageInfo_getAnnotationIndices(covInfo, &len);
    if (len != 1) return false;
    correct &= (y[0] == 0);

    free(y);
    CoverageInfo_destruct(covInfo);
    return correct;
}


int main(int argc, char *argv[]) {
    char bed_path[200] = "tests/test_files/ptBlock/test.bed";

    bool all_tests_passed = true;
    // test 1
    bool test_parse_bed_passed = Test_parse_bed(bed_path);
    all_tests_passed &= test_parse_bed_passed;
    printf("Test parse bed:");
    printf(test_parse_bed_passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");

    // test 2
    bool test_split_blocks_passed = Test_split_blocks(bed_path);
    all_tests_passed &= test_split_blocks_passed;
    printf("Test split blocks:");
    printf(test_split_blocks_passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");

    // test 3
    bool test_CoverageInfo_setAndGetRegionIndex_passed = test_CoverageInfo_setAndGetRegionIndex();
    all_tests_passed &= test_CoverageInfo_setAndGetRegionIndex_passed;
    printf("Test CoverageInfo_setAndGetRegionIndex:");
    printf(test_CoverageInfo_setAndGetRegionIndex_passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");


    // test 4
    bool test_CoverageInfo_getAnnotationFlagFromArray_passed = test_CoverageInfo_getAnnotationFlagFromArray();
    all_tests_passed &= test_CoverageInfo_getAnnotationFlagFromArray_passed;
    printf("Test CoverageInfo_getAnnotationFlagFromArray:");
    printf(test_CoverageInfo_getAnnotationFlagFromArray_passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");


    // test 5
    bool test_CoverageInfo_getAnnotationIndices_passed = test_CoverageInfo_getAnnotationIndices();
    all_tests_passed &= test_CoverageInfo_getAnnotationIndices_passed;
    printf("Test CoverageInfo_getAnnotationIndices:");
    printf(test_CoverageInfo_getAnnotationIndices_passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");

    if (all_tests_passed)
        return 0;
    else
        return 1;
}
