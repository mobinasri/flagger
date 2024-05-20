#include "common.h"
#include "summary_table.h"
#include "ptBlock.h"
#include "chunk.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "sonLib.h"


bool test_SummaryTable_increment() {
    int numberOfRows = 2;
    int numberOfColumns = 2;
    SummaryTable *summaryTable = SummaryTable_construct(numberOfRows, numberOfColumns);

    SummaryTable_increment(summaryTable, 0, 0, 10);
    SummaryTable_increment(summaryTable, 0, 0, 10);
    SummaryTable_increment(summaryTable, 0, 0, 10);
    SummaryTable_increment(summaryTable, 0, 1, 10);
    SummaryTable_increment(summaryTable, 0, 1, 10);

    bool correct = true;
    double **table = summaryTable->table;

    //counts
    correct &= summaryTable->table[0][0] == 30;
    correct &= summaryTable->table[0][1] == 20;
    correct &= summaryTable->table[1][0] == 0;
    correct &= summaryTable->table[1][1] == 0;

    // percentages
    correct &= summaryTable->tablePercentage[0][0] == 60;
    correct &= summaryTable->tablePercentage[0][1] == 40;
    correct &= summaryTable->tablePercentage[1][0] == 0;
    correct &= summaryTable->tablePercentage[1][1] == 0;

    // get total per row
    correct &= summaryTable->totalPerRow[0] == 50;
    correct &= summaryTable->totalPerRow[1] == 0;


    SummaryTable_destruct(summaryTable);
    return correct;
}

bool test_SummaryTable_getRowString() {
    int numberOfRows = 2;
    int numberOfColumns = 2;
    SummaryTable *summaryTable = SummaryTable_construct(numberOfRows, numberOfColumns);

    SummaryTable_increment(summaryTable, 0, 0, 10);
    SummaryTable_increment(summaryTable, 0, 0, 10);
    SummaryTable_increment(summaryTable, 0, 0, 10);
    SummaryTable_increment(summaryTable, 0, 1, 10);
    SummaryTable_increment(summaryTable, 0, 1, 10);

    bool correct = true;


    correct &= strcmp(SummaryTable_getRowString(summaryTable, 0, ','), "30.00,20.00") == 0;
    correct &= strcmp(SummaryTable_getRowString(summaryTable, 1, ','), "0.00,0.00") == 0;
    correct &= strcmp(SummaryTable_getRowStringPercentage(summaryTable, 0, ','), "60.00,40.00") == 0;
    correct &= strcmp(SummaryTable_getRowStringPercentage(summaryTable, 1, ','), "0.00,0.00") == 0;

    SummaryTable_destruct(summaryTable);
    return correct;
}

bool test_SummaryTableList_increment() {
    stList *categoryNames1 = stList_construct3(0, free);
    stList_append(categoryNames1, copyString("cat1_0"));
    stList_append(categoryNames1, copyString("cat1_1"));
    stList *categoryNames2 = stList_construct3(0, free);
    stList_append(categoryNames2, copyString("cat2_0"));
    stList_append(categoryNames2, copyString("cat2_1"));

    int numberOfRows = 2;
    int numberOfColumns = 2;
    SummaryTableList *summaryTableList = SummaryTableList_construct(categoryNames1,
                                                                    categoryNames2,
                                                                    numberOfRows,
                                                                    numberOfColumns);
    int cat1Index = 0;
    int cat2Index = 0;
    SummaryTableList_increment(summaryTableList, cat1Index, cat2Index, 0, 0, 10);
    SummaryTableList_increment(summaryTableList, cat1Index, cat2Index, 0, 0, 10);
    SummaryTableList_increment(summaryTableList, cat1Index, cat2Index, 0, 0, 10);
    SummaryTableList_increment(summaryTableList, cat1Index, cat2Index, 0, 1, 10);
    SummaryTableList_increment(summaryTableList, cat1Index, cat2Index, 0, 1, 10);

    cat1Index = 1;
    cat2Index = 1;
    SummaryTableList_increment(summaryTableList, cat1Index, cat2Index, 0, 0, 40);
    SummaryTableList_increment(summaryTableList, cat1Index, cat2Index, 0, 0, 20);
    SummaryTableList_increment(summaryTableList, cat1Index, cat2Index, 0, 1, 10);
    SummaryTableList_increment(summaryTableList, cat1Index, cat2Index, 0, 1, 10);

    bool correct = true;
    cat1Index = 0;
    cat2Index = 0;
    correct &= SummaryTableList_getValue(summaryTableList, cat1Index, cat2Index, 0, 0) == 30.0;
    correct &= SummaryTableList_getValue(summaryTableList, cat1Index, cat2Index, 0, 1) == 20.0;
    correct &= SummaryTableList_getValue(summaryTableList, cat1Index, cat2Index, 1, 0) == 0.0;
    correct &= SummaryTableList_getValue(summaryTableList, cat1Index, cat2Index, 1, 1) == 0.0;
    correct &= SummaryTableList_getValuePercentage(summaryTableList, cat1Index, cat2Index, 0, 0) == 60.0;
    correct &= SummaryTableList_getValuePercentage(summaryTableList, cat1Index, cat2Index, 0, 1) == 40.0;
    correct &= SummaryTableList_getValuePercentage(summaryTableList, cat1Index, cat2Index, 1, 0) == 0.0;
    correct &= SummaryTableList_getValuePercentage(summaryTableList, cat1Index, cat2Index, 1, 1) == 0.0;

    cat1Index = 1;
    cat2Index = 0;
    correct &= SummaryTableList_getValue(summaryTableList, cat1Index, cat2Index, 0, 0) == 0.0;
    correct &= SummaryTableList_getValue(summaryTableList, cat1Index, cat2Index, 0, 1) == 0.0;
    correct &= SummaryTableList_getValue(summaryTableList, cat1Index, cat2Index, 1, 0) == 0.0;
    correct &= SummaryTableList_getValue(summaryTableList, cat1Index, cat2Index, 1, 1) == 0.0;
    correct &= SummaryTableList_getValuePercentage(summaryTableList, cat1Index, cat2Index, 0, 0) == 0.0;
    correct &= SummaryTableList_getValuePercentage(summaryTableList, cat1Index, cat2Index, 0, 1) == 0.0;
    correct &= SummaryTableList_getValuePercentage(summaryTableList, cat1Index, cat2Index, 1, 0) == 0.0;
    correct &= SummaryTableList_getValuePercentage(summaryTableList, cat1Index, cat2Index, 1, 1) == 0.0;

    cat1Index = 0;
    cat2Index = 1;
    correct &= SummaryTableList_getValue(summaryTableList, cat1Index, cat2Index, 0, 0) == 0.0;
    correct &= SummaryTableList_getValue(summaryTableList, cat1Index, cat2Index, 0, 1) == 0.0;
    correct &= SummaryTableList_getValue(summaryTableList, cat1Index, cat2Index, 1, 0) == 0.0;
    correct &= SummaryTableList_getValue(summaryTableList, cat1Index, cat2Index, 1, 1) == 0.0;
    correct &= SummaryTableList_getValuePercentage(summaryTableList, cat1Index, cat2Index, 0, 0) == 0.0;
    correct &= SummaryTableList_getValuePercentage(summaryTableList, cat1Index, cat2Index, 0, 1) == 0.0;
    correct &= SummaryTableList_getValuePercentage(summaryTableList, cat1Index, cat2Index, 1, 0) == 0.0;
    correct &= SummaryTableList_getValuePercentage(summaryTableList, cat1Index, cat2Index, 1, 1) == 0.0;

    cat1Index = 1;
    cat2Index = 1;
    correct &= SummaryTableList_getValue(summaryTableList, cat1Index, cat2Index, 0, 0) == 60.0;
    correct &= SummaryTableList_getValue(summaryTableList, cat1Index, cat2Index, 0, 1) == 20.0;
    correct &= SummaryTableList_getValue(summaryTableList, cat1Index, cat2Index, 1, 0) == 0.0;
    correct &= SummaryTableList_getValue(summaryTableList, cat1Index, cat2Index, 1, 1) == 0.0;
    correct &= SummaryTableList_getValuePercentage(summaryTableList, cat1Index, cat2Index, 0, 0) == 75.0;
    correct &= SummaryTableList_getValuePercentage(summaryTableList, cat1Index, cat2Index, 0, 1) == 25.0;
    correct &= SummaryTableList_getValuePercentage(summaryTableList, cat1Index, cat2Index, 1, 0) == 0.0;
    correct &= SummaryTableList_getValuePercentage(summaryTableList, cat1Index, cat2Index, 1, 1) == 0.0;

    SummaryTableList_destruct(summaryTableList);
    stList_destruct(categoryNames1);
    stList_destruct(categoryNames2);
    return true;
}


bool test_SummaryTableList_getRowString() {
    stList *categoryNames1 = stList_construct3(0, free);
    stList_append(categoryNames1, copyString("cat1_0"));
    stList_append(categoryNames1, copyString("cat1_1"));
    stList *categoryNames2 = stList_construct3(0, free);
    stList_append(categoryNames2, copyString("cat2_0"));
    stList_append(categoryNames2, copyString("cat2_1"));

    int numberOfRows = 2;
    int numberOfColumns = 2;
    SummaryTableList *summaryTableList = SummaryTableList_construct(categoryNames1,
                                                                    categoryNames2,
                                                                    numberOfRows,
                                                                    numberOfColumns);
    int cat1Index = 0;
    int cat2Index = 0;
    SummaryTableList_increment(summaryTableList, cat1Index, cat2Index, 0, 0, 10);
    SummaryTableList_increment(summaryTableList, cat1Index, cat2Index, 0, 0, 10);
    SummaryTableList_increment(summaryTableList, cat1Index, cat2Index, 0, 0, 10);
    SummaryTableList_increment(summaryTableList, cat1Index, cat2Index, 0, 1, 10);
    SummaryTableList_increment(summaryTableList, cat1Index, cat2Index, 0, 1, 10);

    cat1Index = 1;
    cat2Index = 1;
    SummaryTableList_increment(summaryTableList, cat1Index, cat2Index, 0, 0, 40);
    SummaryTableList_increment(summaryTableList, cat1Index, cat2Index, 0, 0, 20);
    SummaryTableList_increment(summaryTableList, cat1Index, cat2Index, 0, 1, 10);
    SummaryTableList_increment(summaryTableList, cat1Index, cat2Index, 0, 1, 10);

    bool correct = true;
    cat1Index = 0;
    cat2Index = 0;
    correct &= strcmp(
            SummaryTableList_getRowString(summaryTableList, cat1Index, cat2Index, 0, ','),
            "30.00,20.00") == 0;
    correct &= strcmp(
            SummaryTableList_getRowString(summaryTableList, cat1Index, cat2Index, 1, ','),
            "0.00,0.00") == 0;
    correct &= strcmp(
            SummaryTableList_getRowStringPercentage(summaryTableList, cat1Index, cat2Index, 0, ','),
            "60.00,40.00") == 0;
    correct &= strcmp(
            SummaryTableList_getRowStringPercentage(summaryTableList, cat1Index, cat2Index, 1, ','),
            "0.00,0.00") == 0;


    cat1Index = 1;
    cat2Index = 0;
    correct &= strcmp(
            SummaryTableList_getRowString(summaryTableList, cat1Index, cat2Index, 0, ','),
            "0.00,0.00") == 0;
    correct &= strcmp(
            SummaryTableList_getRowString(summaryTableList, cat1Index, cat2Index, 1, ','),
            "0.00,0.00") == 0;
    correct &= strcmp(
            SummaryTableList_getRowStringPercentage(summaryTableList, cat1Index, cat2Index, 0, ','),
            "0.00,0.00") == 0;
    correct &= strcmp(
            SummaryTableList_getRowStringPercentage(summaryTableList, cat1Index, cat2Index, 1, ','),
            "0.00,0.00") == 0;

    cat1Index = 0;
    cat2Index = 1;
    correct &= strcmp(
            SummaryTableList_getRowString(summaryTableList, cat1Index, cat2Index, 0, ','),
            "0.00,0.00") == 0;
    correct &= strcmp(
            SummaryTableList_getRowString(summaryTableList, cat1Index, cat2Index, 1, ','),
            "0.00,0.00") == 0;
    correct &= strcmp(
            SummaryTableList_getRowStringPercentage(summaryTableList, cat1Index, cat2Index, 0, ','),
            "0.00,0.00") == 0;
    correct &= strcmp(
            SummaryTableList_getRowStringPercentage(summaryTableList, cat1Index, cat2Index, 1, ','),
            "0.00,0.00") == 0;

    cat1Index = 1;
    cat2Index = 1;
    correct &= strcmp(
            SummaryTableList_getRowString(summaryTableList, cat1Index, cat2Index, 0, ','),
            "60.00,20.00") == 0;
    correct &= strcmp(
            SummaryTableList_getRowString(summaryTableList, cat1Index, cat2Index, 1, ','),
            "0.00,0.00") == 0;
    correct &= strcmp(
            SummaryTableList_getRowStringPercentage(summaryTableList, cat1Index, cat2Index, 0, ','),
            "75.00,25.00") == 0;
    correct &= strcmp(
            SummaryTableList_getRowStringPercentage(summaryTableList, cat1Index, cat2Index, 1, ','),
            "0.00,0.00") == 0;

    SummaryTableList_destruct(summaryTableList);
    stList_destruct(categoryNames1);
    stList_destruct(categoryNames2);
    return true;
}

int test_ptBlock_updateSummaryTableListWithIterator(const char *covPath, const char *binArrayFilePath,
                                                    bool useChunkIterator, bool isMetricOverlapBased) {
    // whole genome
    double **wholeGenomeBin1 = Double_construct2DArray(4, 5);
    if (isMetricOverlapBased) {
        wholeGenomeBin1[1][3] = 1.0;
    } else {
        wholeGenomeBin1[1][3] = 2.0;
    }

    double **wholeGenomeBin2 = Double_construct2DArray(4, 5);
    if (isMetricOverlapBased) {
        wholeGenomeBin2[0][0] = 1.0;
        wholeGenomeBin2[2][2] = 1.0;
        wholeGenomeBin2[2][4] = 1.0; // undefined label
        wholeGenomeBin2[3][1] = 1.0;
        wholeGenomeBin2[3][3] = 1.0;
    } else {
        wholeGenomeBin2[0][0] = 2.0;
        wholeGenomeBin2[0][4] = 1.0; // undefined label
        wholeGenomeBin2[2][1] = 1.0;
        wholeGenomeBin2[2][2] = 6.0;
        wholeGenomeBin2[2][4] = 6.0; // undefined label
        wholeGenomeBin2[3][0] = 1.0;
        wholeGenomeBin2[3][1] = 4.0;
        wholeGenomeBin2[3][3] = 3.0;
    }

    // annotation 1
    double **annotation1Bin1 = Double_construct2DArray(4, 5);
    if (isMetricOverlapBased) {
        annotation1Bin1[0][0] = 1.0;
    } else {
        annotation1Bin1[0][0] = 2.0;
    }

    double **annotation1Bin2 = Double_construct2DArray(4, 5);
    if (isMetricOverlapBased) {
        annotation1Bin2[2][4] = 1.0; // undefined label
    } else {
        annotation1Bin2[2][1] = 1.0;
        annotation1Bin2[2][4] = 2.0; // undefined label
    }

    // annotation 2
    double **annotation2Bin1 = Double_construct2DArray(4, 5);
    if (isMetricOverlapBased) {
        annotation2Bin1[1][3] = 1.0;
    } else {
        annotation2Bin1[1][3] = 2.0;
    }

    double **annotation2Bin2 = Double_construct2DArray(4, 5);
    if (isMetricOverlapBased) {
        annotation2Bin2[2][2] = 1.0;
        annotation2Bin2[2][4] = 1.0; // undefined label
        annotation2Bin2[3][1] = 1.0;
        annotation2Bin2[3][3] = 1.0;
    } else {
        annotation2Bin2[2][2] = 4.0;
        annotation2Bin2[2][4] = 4.0; // undefined label
        annotation2Bin2[3][1] = 4.0;
        annotation2Bin2[3][3] = 3.0;
    }

    void *iterator;
    ptBlock *(*getNextBlock)(void *, char *);
    void *(*copyIterator)(void *);
    void (*destructIterator)(void *);
    void (*resetIterator)(void *);
    CoverageHeader *header = NULL;
    ChunksCreator *chunksCreator = NULL;
    stHash *blocksPerContig = NULL;

    if (useChunkIterator) {
        int windowLen = 1;
        int chunkCanonicalLen = 10;
        int nThreads = 2;
        chunksCreator = ChunksCreator_constructFromCov(covPath,
                                                       NULL,
                                                       chunkCanonicalLen,
                                                       nThreads,
                                                       windowLen);
        if (ChunksCreator_parseChunks(chunksCreator) != 0) {
            return false;
        }
        iterator = (void *) ChunkIterator_construct(chunksCreator);
        copyIterator = ChunkIterator_copy;
        destructIterator = ChunkIterator_destruct;
        resetIterator = ChunkIterator_reset;
        getNextBlock = ChunkIterator_getNextPtBlock;
        header = chunksCreator->header;
    } else {
        stHash *blocksPerContig = ptBlock_parse_coverage_info_blocks(covPath);
        iterator = (void *) ptBlockItrPerContig_construct(blocksPerContig);
        copyIterator = ptBlockItrPerContig_copy;
        destructIterator = ptBlockItrPerContig_destruct;
        resetIterator = ptBlockItrPerContig_reset;
        getNextBlock = ptBlockItrPerContig_next;
        header = CoverageHeader_construct(covPath);
    }


    IntBinArray *binArray = IntBinArray_constructFromFile(binArrayFilePath);

    int numberOfLabels = header->numberOfLabels;
    stList *categoryNames1 = header->annotationNames;
    stList *categoryNames2 = binArray->names;

    int numberOfRows = numberOfLabels;
    int numberOfColumns = numberOfLabels;
    SummaryTableList *summaryTableList = SummaryTableList_construct(categoryNames1,
                                                                    categoryNames2,
                                                                    numberOfRows,
                                                                    numberOfColumns +
                                                                    1); // +1 for undefined state (label = -1)

    double overlapThreshold = 0.4;
    SummaryTableUpdaterArgs *argsTemplate = SummaryTableUpdaterArgs_construct(iterator,
                                                                              copyIterator,
                                                                              resetIterator,
                                                                              destructIterator,
                                                                              getNextBlock,
                                                                              summaryTableList,
                                                                              binArray,
                                                                              CoverageInfo_overlapAnnotationIndex,
                                                                              -1,
                                                                              get_inference_prediction_label,
                                                                              get_inference_truth_label,
                                                                              isMetricOverlapBased,
                                                                              overlapThreshold);
    int threads = 4;
    // update all tables with multi-threading
    ptBlock_updateSummaryTableListForAllCategory1(argsTemplate, header->numberOfAnnotations, threads);


    bool correct = true;

    SummaryTable *summaryTable;
    // whole genome bin 1
    summaryTable = SummaryTableList_getTable(summaryTableList, 1, 0);
    correct &= Double_equality2DArray(wholeGenomeBin1, summaryTable->table, 4, 5);

    // whole genome bin 2
    summaryTable = SummaryTableList_getTable(summaryTableList, 1, 1);
    correct &= Double_equality2DArray(wholeGenomeBin2, summaryTable->table, 4, 5);

    // annotation 1 bin 1
    summaryTable = SummaryTableList_getTable(summaryTableList, 2, 0);
    correct &= Double_equality2DArray(annotation1Bin1, summaryTable->table, 4, 5);

    // annotation 1 bin 2
    summaryTable = SummaryTableList_getTable(summaryTableList, 2, 1);
    correct &= Double_equality2DArray(annotation1Bin2, summaryTable->table, 4, 5);

    // annotation 2 bin 1
    summaryTable = SummaryTableList_getTable(summaryTableList, 3, 0);
    correct &= Double_equality2DArray(annotation2Bin1, summaryTable->table, 4, 5);

    // annotation 2 bin 2
    summaryTable = SummaryTableList_getTable(summaryTableList, 3, 1);
    correct &= Double_equality2DArray(annotation2Bin2, summaryTable->table, 4, 5);

    SummaryTableList_destruct(summaryTableList);

    if (chunksCreator != NULL) {
        // free iterator
        ChunkIterator_destruct((ChunkIterator *) iterator);
        // free blocks/chunks
        ChunksCreator_destruct(chunksCreator);
    }
    if (blocksPerContig != NULL) {
        // free iterator
        ptBlockItrPerContig_construct((ptBlockItrPerContig *) iterator);
        // free blocks/chunks
        stHash_destruct(blocksPerContig);
        // free header
        CoverageHeader_destruct(header);
    }
    IntBinArray_destruct(binArray);
    Double_destruct2DArray(wholeGenomeBin1, 4);
    Double_destruct2DArray(wholeGenomeBin2, 4);
    Double_destruct2DArray(annotation1Bin1, 4);
    Double_destruct2DArray(annotation2Bin1, 4);
    Double_destruct2DArray(annotation1Bin2, 4);
    Double_destruct2DArray(annotation2Bin2, 4);
    SummaryTableUpdaterArgs_destruct(argsTemplate);


    return correct;
}

int main(int argc, char *argv[]) {

    bool allTestsPassed = true;

    // test 1
    bool test1Passed = test_SummaryTable_increment();
    printf("[summary_table] Test SummaryTable_increment:");
    printf(test1Passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");
    allTestsPassed &= test1Passed;

    // test 2
    bool test2Passed = test_SummaryTable_getRowString();
    printf("[summary_table] Test SummaryTable_getRowString:");
    printf(test2Passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");
    allTestsPassed &= test2Passed;

    // test 3
    bool test3Passed = test_SummaryTableList_increment();
    printf("[summary_table] Test SummaryTableList_increment:");
    printf(test3Passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");
    allTestsPassed &= test3Passed;

    // test 4
    bool test4Passed = test_SummaryTableList_getRowString();
    printf("[summary_table] Test SummaryTableList_getRowString:");
    printf(test4Passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");
    allTestsPassed &= test4Passed;


    // test 5
    bool test5Passed = test_ptBlock_updateSummaryTableListWithIterator("tests/test_files/summary_table/test_1.cov",
                                                                       "tests/test_files/summary_table/test_1_bin_array.txt",
                                                                       true,
                                                                       false);
    printf("[summary_table] Test updating SummaryTableList with ChunkIterator:");
    printf(test5Passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");
    allTestsPassed &= test5Passed;

    // test 6
    bool test6Passed = test_ptBlock_updateSummaryTableListWithIterator("tests/test_files/summary_table/test_1.cov",
                                                                       "tests/test_files/summary_table/test_1_bin_array.txt",
                                                                       false,
                                                                       false);
    printf("[summary_table] Test updating SummaryTableList with ptBlockItrPerContig:");
    printf(test6Passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");
    allTestsPassed &= test6Passed;

    // test 7
    bool test7Passed = test_ptBlock_updateSummaryTableListWithIterator("tests/test_files/summary_table/test_1.cov",
                                                                       "tests/test_files/summary_table/test_1_bin_array.txt",
                                                                       true,
                                                                       true);
    printf("[summary_table] Test updating SummaryTableList with ChunkIterator (overlap-based):");
    printf(test7Passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");
    allTestsPassed &= test7Passed;

    // test 8
    bool test8Passed = test_ptBlock_updateSummaryTableListWithIterator("tests/test_files/summary_table/test_1.cov",
                                                                       "tests/test_files/summary_table/test_1_bin_array.txt",
                                                                       false,
                                                                       true);
    printf("[summary_table] Test updating SummaryTableList with ptBlockItrPerContig (overlap-based):");
    printf(test8Passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");
    allTestsPassed &= test8Passed;

    if (allTestsPassed)
        return 0;
    else
        return 1;
}
