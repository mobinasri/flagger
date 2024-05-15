#include "common.h"
#include "summary_table.h"
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
    stList *categoryNames1 = stList_construct3(0,free);
    stList_append(categoryNames1, copyString("cat1_0"));
    stList_append(categoryNames1, copyString("cat1_1"));
    stList *categoryNames2 = stList_construct3(0,free);
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
    stList *categoryNames1 = stList_construct3(0,free);
    stList_append(categoryNames1, copyString("cat1_0"));
    stList_append(categoryNames1, copyString("cat1_1"));
    stList *categoryNames2 = stList_construct3(0,free);
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

    if (allTestsPassed)
        return 0;
    else
        return 1;
}
