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

    const char *rowOneString = "30.00,20.00";
    const char *rowTwoString = "0.00,0.00";
    const char *rowOneStringPercentage = "60.00,40.00";
    const char *rowTwoStringPercentage = "0.00,0.00";

    correct &= strcmp(SummaryTable_getRowString(summaryTable, 0, ','), rowOneString) == 0;
    correct &= strcmp(SummaryTable_getRowString(summaryTable, 1, ','), rowTwoString) == 0;
    correct &= strcmp(SummaryTable_getRowStringPercentage(summaryTable, 0, ','), rowOneStringPercentage) == 0;
    correct &= strcmp(SummaryTable_getRowStringPercentage(summaryTable, 1, ','), rowTwoStringPercentage) == 0;

    SummaryTable_destruct(summaryTable);
    return correct;
}


int main(int argc, char *argv[]) {

    bool allTestsPassed = true;

    // test 1
    bool test1Passed = test_SummaryTable_increment();
    printf("[summary_table] Test test_SummaryTable_increment:");
    printf(test1Passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");
    allTestsPassed &= test1Passed;

    // test 2
    bool test2Passed = test_SummaryTable_getRowString();
    printf("[summary_table] Test test_SummaryTable_getRowString:");
    printf(test2Passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");
    allTestsPassed &= test2Passed;

    if (allTestsPassed)
        return 0;
    else
        return 1;
}
