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
    correct &= summaryTablePercentage->table[0][0] == 60;
    correct &= summaryTablePercentage->table[0][1] == 40;
    correct &= summaryTablePercentage->table[1][0] == 0;
    correct &= summaryTablePercentage->table[1][1] == 0;

    // get total per row
    correct &= summaryTablePercentage->totalPerRow[0] == 50;
    correct &= summaryTablePercentage->totalPerRow[1] == 0;


    SummaryTable_destruct(summaryTable);
    return correct;
}


int main(int argc, char *argv[]) {

    bool allTestsPassed = true;

    // test 1
    bool test1Passed = test_SummaryTable_increment(test1Path);
    printf("[summary_table] Test test_SummaryTable_increment:");
    printf(test1Passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");
    allTestsPassed &= test1Passed;

    if (allTestsPassed)
        return 0;
    else
        return 1;
}
