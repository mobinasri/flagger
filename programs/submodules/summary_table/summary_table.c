//
// Created by mobin on 5/14/24.
//

#include "summary_table.h"
#include <stdint.h>
#include <string.h>
#include <sys/stat.h>   // stat
#include <stdbool.h>    // bool type
#include <assert.h>
#include <stdio.h>
#include <time.h>
#include "sonLib.h"

SummaryTable *SummaryTable_construct(int numberOfRows, int numberOfColumns) {
    SummaryTable *summaryTable = malloc(sizeof(SummaryTable));
    summaryTable->table = Double_construct2DArray(numberOfRows, numberOfColumns);
    summaryTable->tablePercentage = Double_construct2DArray(numberOfRows, numberOfColumns);
    summaryTable->totalPerRow = Double_construct1DArray(numberOfRows);
    summaryTable->numberOfRows = numberOfRows;
    summaryTable->numberOfColumns = numberOfColumns;
    return summaryTable;
}

void SummaryTable_destruct(SummaryTable *summaryTable) {
    if (summaryTable->table != NULL) {
        Double_destruct2DArray(summaryTable->table, summaryTable->numberOfRows);
    }
    if (summaryTable->tablePercentage != NULL) {
        Double_destruct2DArray(summaryTable->tablePercentage, summaryTable->numberOfRows);
    }
    if (summaryTable->totalPerRow != NULL) {
        Double_destruct1DArray(summaryTable->totalPerRow);
    }
    free(summaryTable);
}

void SummaryTable_increment(SummaryTable *summaryTable, int rowIndex, int columnIndex, double value) {
    if (value < 0) {
        fprintf(stderr, "[%s] Warning: summary table is taking a negative value %.2e\n", get_timestamp(), value);
    }
    summaryTable->table[rowIndex][columnIndex] += value;
    summaryTable->totalPerRow[rowIndex] += value;
    for (int j = 0; j < summaryTable->numberOfColumns; j++) {
        if (0 < summaryTable->totalPerRow[rowIndex]) {
            summaryTable->tablePercentage[rowIndex][j] =
                    summaryTable->table[rowIndex][j] / summaryTable->totalPerRow[rowIndex] * 100.0;
        }
    }
}

char *SummaryTable_getRowString(SummaryTable *summaryTable, int rowIndex, char delimiter) {
    char* rowString = String_joinDoubleArrayWithFormat(summaryTable->table[rowIndex], summaryTable->numberOfColumns, delimiter,
                                            "%.2f");
    strcpy(summaryTable->rowString, rowString);
    free(rowString);
    return summaryTable->rowString;
}

char *SummaryTable_getRowStringPercentage(SummaryTable *summaryTable, int rowIndex, char delimiter) {
    char* rowString = String_joinDoubleArrayWithFormat(summaryTable->tablePercentage[rowIndex], summaryTable->numberOfColumns, delimiter,
                                            "%.2f");
    strcpy(summaryTable->rowString, rowString);
    free(rowString);
    return summaryTable->rowString;
}


int SummaryTableList_getTableIndex(SummaryTableList *summaryTableList, int catIndex1, int catIndex2) {
    return catIndex1 * summaryTableList->numberOfCategories1 + catIndex2;
}

SummaryTable *SummaryTableList_getTable(SummaryTableList *summaryTableList, int catIndex1, int catIndex2) {
    int tableIndex = SummaryTableList_getTableIndex(summaryTableList, catIndex1, catIndex2);
    return stList_get(summaryTableList->summaryTables, tableIndex);
}

void SummaryTableList_increment(SummaryTableList *summaryTableList,
                                int catIndex1,
                                int catIndex2,
                                int rowIndex,
                                int columnIndex,
                                double value) {
    SummaryTable *summaryTable = SummaryTableList_getTable(summaryTableList, catIndex1, catIndex2);
    SummaryTable_increment(summaryTable, rowIndex, columnIndex, value);
}

double SummaryTableList_getValue(SummaryTableList *summaryTableList,
                                 int catIndex1,
                                 int catIndex2,
                                 int rowIndex,
                                 int columnIndex) {
    SummaryTable *summaryTable = SummaryTableList_getTable(summaryTableList, catIndex1, catIndex2);
    return summaryTable->table[rowIndex][columnIndex];
}

double SummaryTableList_getValuePercentage(SummaryTableList *summaryTableList,
                                           int catIndex1,
                                           int catIndex2,
                                           int rowIndex,
                                           int columnIndex) {
    SummaryTable *summaryTable = SummaryTableList_getTable(summaryTableList, catIndex1, catIndex2);
    return summaryTable->tablePercentage[rowIndex][columnIndex];
}

SummaryTableList *SummaryTableList_construct(stList *categoryNames1,
                                             stList *categoryNames2,
                                             int numberOfRows,
                                             int numberOfColumns) {
    SummaryTableList *summaryTableList = malloc(sizeof(SummaryTableList));
    int numberOfCategories1 = stList_length(categoryNames1);
    int numberOfCategories2 = stList_length(categoryNames2);
    int totalNumberOfTables = numberOfCategories1 * numberOfCategories2;
    fprintf(stderr ,"#1\n");
    summaryTableList->numberOfCategories1 = numberOfCategories1;
    summaryTableList->numberOfCategories2 = numberOfCategories2;
    summaryTableList->totalNumberOfTables = totalNumberOfTables;
    summaryTableList->numberOfRows = numberOfRows;
    summaryTableList->numberOfColumns = numberOfColumns;
    // construct tables
    summaryTableList->summaryTables = stList_construct3(0, SummaryTable_destruct);
    for (int c1 = 0; c1 < numberOfCategories1; c1++) {
        for (int c2 = 0; c2 < numberOfCategories2; c2++) {
            int tableIndex = SummaryTableList_getTableIndex(summaryTableList, c1, c2);
            SummaryTable *summaryTable = SummaryTable_construct(numberOfRows, numberOfColumns);
            stList_append(summaryTableList->summaryTables, summaryTable);
        }
    }
    // copy category names 1
    summaryTableList->categoryNames1 = stList_construct3(0, free);
    for (int i = 0; i < stList_length(categoryNames1); i++) {
        stList_append(summaryTableList->categoryNames1, copyString(stList_get(categoryNames1, i)));
    }
    // copy category names 2
    summaryTableList->categoryNames2 = stList_construct3(0, free);
    for (int i = 0; i < stList_length(categoryNames2); i++) {
        stList_append(summaryTableList->categoryNames2, copyString(stList_get(categoryNames2, i)));
    }
    return summaryTableList;
}

char *SummaryTableList_getRowString(SummaryTableList *summaryTableList,
                                    int catIndex1,
                                    int catIndex2,
                                    int rowIndex,
                                    char delimiter) {
    SummaryTable *summaryTable = SummaryTableList_getTable(summaryTableList, catIndex1, catIndex2);
    return SummaryTable_getRowString(summaryTable, rowIndex, delimiter);
}

char *SummaryTableList_getRowStringPercentage(SummaryTableList *summaryTableList,
                                              int catIndex1,
                                              int catIndex2,
                                              int rowIndex,
                                              char delimiter) {
    SummaryTable *summaryTable = SummaryTableList_getTable(summaryTableList, catIndex1, catIndex2);
    return SummaryTable_getRowStringPercentage(summaryTable, rowIndex, delimiter);
}

void SummaryTableList_destruct(SummaryTableList *summaryTableList) {
    if (summaryTableList->summaryTables != NULL) stList_destruct(summaryTableList->summaryTables);
    if (summaryTableList->categoryNames1 != NULL) stList_destruct(summaryTableList->categoryNames1);
    if (summaryTableList->categoryNames2 != NULL) stList_destruct(summaryTableList->categoryNames2);
    free(summaryTableList);
}





