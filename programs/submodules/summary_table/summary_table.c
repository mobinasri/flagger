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

#define HAP_INDEX 2

SummaryTable *SummaryTable_construct(int numberOfRows, int numberOfColumns) {
    SummaryTable *summaryTable = malloc(sizeof(SummaryTable));
    summaryTable->table = Double_construct2DArray(numberOfRows, numberOfColumns);
    summaryTable->tablePercentage = Double_construct2DArray(numberOfRows, numberOfColumns);
    summaryTable->totalPerRow = Double_construct1DArray(numberOfRows);
    summaryTable->totalPerRowPercentage = Double_construct1DArray(numberOfRows);
    summaryTable->totalSum = 0.0;
    summaryTable->numberOfRows = numberOfRows;
    summaryTable->numberOfColumns = numberOfColumns;
    summaryTable->rowNames = NULL;
    summaryTable->columnNames = NULL;
    return summaryTable;
}

SummaryTable *SummaryTable_constructByNames(stList *rowNames, stList *columnNames) {
    SummaryTable *summaryTable = malloc(sizeof(SummaryTable));
    int numberOfRows = stList_length(rowNames);
    int numberOfColumns = stList_length(columnNames);
    summaryTable->table = Double_construct2DArray(numberOfRows, numberOfColumns);
    summaryTable->tablePercentage = Double_construct2DArray(numberOfRows, numberOfColumns);
    summaryTable->totalPerRow = Double_construct1DArray(numberOfRows);
    summaryTable->totalPerRowPercentage = Double_construct1DArray(numberOfRows);
    summaryTable->totalSum = 0.0;
    summaryTable->numberOfRows = numberOfRows;
    summaryTable->numberOfColumns = numberOfColumns;
    summaryTable->rowNames = stList_copyStringList(rowNames);
    summaryTable->columnNames = stList_copyStringList(columnNames);
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
    if (summaryTable->totalPerRowPercentage != NULL) {
        Double_destruct1DArray(summaryTable->totalPerRowPercentage);
    }
    if (summaryTable->rowNames != NULL) {
        stList_destruct(summaryTable->rowNames);
    }
    if (summaryTable->columnNames != NULL) {
        stList_destruct(summaryTable->columnNames);
    }
    free(summaryTable);
}

void SummaryTable_increment(SummaryTable *summaryTable, int rowIndex, int columnIndex, double value) {
    if (value < 0) {
        fprintf(stderr, "[%s] Warning: summary table is taking a negative value %.2e\n", get_timestamp(), value);
    }
    summaryTable->table[rowIndex][columnIndex] += value;
    summaryTable->totalPerRow[rowIndex] += value;
    summaryTable->totalSum += value;
    // update tablePercentage
    for (int j = 0; j < summaryTable->numberOfColumns; j++) {
        if (0 < summaryTable->totalPerRow[rowIndex]) {
            summaryTable->tablePercentage[rowIndex][j] =
                    summaryTable->table[rowIndex][j] / summaryTable->totalPerRow[rowIndex] * 100.0;
        }
    }
    // update totalPerRowPercentage
    if (0 < summaryTable->totalSum) {
        for (int i = 0; i < summaryTable->numberOfRows; i++) {
            summaryTable->totalPerRowPercentage[i] = summaryTable->totalPerRow[i] / summaryTable->totalSum * 100.0;
        }
    }
}

double SummaryTable_getValue(SummaryTable *summaryTable, int rowIndex, int columnIndex) {
    return summaryTable->table[rowIndex][columnIndex];
}

double SummaryTable_getTPCountInRow(SummaryTable *summaryTable, int rowIndex) {
    return summaryTable->table[rowIndex][rowIndex];
}

double SummaryTable_getAllCountInRow(SummaryTable *summaryTable, int rowIndex) {
    return summaryTable->totalPerRow[rowIndex];
}


double SummaryTable_getAccuracyPercentage(SummaryTable *summaryTable, bool excludeLastRowAndColumn) {
    int tp = 0;
    int all = 0;
    for (int i = 0; i < summaryTable->numberOfRows; i++) {
        if (excludeLastRowAndColumn && i == summaryTable->numberOfRows - 1) break;
        tp += summaryTable->table[i][i];
        all += summaryTable->totalPerRow[i];
        // remove the value of last column if excludeLastRowAndColumn was true
        all -= excludeLastRowAndColumn ? summaryTable->table[i][summaryTable->numberOfColumns - 1] : 0;
    }
    if (all <= 0.0) {
        return 0.0;
    }
    double acc = (double) tp / all * 100.0;
    return acc;
}

double SummaryTable_getMicroAverageAcrossRowsPercentage(SummaryTable *summaryTable, bool excludeLastRow) {
    double tp = 0;
    double all = 0;
    for (int i = 0; i < summaryTable->numberOfRows; i++) {
        if (excludeLastRow && i == summaryTable->numberOfRows - 1) break;
        tp += summaryTable->table[i][i];
        all += summaryTable->totalPerRow[i];
    }
    if (all <= 0.0) {
        return 0.0;
    }
    double microAvg = tp / all * 100.0;
    return microAvg;
}

double SummaryTable_getMacroAverageAcrossRowsPercentage(SummaryTable *summaryTable, bool excludeLastRow) {
    double sum = 0.0;
    int n = 0;
    for (int i = 0; i < summaryTable->numberOfRows; i++) {
        if (excludeLastRow && i == summaryTable->numberOfRows - 1) break;
        if (0 < summaryTable->totalPerRow[i]) {
            sum += summaryTable->table[i][i] / summaryTable->totalPerRow[i];
            n += 1;
        }
    }
    double macroAvg = (double) sum / n * 100.0;
    return macroAvg;
}

char *SummaryTable_getRowName(SummaryTable *summaryTable, int rowIndex) {
    if (summaryTable->numberOfRows <= rowIndex) {
        fprintf(stderr, "[%s] Error: row index %d cannot be greater than %d.\n",
                get_timestamp(),
                rowIndex,
                summaryTable->numberOfRows - 1);
        exit(EXIT_FAILURE);
    }
    //reset row string
    summaryTable->rowName[0] = '\0';
    // add row name first if row names exist
    if (summaryTable->rowNames != NULL) {
        char *rowName = stList_get(summaryTable->rowNames, rowIndex);
        // copy name + delimiter
        sprintf(summaryTable->rowName, "%s", rowName);
    } else {
        sprintf(summaryTable->rowName, "%d", rowIndex);
    }
    return summaryTable->rowName;
}

char *SummaryTable_getRowString(SummaryTable *summaryTable, int rowIndex, char delimiter, bool addRowIndex) {
    if (summaryTable->numberOfRows <= rowIndex) {
        fprintf(stderr, "[%s] Error: row index %d cannot be greater than %d.\n",
                get_timestamp(),
                rowIndex,
                summaryTable->numberOfRows - 1);
        exit(EXIT_FAILURE);
    }
    //reset row string
    summaryTable->rowString[0] = '\0';
    // add row name first if row names exist
    if (summaryTable->rowNames != NULL) {
        char *rowName = stList_get(summaryTable->rowNames, rowIndex);
        // copy name + delimiter
        sprintf(summaryTable->rowString, "%s%c", rowName, delimiter);
    } else if (addRowIndex) {
        sprintf(summaryTable->rowString, "%d%c", rowIndex, delimiter);
    }
    char *rowString = String_joinDoubleArrayWithFormat(summaryTable->table[rowIndex],
                                                       summaryTable->numberOfColumns,
                                                       delimiter,
                                                       "%.2f");
    strcpy(summaryTable->rowString + strlen(summaryTable->rowString), rowString);
    free(rowString);
    return summaryTable->rowString;
}

char *SummaryTable_getRowStringPercentage(SummaryTable *summaryTable, int rowIndex, char delimiter, bool addRowIndex) {
    if (summaryTable->numberOfRows <= rowIndex) {
        fprintf(stderr, "[%s] Error: row index %d cannot be greater than %d.\n",
                get_timestamp(),
                rowIndex,
                summaryTable->numberOfRows - 1);
        exit(EXIT_FAILURE);
    }
    //reset row string
    summaryTable->rowString[0] = '\0';
    // add row name first if row names exist
    if (summaryTable->rowNames != NULL) {
        char *rowName = stList_get(summaryTable->rowNames, rowIndex);
        // copy name + delimiter
        sprintf(summaryTable->rowString, "%s%c", rowName, delimiter);
    } else if (addRowIndex) {
        sprintf(summaryTable->rowString, "%d%c", rowIndex, delimiter);
    }
    char *rowString = String_joinDoubleArrayWithFormat(summaryTable->tablePercentage[rowIndex],
                                                       summaryTable->numberOfColumns,
                                                       delimiter,
                                                       "%.2f");
    strcpy(summaryTable->rowString + strlen(summaryTable->rowString), rowString);
    free(rowString);
    return summaryTable->rowString;
}


char *SummaryTable_getTotalPerRowString(SummaryTable *summaryTable, char delimiter, const char *prefixEntry) {
    //reset row string
    summaryTable->rowString[0] = '\0';
    if (prefixEntry != NULL) {
        // copy entry + delimiter
        sprintf(summaryTable->rowString, "%s%c", prefixEntry, delimiter);
    }
    char *rowString = String_joinDoubleArrayWithFormat(summaryTable->totalPerRow,
                                                       summaryTable->numberOfRows,
                                                       delimiter,
                                                       "%.2f");
    strcpy(summaryTable->rowString + strlen(summaryTable->rowString), rowString);
    free(rowString);
    return summaryTable->rowString;
}

char *SummaryTable_getTotalPerRowStringPercentage(SummaryTable *summaryTable, char delimiter, const char *prefixEntry) {
    //reset row string
    summaryTable->rowString[0] = '\0';
    if (prefixEntry != NULL) {
        // copy entry + delimiter
        sprintf(summaryTable->rowString, "%s%c", prefixEntry, delimiter);
    }
    char *rowString = String_joinDoubleArrayWithFormat(summaryTable->totalPerRowPercentage,
                                                       summaryTable->numberOfRows,
                                                       delimiter,
                                                       "%.2f");
    strcpy(summaryTable->rowString + strlen(summaryTable->rowString), rowString);
    free(rowString);
    return summaryTable->rowString;
}

int SummaryTableList_getTableIndex(SummaryTableList *summaryTableList, int catIndex1, int catIndex2) {
    return catIndex1 * summaryTableList->numberOfCategories2 + catIndex2;
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
    summaryTableList->numberOfCategories1 = numberOfCategories1;
    summaryTableList->numberOfCategories2 = numberOfCategories2;
    summaryTableList->totalNumberOfTables = totalNumberOfTables;
    summaryTableList->numberOfRows = numberOfRows;
    summaryTableList->numberOfColumns = numberOfColumns;
    // construct tables
    summaryTableList->summaryTables = stList_construct3(summaryTableList->totalNumberOfTables, SummaryTable_destruct);
    for (int c1 = 0; c1 < numberOfCategories1; c1++) {
        for (int c2 = 0; c2 < numberOfCategories2; c2++) {
            int tableIndex = SummaryTableList_getTableIndex(summaryTableList, c1, c2);
            SummaryTable *summaryTable = SummaryTable_construct(numberOfRows, numberOfColumns);
            stList_set(summaryTableList->summaryTables, tableIndex, summaryTable);
        }
    }
    // copy category names 1
    summaryTableList->categoryNames1 = stList_copyStringList(categoryNames1);
    // copy category names 2
    summaryTableList->categoryNames2 = stList_copyStringList(categoryNames2);
    return summaryTableList;
}

SummaryTableList *SummaryTableList_constructByNames(stList *categoryNames1,
                                                    stList *categoryNames2,
                                                    stList *rowNames,
                                                    stList *columnNames) {
    SummaryTableList *summaryTableList = malloc(sizeof(SummaryTableList));
    int numberOfRows = stList_length(rowNames);
    int numberOfColumns = stList_length(columnNames);
    int numberOfCategories1 = stList_length(categoryNames1);
    int numberOfCategories2 = stList_length(categoryNames2);
    int totalNumberOfTables = numberOfCategories1 * numberOfCategories2;
    summaryTableList->numberOfCategories1 = numberOfCategories1;
    summaryTableList->numberOfCategories2 = numberOfCategories2;
    summaryTableList->totalNumberOfTables = totalNumberOfTables;
    summaryTableList->numberOfRows = numberOfRows;
    summaryTableList->numberOfColumns = numberOfColumns;
    // construct tables
    summaryTableList->summaryTables = stList_construct3(summaryTableList->totalNumberOfTables, SummaryTable_destruct);
    for (int c1 = 0; c1 < numberOfCategories1; c1++) {
        for (int c2 = 0; c2 < numberOfCategories2; c2++) {
            int tableIndex = SummaryTableList_getTableIndex(summaryTableList, c1, c2);
            SummaryTable *summaryTable = SummaryTable_constructByNames(rowNames, columnNames);
            stList_set(summaryTableList->summaryTables, tableIndex, summaryTable);
        }
    }
    // copy category names 1
    summaryTableList->categoryNames1 = stList_copyStringList(categoryNames1);
    // copy category names 2
    summaryTableList->categoryNames2 = stList_copyStringList(categoryNames2);
    return summaryTableList;
}

char *SummaryTableList_getRowString(SummaryTableList *summaryTableList,
                                    int catIndex1,
                                    int catIndex2,
                                    int rowIndex,
                                    char delimiter,
                                    bool addRowIndex) {
    SummaryTable *summaryTable = SummaryTableList_getTable(summaryTableList, catIndex1, catIndex2);
    return SummaryTable_getRowString(summaryTable, rowIndex, delimiter, addRowIndex);
}

char *SummaryTableList_getRowStringPercentage(SummaryTableList *summaryTableList,
                                              int catIndex1,
                                              int catIndex2,
                                              int rowIndex,
                                              char delimiter,
                                              bool addRowIndex) {
    SummaryTable *summaryTable = SummaryTableList_getTable(summaryTableList, catIndex1, catIndex2);
    return SummaryTable_getRowStringPercentage(summaryTable, rowIndex, delimiter, addRowIndex);
}

char *SummaryTableList_getTotalPerRowString(SummaryTableList *summaryTableList,
                                            int catIndex1,
                                            int catIndex2,
                                            char delimiter,
                                            const char *prefixEntry) {
    SummaryTable *summaryTable = SummaryTableList_getTable(summaryTableList, catIndex1, catIndex2);
    return SummaryTable_getTotalPerRowString(summaryTable, delimiter, prefixEntry);
}

char *SummaryTableList_getTotalPerRowStringPercentage(SummaryTableList *summaryTableList,
                                                      int catIndex1,
                                                      int catIndex2,
                                                      char delimiter,
                                                      const char *prefixEntry) {
    SummaryTable *summaryTable = SummaryTableList_getTable(summaryTableList, catIndex1, catIndex2);
    return SummaryTable_getTotalPerRowStringPercentage(summaryTable, delimiter, prefixEntry);
}

void SummaryTableList_writeIntoFile(SummaryTableList *summaryTableList, FILE *fp, const char *linePrefix) {
    for (int c1 = 0; c1 < summaryTableList->numberOfCategories1; c1++) {
        char *c1Name = stList_get(summaryTableList->categoryNames1, c1);
        for (int c2 = 0; c2 < summaryTableList->numberOfCategories2; c2++) {
            char *c2Name = stList_get(summaryTableList->categoryNames2, c2);
            for (int rowIndex = 0; rowIndex < summaryTableList->numberOfRows; rowIndex++) {
                char *tableRowString = SummaryTableList_getRowString(summaryTableList,
                                                                     c1,
                                                                     c2,
                                                                     rowIndex,
                                                                     '\t',
                                                                     true);
                fprintf(fp, "%s\t%s\t%s\t%s\n", linePrefix, c1Name, c2Name, tableRowString);
            }
        }
    }
}

void SummaryTableList_writePercentageIntoFile(SummaryTableList *summaryTableList, FILE *fp, const char *linePrefix) {
    for (int c1 = 0; c1 < summaryTableList->numberOfCategories1; c1++) {
        char *c1Name = stList_get(summaryTableList->categoryNames1, c1);
        for (int c2 = 0; c2 < summaryTableList->numberOfCategories2; c2++) {
            char *c2Name = stList_get(summaryTableList->categoryNames2, c2);
            for (int rowIndex = 0; rowIndex < summaryTableList->numberOfRows; rowIndex++) {
                char *tableRowString = SummaryTableList_getRowStringPercentage(summaryTableList,
                                                                               c1,
                                                                               c2,
                                                                               rowIndex,
                                                                               '\t',
                                                                               true);
                fprintf(fp, "%s\t%s\t%s\t%s\n", linePrefix, c1Name, c2Name, tableRowString);
            }
        }
    }
}

void SummaryTableList_writeTotalPerRowIntoFile(SummaryTableList *summaryTableList, FILE *fp, const char *linePrefix) {
    for (int c1 = 0; c1 < summaryTableList->numberOfCategories1; c1++) {
        char *c1Name = stList_get(summaryTableList->categoryNames1, c1);
        for (int c2 = 0; c2 < summaryTableList->numberOfCategories2; c2++) {
            char *c2Name = stList_get(summaryTableList->categoryNames2, c2);
            char *tableRowString = SummaryTableList_getTotalPerRowString(summaryTableList,
                                                                         c1,
                                                                         c2,
                                                                         '\t',
                                                                         "ALL_LABELS");
            fprintf(fp, "%s\t%s\t%s\t%s\n", linePrefix, c1Name, c2Name, tableRowString);
        }
    }
}

void SummaryTableList_writeTotalPerRowPercentageIntoFile(SummaryTableList *summaryTableList, FILE *fp,
                                                         const char *linePrefix) {
    for (int c1 = 0; c1 < summaryTableList->numberOfCategories1; c1++) {
        char *c1Name = stList_get(summaryTableList->categoryNames1, c1);
        for (int c2 = 0; c2 < summaryTableList->numberOfCategories2; c2++) {
            char *c2Name = stList_get(summaryTableList->categoryNames2, c2);
            char *tableRowString = SummaryTableList_getTotalPerRowStringPercentage(summaryTableList,
                                                                                   c1,
                                                                                   c2,
                                                                                   '\t',
                                                                                   "ALL_LABELS");
            fprintf(fp, "%s\t%s\t%s\t%s\n", linePrefix, c1Name, c2Name, tableRowString);
        }
    }

}

void SummaryTableList_writeFinalStatisticsIntoFile(SummaryTableList *recallTableList,
                                                   SummaryTableList *precisionTableList,
                                                   FILE *fout,
                                                   const char *linePrefix) {
    // last row is reserved for unknown labels
    int numberOfLabels = recallTableList->numberOfRows - 1;

    // iterate over category 1 (annotation/region indices)
    for (int c1 = 0; c1 < recallTableList->numberOfCategories1; c1++) {
        char *c1Name = stList_get(recallTableList->categoryNames1, c1);
        // iterate over category 2 (size bins)
        for (int c2 = 0; c2 < recallTableList->numberOfCategories2; c2++) {
            char *c2Name = stList_get(recallTableList->categoryNames2, c2);
            // get recall table
            SummaryTable *recallTable = SummaryTableList_getTable(recallTableList, c1, c2);
            // get precision table
            SummaryTable *precisionTable = SummaryTableList_getTable(precisionTableList, c1, c2);
            double totalTpRecall = 0.0;
            double totalTpPrecision = 0.0;
            double totalRecall = 0.0;
            double totalPrecision = 0.0;

            // macro mean
            double sumOfRecall = 0.0;
            double sumOfPrecision = 0.0;
            // macro mean no hap
            double sumOfRecallNoHap = 0.0;
            double sumOfPrecisionNoHap = 0.0;

            // harmonic mean
            double sumOfReciprocalRecall = 0.0;
            double sumOfReciprocalPrecision = 0.0;
            // harmonic mean no hap
            double sumOfReciprocalRecallNoHap = 0.0;
            double sumOfReciprocalPrecisionNoHap = 0.0;

            // number of states with non-zero denominator
            int numberOfNonZeroDenomRecall = 0;
            int numberOfNonZeroDenomPrecision = 0;
            // no hap
            int numberOfNonZeroDenomRecallNoHap = 0;
            int numberOfNonZeroDenomPrecisionNoHap = 0;
            // iterate over row indices
            for (int rowIndex = 0; rowIndex < numberOfLabels; rowIndex++) {
                // recall is used as suffix for the variables related to the table whose reference label is truth
                // precision is used as suffix for the variables related to the table whose reference label is prediction
                double tpRecall = SummaryTable_getTPCountInRow(recallTable, rowIndex);
                totalTpRecall += tpRecall;
                double tpPrecision = SummaryTable_getTPCountInRow(precisionTable, rowIndex);
                totalTpPrecision += tpPrecision;
                double fn = SummaryTable_getAllCountInRow(recallTable, rowIndex) - tpRecall;
                double fp = SummaryTable_getAllCountInRow(precisionTable, rowIndex) - tpPrecision;
                totalRecall += tpRecall + fn;
                totalPrecision += tpPrecision + fp;
                char *rowName = SummaryTable_getRowName(recallTable, rowIndex);
                double recallPercent = tpRecall / (tpRecall + fn + 1.0e-9) * 100.0;
                double precisionPercent = tpPrecision / (tpPrecision + fp + 1.0e-9) * 100.0;
                numberOfNonZeroDenomRecall += 1e-9 < (tpRecall + fn) ? 1 : 0;
                numberOfNonZeroDenomPrecision += 1e-9 < (tpPrecision + fp) ? 1 : 0;

                if(rowIndex != HAP_INDEX){ // hap index is 2
                    numberOfNonZeroDenomRecallNoHap += 1e-9 < (tpRecall + fn) ? 1 : 0;
                    numberOfNonZeroDenomPrecisionNoHap += 1e-9 < (tpPrecision + fp) ? 1 : 0;
                }
                // for calculating macro-average later
                sumOfRecall += recallPercent;
                sumOfPrecision += precisionPercent;
                if(rowIndex != HAP_INDEX) { // hap index is 2
                    sumOfRecallNoHap += recallPercent;
                    sumOfPrecisionNoHap += precisionPercent;
                }

                // for calculating harmonic mean later
                // 1.0e9 is sufficiently large
                if (1e-9 < (tpRecall + fn)) {
                    sumOfReciprocalRecall += 0.0 < recallPercent ? 1.0 / recallPercent : 1.0e9;
                    if(rowIndex != HAP_INDEX) { // hap index is 2
                        sumOfReciprocalRecallNoHap += 0.0 < recallPercent ? 1.0 / recallPercent : 1.0e9;
                    }
                }
                if (1e-9 < (tpPrecision + fp)) {
                    sumOfReciprocalPrecision += 0.0 < precisionPercent ? 1.0 / precisionPercent : 1.0e9;
                    if(rowIndex != HAP_INDEX) { // hap index is 2
                        sumOfReciprocalPrecisionNoHap += 0.0 < precisionPercent ? 1.0 / precisionPercent : 1.0e9;
                    }
                }
                double f1Score = 2 * precisionPercent * recallPercent / (precisionPercent + recallPercent + 1.0e-9);

                char recallPercentStr[20];
                if (1e-9 < (tpRecall + fn)) {
                    sprintf(recallPercentStr, "%.2f", recallPercent);
                } else {
                    sprintf(recallPercentStr, "NA");
                }

                char precisionPercentStr[20];
                if (1e-9 < (tpPrecision + fp)) {
                    sprintf(precisionPercentStr, "%.2f", precisionPercent);
                } else {
                    sprintf(precisionPercentStr, "NA");
                }

                char f1ScoreStr[20];
                if (1e-9 < (tpRecall + fn) && 1e-9 < (tpPrecision + fp)) {
                    sprintf(f1ScoreStr, "%.2f", f1Score);
                } else {
                    sprintf(f1ScoreStr, "NA");
                }
                // TP_Prediction_Ref
                // TP_Truth_Ref
                // FP
                // FN
                // Total_Prediction_Ref
                // Total_Truth_Ref
                // Precision
                // Recall
                // F1-Score
                // Accuracy_Prediction_Ref
                // Accuracy_Truth_Ref
                fprintf(fout, "%s\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%s\t%s\t%s\t%s\t%s\n",
                        linePrefix,
                        c1Name, c2Name,
                        rowName,
                        tpPrecision, tpRecall, fp, fn,
                        tpPrecision + fp, tpRecall + fn,
                        precisionPercentStr, recallPercentStr, f1ScoreStr, "NA", "NA");
            }
            double macroAverageRecall = 0.0;
            double macroAveragePrecision = 0.0;
            double macroAverageRecallNoHap = 0.0;
            double macroAveragePrecisionNoHap = 0.0;

            char macroAverageRecallStr[20];
            if (0 < numberOfNonZeroDenomRecall) {
                macroAverageRecall = sumOfRecall / numberOfNonZeroDenomRecall;
                sprintf(macroAverageRecallStr, "%.2f", macroAverageRecall);
            } else {
                sprintf(macroAverageRecallStr, "NA");
            }

            char macroAveragePrecisionStr[20];
            if (0 < numberOfNonZeroDenomPrecision) {
                macroAveragePrecision = sumOfPrecision / numberOfNonZeroDenomPrecision;
                sprintf(macroAveragePrecisionStr, "%.2f", macroAveragePrecision);
            } else {
                sprintf(macroAveragePrecisionStr, "NA");
            }

            char macroAverageRecallStrNoHap[20];
            if (0 < numberOfNonZeroDenomRecallNoHap) {
                macroAverageRecallNoHap = sumOfRecallNoHap / numberOfNonZeroDenomRecallNoHap;
                sprintf(macroAverageRecallStrNoHap, "%.2f", macroAverageRecallNoHap);
            } else {
                sprintf(macroAverageRecallStrNoHap, "NA");
            }

            char macroAveragePrecisionStrNoHap[20];
            if (0 < numberOfNonZeroDenomPrecisionNoHap) {
                macroAveragePrecisionNoHap = sumOfPrecisionNoHap / numberOfNonZeroDenomPrecisionNoHap;
                sprintf(macroAveragePrecisionStrNoHap, "%.2f", macroAveragePrecisionNoHap);
            } else {
                sprintf(macroAveragePrecisionStrNoHap, "NA");
            }

            double harmonicMeanRecall = 0.0;
            double harmonicMeanPrecision = 0.0;
            double harmonicMeanRecallNoHap = 0.0;
            double harmonicMeanPrecisionNoHap = 0.0;

            char harmonicMeanRecallStr[20];
            if (0 < numberOfNonZeroDenomRecall) {
                harmonicMeanRecall = (double) numberOfNonZeroDenomRecall / sumOfReciprocalRecall;
                sprintf(harmonicMeanRecallStr, "%.2f", harmonicMeanRecall);
            } else {
                sprintf(harmonicMeanRecallStr, "NA");
            }

            char harmonicMeanPrecisionStr[20];
            if (0 < numberOfNonZeroDenomPrecision) {
                harmonicMeanPrecision = (double) numberOfNonZeroDenomPrecision / sumOfReciprocalPrecision;
                sprintf(harmonicMeanPrecisionStr, "%.2f", harmonicMeanPrecision);
            } else {
                sprintf(harmonicMeanPrecisionStr, "NA");
            }

            char harmonicMeanRecallStrNoHap[20];
            if (0 < numberOfNonZeroDenomRecallNoHap) {
                harmonicMeanRecallNoHap = (double) numberOfNonZeroDenomRecallNoHap / sumOfReciprocalRecallNoHap;
                sprintf(harmonicMeanRecallStrNoHap, "%.2f", harmonicMeanRecallNoHap);
            } else {
                sprintf(harmonicMeanRecallStrNoHap, "NA");
            }

            char harmonicMeanPrecisionStrNoHap[20];
            if (0 < numberOfNonZeroDenomPrecisionNoHap) {
                harmonicMeanPrecisionNoHap = (double) numberOfNonZeroDenomPrecisionNoHap / sumOfReciprocalPrecisionNoHap;
                sprintf(harmonicMeanPrecisionStrNoHap, "%.2f", harmonicMeanPrecisionNoHap);
            } else {
                sprintf(harmonicMeanPrecisionStrNoHap, "NA");
            }

            char macroF1ScoreStr[20];
            if (0 < numberOfNonZeroDenomRecall && 0 < numberOfNonZeroDenomPrecision) {
                double macroF1Score = 2 * macroAverageRecall * macroAveragePrecision /
                                      (macroAverageRecall + macroAveragePrecision + 1.0e-9);
                sprintf(macroF1ScoreStr, "%.2f", macroF1Score);
            } else {
                sprintf(macroF1ScoreStr, "NA");
            }

            char macroF1ScoreStrNoHap[20];
            if (0 < numberOfNonZeroDenomRecallNoHap && 0 < numberOfNonZeroDenomPrecisionNoHap) {
                double macroF1ScoreNoHap = 2 * macroAverageRecallNoHap * macroAveragePrecisionNoHap /
                                      (macroAverageRecallNoHap + macroAveragePrecisionNoHap + 1.0e-9);
                sprintf(macroF1ScoreStrNoHap, "%.2f", macroF1ScoreNoHap);
            } else {
                sprintf(macroF1ScoreStrNoHap, "NA");
            }

            char harmonicF1ScoreStr[20];
            if (0 < numberOfNonZeroDenomRecall && 0 < numberOfNonZeroDenomPrecision) {
                double harmonicF1Score = 2 * harmonicMeanRecall * harmonicMeanPrecision /
                                         (harmonicMeanRecall + harmonicMeanPrecision + 1.0e-9);
                sprintf(harmonicF1ScoreStr, "%.2f", harmonicF1Score);
            } else {
                sprintf(harmonicF1ScoreStr, "NA");
            }

            char harmonicF1ScoreStrNoHap[20];
            if (0 < numberOfNonZeroDenomRecallNoHap && 0 < numberOfNonZeroDenomPrecisionNoHap) {
                double harmonicF1ScoreNoHap = 2 * harmonicMeanRecallNoHap * harmonicMeanPrecisionNoHap /
                                         (harmonicMeanRecallNoHap + harmonicMeanPrecisionNoHap + 1.0e-9);
                sprintf(harmonicF1ScoreStrNoHap, "%.2f", harmonicF1ScoreNoHap);
            } else {
                sprintf(harmonicF1ScoreStrNoHap, "NA");
            }

            fprintf(fout, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                    linePrefix,
                    c1Name, c2Name,
                    "MACRO_AVERAGE",
                    "NA", "NA", "NA", "NA",
                    "NA", "NA",
                    macroAveragePrecisionStr, macroAverageRecallStr, macroF1ScoreStr, "NA", "NA");

            fprintf(fout, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                    linePrefix,
                    c1Name, c2Name,
                    "MACRO_AVERAGE_NO_HAP",
                    "NA", "NA", "NA", "NA",
                    "NA", "NA",
                    macroAveragePrecisionStrNoHap, macroAverageRecallStrNoHap, macroF1ScoreStrNoHap, "NA", "NA");

            fprintf(fout, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                    linePrefix,
                    c1Name, c2Name,
                    "HARMONIC_MEAN",
                    "NA", "NA", "NA", "NA",
                    "NA", "NA",
                    harmonicMeanPrecisionStr, harmonicMeanRecallStr, harmonicF1ScoreStr, "NA", "NA");

            fprintf(fout, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                    linePrefix,
                    c1Name, c2Name,
                    "HARMONIC_MEAN_NO_HAP",
                    "NA", "NA", "NA", "NA",
                    "NA", "NA",
                    harmonicMeanPrecisionStrNoHap, harmonicMeanRecallStrNoHap, harmonicF1ScoreStrNoHap, "NA", "NA");

            double accuracyPredictionRef = totalTpPrecision / (totalPrecision + 1.0e-9) * 100.0;
            double accuracyTruthRef = totalTpRecall / (totalRecall + 1e-9) * 100.0;
            fprintf(fout, "%s\t%s\t%s\t%s\t%.2f\t%.2f\t%s\t%s\t%.2f\t%.2f\t%s\t%s\t%s\t%.2f\t%.2f\n",
                    linePrefix,
                    c1Name, c2Name,
                    "ACCURACY",
                    totalTpPrecision, totalTpRecall, "NA", "NA",
                    totalPrecision, totalRecall,
                    "NA", "NA", "NA", accuracyPredictionRef, accuracyTruthRef);
        }
    }
}


void SummaryTableList_writeFinalAunStatisticsIntoFile(SummaryTableList *numeratorTableList,
                                                      SummaryTableList *denominatorTableList,
                                                      FILE *fout,
                                                      const char *linePrefix) {
    // last row is reserved for unknown labels
    int numberOfLabels = numeratorTableList->numberOfRows - 1;

    // iterate over category 1 (annotation/region indices)
    for (int c1 = 0; c1 < numeratorTableList->numberOfCategories1; c1++) {
        char *c1Name = stList_get(numeratorTableList->categoryNames1, c1);
        // iterate over category 2 (size bins)
        for (int c2 = 0; c2 < numeratorTableList->numberOfCategories2; c2++) {
            char *c2Name = stList_get(numeratorTableList->categoryNames2, c2);
            // get numerator table
            SummaryTable *numeratorTable = SummaryTableList_getTable(numeratorTableList, c1, c2);
            // get denominator table
            SummaryTable *denominatorTable = SummaryTableList_getTable(denominatorTableList, c1, c2);
            double aunSum = 0.0;
            double aunSumReciprocal = 0.0;
            int numberOfNonZeroDenom = 0;
            // iterate over row indices
            for (int rowIndex = 0; rowIndex < numberOfLabels; rowIndex++) {
                double aunDenom = SummaryTable_getValue(denominatorTable, rowIndex, rowIndex);
                double aunNum = SummaryTable_getValue(numeratorTable, rowIndex, rowIndex);
                double aun = aunNum / (aunDenom + 1e-9);
                numberOfNonZeroDenom += 0 < aunDenom ? 1 : 0;
                aunSum += aun;
                if (0 < aunDenom) {
                    //1.0e9 is sufficiently large
                    aunSumReciprocal += 0.0 < aun ? 1.0 / aun : 1.0e9;
                }
                char *rowName = SummaryTable_getRowName(numeratorTable, rowIndex);
                fprintf(fout, "%s\t%s\t%s\t%s\t%.2f\n",
                        linePrefix,
                        c1Name, c2Name,
                        rowName,
                        aun);
            }
            double aunAvg = 0.0;
            char aunAvgStr[20];
            if (0 < numberOfNonZeroDenom) {
                aunAvg = aunSum / numberOfNonZeroDenom;
                sprintf(aunAvgStr, "%.2f", aunAvg);
            } else {
                sprintf(aunAvgStr, "NA");
            }

            double aunHarmonicMean = 0.0;
            char aunHarmonicMeanStr[20];
            if (0 < numberOfNonZeroDenom) {
                aunHarmonicMean = (double) numberOfNonZeroDenom / aunSumReciprocal;
                sprintf(aunHarmonicMeanStr, "%.2f", aunHarmonicMean);
            } else {
                sprintf(aunHarmonicMeanStr, "NA");
            }

            fprintf(fout, "%s\t%s\t%s\t%s\t%s\n",
                    linePrefix,
                    c1Name, c2Name,
                    "AVERAGE",
                    aunAvgStr);

            fprintf(fout, "%s\t%s\t%s\t%s\t%s\n",
                    linePrefix,
                    c1Name, c2Name,
                    "HARMONIC_MEAN",
                    aunHarmonicMeanStr);
        }
    }
}


void SummaryTableList_destruct(SummaryTableList *summaryTableList) {
    if (summaryTableList->summaryTables != NULL) stList_destruct(summaryTableList->summaryTables);
    if (summaryTableList->categoryNames1 != NULL) stList_destruct(summaryTableList->categoryNames1);
    if (summaryTableList->categoryNames2 != NULL) stList_destruct(summaryTableList->categoryNames2);
    free(summaryTableList);
}


// It will modify confusion row in place
void convertBaseLevelToOverlapBased(double *refLabelConfusionRow,
                                    int columnSize,
                                    int refLabelBlockLength,
                                    double overlapThreshold) {
    bool atLeastOneHit = false;
    for (int i = 0; i < columnSize; i++) {
        double overlapRatio = (double) refLabelConfusionRow[i] / refLabelBlockLength;
        if (overlapThreshold < overlapRatio) atLeastOneHit = true;
        refLabelConfusionRow[i] = overlapThreshold < overlapRatio ? 1 : 0;
    }
    // if no hit was found set the last column to 1
    // last column is reserved for not defined labels
    // it should rarely happen that no hit is found
    if (atLeastOneHit == false) {
        refLabelConfusionRow[columnSize - 1] = 1;
    }
}

SummaryTableUpdaterArgs *SummaryTableUpdaterArgs_construct(void *blockIterator,
                                                           void *(*copyBlockIterator)(void *),
                                                           void *(*resetBlockIterator)(void *),
                                                           void (*destructBlockIterator)(void *),
                                                           ptBlock *(*getNextBlock)(void *, char *),
                                                           SummaryTableList *summaryTableList,
                                                           IntBinArray *sizeBinArray,
                                                           int (*overlapFuncCategoryIndex1)(CoverageInfo *, int),
                                                           int categoryIndex1,
                                                           int8_t (*getRefLabelFunction)(Inference *),
                                                           int8_t (*getQueryLabelFunction)(Inference *),
                                                           MetricType metricType,
                                                           double overlapThreshold,
                                                           SummaryTableList *auxiliarySummaryTableList) {
    SummaryTableUpdaterArgs *args = (SummaryTableUpdaterArgs *) malloc(1 * sizeof(SummaryTableUpdaterArgs));
    args->blockIterator = blockIterator;
    args->copyBlockIterator = copyBlockIterator;
    args->resetBlockIterator = resetBlockIterator;
    args->destructBlockIterator = destructBlockIterator;
    args->getNextBlock = getNextBlock;
    args->summaryTableList = summaryTableList;
    args->sizeBinArray = sizeBinArray;
    args->overlapFuncCategoryIndex1 = overlapFuncCategoryIndex1;
    args->categoryIndex1 = categoryIndex1;
    args->getRefLabelFunction = getRefLabelFunction;
    args->getQueryLabelFunction = getQueryLabelFunction;
    args->metricType = metricType;
    args->overlapThreshold = overlapThreshold;
    args->auxiliarySummaryTableList = auxiliarySummaryTableList;
    return args;
}

SummaryTableUpdaterArgs *SummaryTableUpdaterArgs_copy(SummaryTableUpdaterArgs *src) {
    SummaryTableUpdaterArgs *dest = (SummaryTableUpdaterArgs *) malloc(1 * sizeof(SummaryTableUpdaterArgs));
    dest->blockIterator = src->blockIterator;
    dest->copyBlockIterator = src->copyBlockIterator;
    dest->resetBlockIterator = src->resetBlockIterator;
    dest->destructBlockIterator = src->destructBlockIterator;
    dest->getNextBlock = src->getNextBlock;
    dest->summaryTableList = src->summaryTableList;
    dest->sizeBinArray = src->sizeBinArray;
    dest->overlapFuncCategoryIndex1 = src->overlapFuncCategoryIndex1;
    dest->categoryIndex1 = src->categoryIndex1;
    dest->getRefLabelFunction = src->getRefLabelFunction;
    dest->getQueryLabelFunction = src->getQueryLabelFunction;
    dest->metricType = src->metricType;
    dest->overlapThreshold = src->overlapThreshold;
    dest->auxiliarySummaryTableList = src->auxiliarySummaryTableList;
    return dest;
}

void SummaryTableUpdaterArgs_destruct(SummaryTableUpdaterArgs *args) {
    free(args);
}


void
SummaryTableList_addCreationJobsForAllCategory1(SummaryTableUpdaterArgs *argsTemplate, int sizeOfCategory1,
                                                tpool_t *threadPool) {
    for (int categoryIndex1 = 0; categoryIndex1 < sizeOfCategory1; categoryIndex1++) {
        // make a copy of args
        SummaryTableUpdaterArgs *argsToRun = SummaryTableUpdaterArgs_copy(argsTemplate);
        // set the category index 1 whose related summary table is going to be updated
        argsToRun->categoryIndex1 = categoryIndex1;
        // create arg struct for tpool
        work_arg_t *argWork = malloc(sizeof(work_arg_t));
        argWork->data = (void *) argsToRun;
        // Add a new job to the thread pool
        tpool_add_work(threadPool,
                       SummaryTableList_updateByUpdaterArgsForThreadPool,
                       (void *) argWork);
        //fprintf(stderr, "[%s] Created thread for updating summary table for category1 index %d (out of range [0-%d])\n",
        //        get_timestamp(), categoryIndex1, sizeOfCategory1 - 1);
    }
}

void SummaryTableList_updateByUpdaterArgsForThreadPool(void *argWork_) {
    // get the arguments
    work_arg_t *argWork = argWork_;
    SummaryTableUpdaterArgs *args = argWork->data;
    SummaryTableList_updateByUpdaterArgs(args);
    SummaryTableUpdaterArgs_destruct(args);
}


// blockIterator can be created by either
// ChunkIterator_construct or ptBlockItrPerContig_construct
void SummaryTableList_updateByUpdaterArgs(SummaryTableUpdaterArgs *args) {
    // fetch the arguments
    void *(*copyBlockIterator)(void *) = args->copyBlockIterator;
    void (*resetBlockIterator)(void *) = args->resetBlockIterator;
    void (*destructBlockIterator)(void *) = args->destructBlockIterator;
    ptBlock *(*getNextBlock)(void *, char *) = args->getNextBlock;
    SummaryTableList *summaryTableList = args->summaryTableList;
    // in the current implementation aux table can only be the table that contains total lengths
    // of truth labels. It is required for computing AuN ratio values.
    SummaryTableList *auxTableList = args->auxiliarySummaryTableList;
    IntBinArray *sizeBinArray = args->sizeBinArray;
    int (*overlapFuncCategoryIndex1)(CoverageInfo *, int) = args->overlapFuncCategoryIndex1;
    int categoryIndex1 = args->categoryIndex1;
    int8_t(*getRefLabelFunction)(Inference * ) = args->getRefLabelFunction;
    int8_t(*getQueryLabelFunction)(Inference * ) = args->getQueryLabelFunction;
    MetricType metricType = args->metricType;
    double overlapThreshold = args->overlapThreshold;

    // will be used only for computing AuN metric
    stList **queryLengthsPerLabel = (stList **) malloc(summaryTableList->numberOfColumns * sizeof(stList * ));
    for (int labelIndex = 0; labelIndex < summaryTableList->numberOfColumns; labelIndex++) {
        queryLengthsPerLabel[labelIndex] = stList_construct3(0, free);
    }

    // make a copy of iterator and reset it
    void *blockIterator = copyBlockIterator(args->blockIterator);
    resetBlockIterator(blockIterator);

    if (summaryTableList->numberOfCategories1 <= categoryIndex1) {
        fprintf(stderr,
                "[%s] Error: annotation index (%d) for updating category 1 is greater than or equal to the size of the category (%d).",
                get_timestamp(),
                categoryIndex1,
                summaryTableList->numberOfCategories1);
        exit(EXIT_FAILURE);
    }

    // +1 for undefined query labels
    double *refLabelConfusionRow = Double_construct1DArray(summaryTableList->numberOfColumns);
    Double_fill1DArray(refLabelConfusionRow, summaryTableList->numberOfColumns, 0.0);

    int refLabelStart = -1;
    int queryLabelStart = -1;
    int preRefLabel = -1;
    int preQueryLabel = -1;
    int preBlockEnd = -1;
    CoverageInfo *preCoverageInfo = NULL;

    char ctg[200];
    char preCtg[200];
    preCtg[0] = '\0';

    ptBlock *block = NULL;
    int *queryLabelBlockLenPtr = NULL;
    int queryLabelBlockLen = 0;

    while ((block = getNextBlock(blockIterator, ctg)) != NULL) {

        // get coverage info for this block
        CoverageInfo *coverageInfo = (CoverageInfo *) block->data;
        int start = block->rfs;
        int end = block->rfe;

        // get labels
        if (coverageInfo->data == NULL) {
            fprintf(stderr, "[%s] Warning: inference data does not exist for updating summary tables.\n",
                    get_timestamp());
        }
        Inference *inference = coverageInfo->data;
        int refLabel = getRefLabelFunction(inference);
        // if refLabel is -1 change it to numberOfRows - 1 since last row is for "Unk"
        refLabel = refLabel == -1 ? summaryTableList->numberOfRows - 1 : refLabel;
        int queryLabel = getQueryLabelFunction(inference);
        // if queryLabel is -1 change it to numberOfColumns - 1 since last column is for "Unk"
        queryLabel = queryLabel == -1 ? summaryTableList->numberOfColumns - 1 : queryLabel;

        // set event flags
        bool contigChanged = (preCtg[0] != '\0') && (strcmp(preCtg, ctg) != 0);
        bool refLabelChanged = refLabel != preRefLabel;
        bool queryLabelChanged = queryLabel != preQueryLabel;

        bool annotationInCurrent = overlapFuncCategoryIndex1(coverageInfo, categoryIndex1);
        bool annotationInPrevious = overlapFuncCategoryIndex1(preCoverageInfo, categoryIndex1);
        bool annotationContinued = annotationInCurrent && annotationInPrevious;
        bool annotationStarted = annotationInCurrent && !annotationInPrevious;
        bool annotationEnded = !annotationInCurrent && annotationInPrevious;

        bool preRefLabelIsValid = preRefLabel != -1;
        bool preQueryLabelIsValid = preQueryLabel != -1;
        /*
        fprintf(stderr, "%s %d-%d annot = %d \n", ctg, start, end, categoryIndex1);
        fprintf(stderr, "\t rlabel = %d\n", refLabel);
        fprintf(stderr, "\t qlabel = %d\n", queryLabel);
        fprintf(stderr, "\t contigChanged = %d\n", contigChanged);
        fprintf(stderr, "\t refLabelChanged = %d\n", refLabelChanged);
        fprintf(stderr, "\t annotationInCurrent = %d\n", annotationInCurrent);
        fprintf(stderr, "\t annotationInPrevious = %d\n", annotationInPrevious);
        fprintf(stderr, "\t annotationContinued = %d\n", annotationContinued);
        fprintf(stderr, "\t annotationStarted = %d\n", annotationStarted);
        fprintf(stderr, "\t annotationEnded = %d\n", annotationEnded);
        fprintf(stderr, "\t preRefLabelIsValid = %d\n", preRefLabelIsValid);
         */

        // one annotation-ref-label block has ended
        // update summary table
        if (preRefLabelIsValid && (
                (annotationContinued && refLabelChanged) ||
                (annotationInPrevious && contigChanged) ||
                annotationEnded)
                ) {
            int refLabelBlockLen = preBlockEnd - refLabelStart + 1;
            int binIndicesLength = 0;
            int *binIndices = IntBinArray_getBinIndices(sizeBinArray, refLabelBlockLen, &binIndicesLength);
            if (metricType == METRIC_OVERLAP_BASED) {
                convertBaseLevelToOverlapBased(refLabelConfusionRow,
                                               summaryTableList->numberOfColumns,
                                               refLabelBlockLen,
                                               overlapThreshold);
            }

            if (metricType == METRIC_AUN) {
                if (preQueryLabelIsValid) {
                    // add last block (with a contiguous ref and query label)
                    queryLabelBlockLenPtr = (int *) malloc(sizeof(int));
                    *queryLabelBlockLenPtr = preBlockEnd - queryLabelStart + 1;
                    stList_append(queryLengthsPerLabel[preQueryLabel], queryLabelBlockLenPtr);
                }
                //fprintf(stderr, "new ref block\n");
                // iterating over query labels
                for (int q = 0; q < summaryTableList->numberOfColumns; q++) {
                    for (int queryBlockIndex = 0;
                         queryBlockIndex < stList_length(queryLengthsPerLabel[q]);
                         queryBlockIndex++) {
                        queryLabelBlockLenPtr = (int *) stList_get(queryLengthsPerLabel[q], queryBlockIndex);
                        queryLabelBlockLen = (int) *queryLabelBlockLenPtr;
                        //fprintf(stderr, "annot_idx=%d\tref=%d\tquery=%d\trefBlockLen=%d\tqueryBlockLen=%d\n", categoryIndex1, preRefLabel, q, refLabelBlockLen, queryLabelBlockLen);
                        // update ref label confusion row for computing AuN ratio values
                        refLabelConfusionRow[q] += (double) queryLabelBlockLen * queryLabelBlockLen;
                    }
                }
            }


            // iterating over size bin indices
            for (int bi = 0; bi < binIndicesLength; bi++) {
                int binIndex = binIndices[bi];
                // each bin index as its own total ref length
                double totalSizeOfPreRefLabel = metricType == METRIC_AUN ? SummaryTableList_getValue(auxTableList,
                                                                                                     categoryIndex1,
                                                                                                     binIndex,
                                                                                                     preRefLabel,
                                                                                                     preRefLabel) : 1.0;
                // iterating over query labels
                for (int q = 0; q < summaryTableList->numberOfColumns; q++) {
                    int catIndex1 = categoryIndex1;
                    int catIndex2 = binIndex;
                    int rowIndex = preRefLabel;
                    int columnIndex = q;
                    //fprintf(stderr, "block len = %d, [%d][%d] [%d][%d] += %d\n", blockLen, catIndex1, catIndex2, rowIndex, columnIndex, refLabelConfusionRow[q]);
                    SummaryTableList_increment(summaryTableList,
                                               catIndex1,
                                               catIndex2,
                                               rowIndex,
                                               columnIndex,
                                               refLabelConfusionRow[q] / totalSizeOfPreRefLabel);
                }
            }
            free(binIndices);
        }

        // query label changed within a ref block
        if (annotationInCurrent &&
            metricType == METRIC_AUN &&
            preQueryLabelIsValid &&
            queryLabelChanged &&
            (annotationContinued && !refLabelChanged) &&
            !contigChanged) {
            queryLabelBlockLenPtr = (int *) malloc(sizeof(int));
            *queryLabelBlockLenPtr = preBlockEnd - queryLabelStart + 1;
            stList_append(queryLengthsPerLabel[preQueryLabel], queryLabelBlockLenPtr);
        }


        // reset confusion row and update start location
        if ((!annotationInCurrent && contigChanged) ||
            annotationEnded) {
            refLabelStart = -1;
            queryLabelStart = -1;
            Double_fill1DArray(refLabelConfusionRow, summaryTableList->numberOfColumns, 0);
        }

        // reset confusion row and update start location
        if ((annotationContinued && refLabelChanged) ||
            (annotationInCurrent && contigChanged) ||
            annotationStarted) {
            refLabelStart = start;
            Double_fill1DArray(refLabelConfusionRow, summaryTableList->numberOfColumns, 0);
            for (int q = 0; q < summaryTableList->numberOfColumns; q++) {
                // reset the list of query block lengths
                stList_destruct(queryLengthsPerLabel[q]);
                queryLengthsPerLabel[q] = stList_construct3(0, free);
            }
        }

        // update start location for query label
        if ((annotationContinued && refLabelChanged) ||
            (annotationContinued && queryLabelChanged) ||
            (annotationInCurrent && contigChanged) ||
            annotationStarted) {
            queryLabelStart = start;
        }

        // update confusion row
        if (annotationInCurrent && metricType != METRIC_AUN) {
            refLabelConfusionRow[queryLabel] += end - start + 1;
        }


        preCoverageInfo = coverageInfo;
        preRefLabel = refLabel;
        preQueryLabel = queryLabel;
        strcpy(preCtg, ctg);
        preBlockEnd = end;
    }

    bool preRefLabelIsValid = preRefLabel != -1;
    bool preQueryLabelIsValid = preQueryLabel != -1;
    // check last window and update summary tables if it had overlap with annotation
    bool annotationInLastWindow = overlapFuncCategoryIndex1(preCoverageInfo, categoryIndex1);
    if (annotationInLastWindow && preRefLabelIsValid) {
        int refLabelBlockLen = preBlockEnd - refLabelStart + 1;
        int binIndicesLength = 0;
        int *binIndices = IntBinArray_getBinIndices(sizeBinArray, refLabelBlockLen, &binIndicesLength);
        if (metricType == METRIC_OVERLAP_BASED) {
            convertBaseLevelToOverlapBased(refLabelConfusionRow,
                                           summaryTableList->numberOfColumns,
                                           refLabelBlockLen,
                                           overlapThreshold);
        }
        if (metricType == METRIC_AUN) {
            if (preQueryLabelIsValid) {
                // add last block (with a contiguous ref and query label)
                queryLabelBlockLenPtr = (int *) malloc(sizeof(int));
                *queryLabelBlockLenPtr = preBlockEnd - queryLabelStart + 1;
                stList_append(queryLengthsPerLabel[preQueryLabel], queryLabelBlockLenPtr);
            }
            //fprintf(stderr, "new ref block final\n");
            // iterating over query labels
            for (int q = 0; q < summaryTableList->numberOfColumns; q++) {
                for (int queryBlockIndex = 0;
                     queryBlockIndex < stList_length(queryLengthsPerLabel[q]);
                     queryBlockIndex++) {
                    queryLabelBlockLenPtr = (int *) stList_get(queryLengthsPerLabel[q], queryBlockIndex);
                    queryLabelBlockLen = (int) *queryLabelBlockLenPtr;
                    //fprintf(stderr, "annot_idx=%d\tref=%d\tquery=%d\trefBlockLen=%d\tqueryBlockLen=%d\n", categoryIndex1, preRefLabel, q, refLabelBlockLen, queryLabelBlockLen);
                    // update ref label confusion row for computing AuN ratio values
                    refLabelConfusionRow[q] += (double) queryLabelBlockLen * queryLabelBlockLen;
                }
            }
        }
        // iterating over size bin indices
        for (int bi = 0; bi < binIndicesLength; bi++) {
            int binIndex = binIndices[bi];
            // each bin index as its own total ref length
            double totalSizeOfPreRefLabel = metricType == METRIC_AUN ? SummaryTableList_getValue(auxTableList,
                                                                                                 categoryIndex1,
                                                                                                 binIndex,
                                                                                                 preRefLabel,
                                                                                                 preRefLabel) : 1.0;

            // iterating over query labels
            for (int q = 0; q < summaryTableList->numberOfColumns; q++) {
                int catIndex1 = categoryIndex1;
                int catIndex2 = binIndex;
                int rowIndex = preRefLabel;
                int columnIndex = q;
                //fprintf(stderr, "[%d][%d] [%d][%d] += %d\n", catIndex1, catIndex2, rowIndex, columnIndex, refLabelConfusionRow[q]);
                SummaryTableList_increment(summaryTableList,
                                           catIndex1,
                                           catIndex2,
                                           rowIndex,
                                           columnIndex,
                                           refLabelConfusionRow[q] / totalSizeOfPreRefLabel);
            }
        }
        free(binIndices);
    }

    Double_destruct1DArray(refLabelConfusionRow);
    destructBlockIterator(blockIterator);
    for (int labelIndex = 0; labelIndex < summaryTableList->numberOfColumns; labelIndex++) {
        stList_destruct(queryLengthsPerLabel[labelIndex]);
    }
    free(queryLengthsPerLabel);
}

void SummaryTableList_constructAndFillByIterator(void *blockIterator,
                                                 BlockIteratorType blockIteratorType,
                                                 stList *categoryNames,
                                                 CategoryType categoryType,
                                                 IntBinArray *sizeBinArray,
                                                 MetricType metricType,
                                                 double overlapRatioThreshold,
                                                 int numberOfLabelsWithUnknown,
                                                 stList *labelNamesWithUnknown,
                                                 ComparisonType comparisonType,
                                                 SummaryTableList *auxiliarySummaryTableList,
                                                 SummaryTableListFullCatalog *catalog,
                                                 tpool_t *threadPool) {

    ptBlock *(*getNextBlock)(void *, char *);
    void *(*copyIterator)(void *);
    void (*destructIterator)(void *);
    void (*resetIterator)(void *);
    if (blockIteratorType == ITERATOR_BY_CHUNK) {
        copyIterator = ChunkIterator_copy;
        destructIterator = ChunkIterator_destruct;
        resetIterator = ChunkIterator_reset;
        getNextBlock = ChunkIterator_getNextPtBlock;
    } else if (blockIteratorType == ITERATOR_BY_COV_BLOCK) {
        copyIterator = ptBlockItrPerContig_copy;
        destructIterator = ptBlockItrPerContig_destruct;
        resetIterator = ptBlockItrPerContig_reset;
        getNextBlock = ptBlockItrPerContig_next;
    } else {
        fprintf(stderr, "[%s] block iterator type is not valid.\n", get_timestamp());
        exit(EXIT_FAILURE);
    }


    bool (*overlapFuncCategoryIndex1)(CoverageInfo *, int);
    if (categoryType == CATEGORY_ANNOTATION) {
        overlapFuncCategoryIndex1 = CoverageInfo_overlapAnnotationIndex;
    } else if (categoryType == CATEGORY_REGION) {
        overlapFuncCategoryIndex1 = CoverageInfo_overlapRegionIndex;
    } else {
        fprintf(stderr, "[%s] category type is not valid.\n", get_timestamp());
        exit(EXIT_FAILURE);
    }

    int8_t(*getRefLabelFunction)(Inference * );
    if (comparisonType == COMPARISON_TRUTH_VS_PREDICTION || comparisonType == COMPARISON_TRUTH_VS_TRUTH) {
        getRefLabelFunction = get_inference_truth_label;
    } else {
        getRefLabelFunction = get_inference_prediction_label;
    }


    int8_t(*getQueryLabelFunction)(Inference * );
    if (comparisonType == COMPARISON_TRUTH_VS_PREDICTION || comparisonType == COMPARISON_PREDICTION_VS_PREDICTION) {
        getQueryLabelFunction = get_inference_prediction_label;
    } else {
        getQueryLabelFunction = get_inference_truth_label;
    }

    // create summary table list
    stList *categoryNames1 = categoryNames;
    stList *categoryNames2 = sizeBinArray->names;

    int sizeOfCategory1 = stList_length(categoryNames1);
    int numberOfRows = numberOfLabelsWithUnknown;
    int numberOfColumns = numberOfLabelsWithUnknown;
    SummaryTableList *summaryTableList;
    if (labelNamesWithUnknown == NULL) {
        summaryTableList = SummaryTableList_construct(categoryNames1,
                                                      categoryNames2,
                                                      numberOfRows,
                                                      numberOfColumns);
    } else {
        summaryTableList = SummaryTableList_constructByNames(categoryNames1,
                                                             categoryNames2,
                                                             labelNamesWithUnknown,
                                                             labelNamesWithUnknown);
    }
    // insert summary table list to the catalog
    SummaryTableListFullCatalog_update(catalog, summaryTableList, categoryType, metricType, comparisonType);

    // create a template of update args with category 1 index set to -1
    SummaryTableUpdaterArgs *argsTemplate = SummaryTableUpdaterArgs_construct((void *) blockIterator,
                                                                              copyIterator,
                                                                              resetIterator,
                                                                              destructIterator,
                                                                              getNextBlock,
                                                                              summaryTableList,
                                                                              sizeBinArray,
                                                                              overlapFuncCategoryIndex1,
                                                                              -1,
                                                                              getRefLabelFunction,
                                                                              getQueryLabelFunction,
                                                                              metricType,
                                                                              overlapRatioThreshold,
                                                                              auxiliarySummaryTableList);
    // update all tables with multi-threading
    SummaryTableList_addCreationJobsForAllCategory1(argsTemplate, sizeOfCategory1, threadPool);

    SummaryTableUpdaterArgs_destruct(argsTemplate);
}


SummaryTableListFullCatalog *
SummaryTableListFullCatalog_constructNull(int dimCategoryType, int dimMetricType, int dimComparisonType) {
    SummaryTableListFullCatalog *catalog = malloc(sizeof(SummaryTableListFullCatalog));
    // make a 3D array of summary table list with the dimension of dimCategoryType x dimMetricType x dimComparisonType
    catalog->array = (SummaryTableList ****) malloc(dimCategoryType * sizeof(SummaryTableList ***));
    for (int categoryType = 0; categoryType < dimCategoryType; categoryType++) {
        // 1st dimension
        catalog->array[categoryType] = (SummaryTableList ***) malloc(dimMetricType * sizeof(SummaryTableList **));
        for (int metricType = 0; metricType < dimMetricType; metricType++) {
            // 2nd dimension
            catalog->array[categoryType][metricType] = (SummaryTableList **) malloc(dimComparisonType * sizeof(SummaryTableList *));
            // 3rd dimension
            for (int comparisonType = 0; comparisonType < dimComparisonType; comparisonType++) {
                catalog->array[categoryType][metricType][comparisonType] = NULL;
            } // end loop for category type
        }// end loop for comparison type
    }// end loop for metric type
    catalog->dimMetricType = dimMetricType;
    catalog->dimCategoryType = dimCategoryType;
    catalog->dimComparisonType = dimComparisonType;
    return catalog;
}

void SummaryTableListFullCatalog_update(SummaryTableListFullCatalog *catalog,
                                        SummaryTableList *summaryTableList,
                                        CategoryType categoryType,
                                        MetricType metricType,
                                        ComparisonType comparisonType) {
    // destruct old summary list
    if (catalog->array[categoryType][metricType][comparisonType] != NULL) {
        SummaryTableList_destruct(catalog->array[categoryType][metricType][comparisonType]);
    }
    // insert new one
    catalog->array[categoryType][metricType][comparisonType] = summaryTableList;
}

SummaryTableList *SummaryTableListFullCatalog_get(SummaryTableListFullCatalog *catalog,
                                                  CategoryType categoryType,
                                                  MetricType metricType,
                                                  ComparisonType comparisonType) {
    return catalog->array[categoryType][metricType][comparisonType];
}

void SummaryTableListFullCatalog_destruct(SummaryTableListFullCatalog *catalog) {
    for (int metricType = 0; metricType < catalog->dimMetricType; metricType++) {
        for (int comparisonType = 0; comparisonType < catalog->dimComparisonType; comparisonType++) {
            for (int categoryType = 0; categoryType < catalog->dimCategoryType; categoryType++) {
                if (catalog->array[categoryType][metricType][comparisonType] != NULL) {
                    SummaryTableList_destruct(catalog->array[categoryType][metricType][comparisonType]);
                }
            } // end loop for category type
        }// end loop for comparison type
    }// end loop for metric type
    free(catalog);
}

void SummaryTableListFullCatalog_write(SummaryTableListFullCatalog *catalog,
                                       CoverageHeader *header,
                                       stList *labelNamesWithUnknown,
                                       char *outputPath){

    // open file for writing summary tables
    FILE *fout = fopen(outputPath, "w");
    if (fout == NULL) {
        fprintf(stderr, "[%s] Error: %s cannot be opened.\n", get_timestamp(), outputPath);
        exit(EXIT_FAILURE);
    }

    // write column names in the header line
    char linePrefix[1000];
    sprintf(linePrefix,
            "#Statistic\tMetric_Type\tEntry_Type\tCategory_Type\tCategory_Name\tSize_Bin_Name\tRef_Label");
    int numberOfLabelsWithUnknown = stList_length(labelNamesWithUnknown);
    for (int i = 0; i < numberOfLabelsWithUnknown; i++) {
        // use label names if they are given
        if (labelNamesWithUnknown != NULL) {
            char *labelName = stList_get(labelNamesWithUnknown, i);
            sprintf(linePrefix + strlen(linePrefix), "\t%s", labelName);
        } else { // use label indices
            if (i == numberOfLabelsWithUnknown - 1) { // last index is reserved for unknown label
                sprintf(linePrefix + strlen(linePrefix), "\tlabel_unk");
            } else {
                sprintf(linePrefix + strlen(linePrefix), "\tlabel_%d", i);
            }
        }
    }
    fprintf(fout, "%s\n", linePrefix);

    //reset line
    linePrefix[0] = '\0';


    // open files for writing prediction vs truth stats (like precision and recall) if both labels are available

    char outputPathFinalStats[1000];
    char outputPathFinalAunStats[1000];
    char prefix[1000];
    memcpy(prefix, outputPath, strlen(outputPath) - 4); // ".tsv" has 4 characters
    prefix[strlen(outputPath) - 4] = '\0'; // ".tsv" has 4 characters
    sprintf(outputPathFinalStats, "%s.benchmarking.tsv", prefix);
    sprintf(outputPathFinalAunStats, "%s.benchmarking.auN_ratio.tsv", prefix);

    FILE *foutFinalStats = NULL;
    FILE *foutFinalAunStats = NULL;
    if (header->isTruthAvailable && header->isPredictionAvailable) {
        // tsv file for dumping precision/recall/accuracy stats
        foutFinalStats = fopen(outputPathFinalStats, "w");
        if (foutFinalStats == NULL) {
            fprintf(stderr,
                    "[%s] Error: %s cannot be opened for making a tsv file with final benchmarking stats.\n",
                    get_timestamp(), outputPathFinalStats);
            exit(EXIT_FAILURE);
        }
        fprintf(foutFinalStats,
                "#Metric_Type\tCategory_Type\tCategory_Name\tSize_Bin_Name\tLabel\tTP_Prediction_Ref\tTP_Truth_Ref\tFP\tFN\tTotal_Prediction_Ref\tTotal_Truth_Ref\tPrecision\tRecall\tF1-Score\tAccuracy_Prediction_Ref\tAccuracy_Truth_Ref\n");

        // tsv file for dumping aun stats
        foutFinalAunStats = fopen(outputPathFinalAunStats, "w");
        if (foutFinalAunStats == NULL) {
            fprintf(stderr,
                    "[%s] Error: %s cannot be opened for making a tsv file with benchmarking auN ratio values.\n",
                    get_timestamp(), outputPathFinalAunStats);
            exit(EXIT_FAILURE);
        }
        fprintf(foutFinalAunStats,
                "#Category_Type\tCategory_Name\tSize_Bin_Name\tLabel\tauN_Ratio\n");

    }


    for (int categoryType = 0; categoryType < catalog->dimCategoryType; categoryType++) {
        for (int metricType = 0; metricType < catalog->dimMetricType; metricType++) {
            for (int comparisonType = 0; comparisonType < catalog->dimComparisonType; comparisonType++) {
                SummaryTableList *summaryTableList = SummaryTableListFullCatalog_get(catalog,
                                                                                     categoryType,
                                                                                     metricType,
                                                                                     comparisonType);
                if (summaryTableList == NULL) continue;

                // write count values
                sprintf(linePrefix, "%s\t%s\tcount\t%s",
                        ComparisonTypeToString[comparisonType],
                        MetricTypeToString[metricType],
                        CategoryTypeToString[categoryType]);
                // for these comparisons a single row is enough
                if (comparisonType == COMPARISON_TRUTH_VS_TRUTH ||
                    comparisonType == COMPARISON_PREDICTION_VS_PREDICTION) {
                    SummaryTableList_writeTotalPerRowIntoFile(summaryTableList, fout, linePrefix);
                } else {
                    SummaryTableList_writeIntoFile(summaryTableList, fout, linePrefix);
                }


                // write percentages
                sprintf(linePrefix, "%s\t%s\tpercentage\t%s",
                        ComparisonTypeToString[comparisonType],
                        MetricTypeToString[metricType],
                        CategoryTypeToString[categoryType]);
                // for these comparisons a single row is enough
                if (comparisonType == COMPARISON_TRUTH_VS_TRUTH ||
                    comparisonType == COMPARISON_PREDICTION_VS_PREDICTION) {
                    SummaryTableList_writeTotalPerRowPercentageIntoFile(summaryTableList, fout, linePrefix);
                } else {
                    SummaryTableList_writePercentageIntoFile(summaryTableList, fout, linePrefix);
                }

                fprintf(stderr,
                        "[%s] Writing summary tables is done for categoryType = %s, metricType = %s, comparisonType = %s .\n",
                        get_timestamp(),
                        CategoryTypeToString[categoryType],
                        MetricTypeToString[metricType],
                        ComparisonTypeToString[comparisonType]);

            }// end comparison loop


            //
            // METRIC_BASE_LEVEL and METRIC_OVERLAP_BASED :
            //
            // if both labels are present for printing benchmarking stats
            // precision and recall are not defined for METRIC_AUN
            if (header->isTruthAvailable && header->isPredictionAvailable && metricType != METRIC_AUN) {
                SummaryTableList *precisionTables = SummaryTableListFullCatalog_get(catalog,
				                                                    categoryType,
				                                                    metricType,
				                                                    COMPARISON_PREDICTION_VS_TRUTH);
                SummaryTableList *recallTables = SummaryTableListFullCatalog_get(catalog,
                                                                                    categoryType,
                                                                                    metricType,
										    COMPARISON_TRUTH_VS_PREDICTION);

                // write percentages
                sprintf(linePrefix, "%s\t%s",
                        MetricTypeToString[metricType],
                        CategoryTypeToString[categoryType]);
                SummaryTableList_writeFinalStatisticsIntoFile(recallTables,
                                                              precisionTables,
                                                              foutFinalStats,
                                                              linePrefix);
                fprintf(stderr,
                        "[%s] Writing benchmarking stats (precision/recall/f1score) is done for categoryType = %s, metricType = %s .\n",
                        get_timestamp(),
                        CategoryTypeToString[categoryType],
                        MetricTypeToString[metricType]);

            }//end if

        } // end metric loop

        //
        // METRIC_AUN
        //
        // if both labels are present for printing benchmarking auN stats
        if (header->isTruthAvailable && header->isPredictionAvailable) {
            // auN ratio is only defined once we have truth labels as reference
            SummaryTableList *aunDenominatorTables = SummaryTableListFullCatalog_get(catalog,
                                                                                     categoryType,
                                                                                     METRIC_AUN,
                                                                                     COMPARISON_TRUTH_VS_TRUTH);
            SummaryTableList *aunNumeratorTables = SummaryTableListFullCatalog_get(catalog,
                                                                                   categoryType,
                                                                                   METRIC_AUN,
                                                                                   COMPARISON_TRUTH_VS_PREDICTION);

            // write auN ratio values
            sprintf(linePrefix,"%s",
                    CategoryTypeToString[categoryType]);
            SummaryTableList_writeFinalAunStatisticsIntoFile(aunNumeratorTables,
                                                             aunDenominatorTables,
                                                             foutFinalAunStats,
                                                             linePrefix);
        } // end if

    } // end category loop


    fprintf(stderr,
            "[%s] Writing tables to file %s is done.\n",
            get_timestamp(), outputPath);
    if (foutFinalStats != NULL) {
        fprintf(stderr,
                "[%s] Writing tables to file %s is done.\n",
                get_timestamp(), outputPathFinalStats);
    }
    if (foutFinalAunStats != NULL) {
        fprintf(stderr,
                "[%s] Writing tables to file %s is done.\n",
                get_timestamp(), outputPathFinalAunStats);
    }

    // close files
    fclose(fout);
    if (foutFinalStats != NULL) {
        fclose(foutFinalStats);
    }
    if (foutFinalAunStats != NULL) {
        fclose(foutFinalAunStats);
    }
}


void SummaryTableList_addCreationJobsForOneMetricType(CoverageHeader *header,
                                                      void *blockIterator,
                                                      BlockIteratorType blockIteratorType,
                                                      IntBinArray *sizeBinArray,
                                                      MetricType metricType,
                                                      double overlapRatioThreshold,
                                                      int numberOfLabelsWithUnknown,
                                                      stList *labelNamesWithUnknown,
                                                      SummaryTableListFullCatalog *catalog,
                                                      tpool_t *threadPool) {
    SummaryTableList *auxiliarySummaryTableList = NULL;

    // iterating over category types; annotation and region
    for (int categoryType = 0; categoryType < NUMBER_OF_CATEGORY_TYPES; categoryType++) {
        stList *categoryNames =
                categoryType == CATEGORY_REGION ? header->regionNames : header->annotationNames;
        // iterate over comparison types such as precision and recall
        for (int comparisonType = 0; comparisonType < NUMBER_OF_COMPARISON_TYPES; comparisonType++) {
            bool truthLabelIsNeeded = comparisonType == COMPARISON_TRUTH_VS_PREDICTION ||
                                      comparisonType == COMPARISON_PREDICTION_VS_TRUTH ||
                                      comparisonType == COMPARISON_TRUTH_VS_TRUTH;
            bool predictionLabelIsNeeded = comparisonType == COMPARISON_TRUTH_VS_PREDICTION ||
                                           comparisonType == COMPARISON_PREDICTION_VS_TRUTH ||
                                           comparisonType == COMPARISON_PREDICTION_VS_PREDICTION;
            if (header->isTruthAvailable == false && truthLabelIsNeeded) continue;
            if (header->isPredictionAvailable == false && predictionLabelIsNeeded) continue;

            // auN ratio stats are not defined for prediction as reference
            // so skip these comparisons
            if (metricType == METRIC_AUN &&
                (comparisonType == COMPARISON_PREDICTION_VS_PREDICTION ||
                 comparisonType == COMPARISON_PREDICTION_VS_TRUTH)) {
                continue;
            }

	    // now create jobs for METRIC_AUN
            // get auxiliary table necessary for getting auN metric values
            auxiliarySummaryTableList = metricType == METRIC_AUN ? SummaryTableListFullCatalog_get(catalog,
                                                                    categoryType,
                                                                    METRIC_BASE_LEVEL,
                                                                    COMPARISON_TRUTH_VS_TRUTH) : NULL;
        if (metricType == METRIC_AUN && auxiliarySummaryTableList == NULL) {
            fprintf(stderr,
                    "[%s] Error: For creating summary tables for the metric type of %s it is required to have the tables for the metric type of %s constructed and accessible beforehand.\n ",
                    get_timestamp(),
                    MetricTypeToString[METRIC_AUN],
                    MetricTypeToString[METRIC_BASE_LEVEL]);
            exit(EXIT_FAILURE);
        }

            fprintf(stderr,
                    "[%s] Initiating jobs for creating summary tables for categoryType = %s, metricType = %s, comparisonType = %s .\n",
                    get_timestamp(),
                    CategoryTypeToString[categoryType],
                    MetricTypeToString[metricType],
                    ComparisonTypeToString[comparisonType]);
            SummaryTableList_constructAndFillByIterator(blockIterator,
                                                        blockIteratorType,
                                                        categoryNames,
                                                        categoryType,
                                                        sizeBinArray,
                                                        metricType,
                                                        overlapRatioThreshold,
                                                        numberOfLabelsWithUnknown,
                                                        labelNamesWithUnknown,
                                                        comparisonType,
                                                        auxiliarySummaryTableList,
                                                        catalog,
                                                        threadPool);
        }
    }
}

void SummaryTableList_createAndWriteAllTables(void *iterator,
                                              BlockIteratorType blockIteratorType,
                                              CoverageHeader *header,
                                              const char *outputPath,
                                              const char *binArrayFilePath,
                                              stList *labelNamesWithUnknown,
                                              double overlapRatioThreshold,
                                              int threads) {

    IntBinArray *binArray;
    if (binArrayFilePath != NULL) {
        // parse bin intervals
        binArray = IntBinArray_constructFromFile(binArrayFilePath);
    } else {
        binArray = IntBinArray_constructSingleBin(0, 1e9, "ALL_SIZES");
    }


    // The list of label names has an additional name "Unk" for handling labels with the value of -1
    if (labelNamesWithUnknown != NULL && stList_length(labelNamesWithUnknown) - 1 != header->numberOfLabels) {
        fprintf(stderr,
                "[%s] Error: Number of label names %d  does not match the number of labels in the header %d.\n",
                get_timestamp(),
                stList_length(labelNamesWithUnknown) - 1,
                header->numberOfLabels);
        exit(EXIT_FAILURE);
    }


    int numberOfLabelsWithUnknown = header->numberOfLabels + 1;


    // make a 3D array of summary table list with the dimension of
    // NUMBER_OF_METRIC_TYPES x NUMBER_OF_COMPARISON_TYPES x NUMBER_OF_CATEGORY_TYPES
    SummaryTableListFullCatalog *catalog = SummaryTableListFullCatalog_constructNull(NUMBER_OF_CATEGORY_TYPES,
                                                                                     NUMBER_OF_METRIC_TYPES,
                                                                                     NUMBER_OF_COMPARISON_TYPES);

    // create a thread pool to parallelize table list creation
    tpool_t *threadPool = tpool_create(threads);
    if (header->isTruthAvailable || header->isPredictionAvailable) {
        // first add jobs are submitted except the ones for METRIC_AUN since
        // this metric is dependent on METRIC_BASE_LEVEL
        for (int metricType = 0; metricType < NUMBER_OF_METRIC_TYPES; metricType++) {
            if (metricType == METRIC_AUN) continue;
            // submit table creation jobs
            SummaryTableList_addCreationJobsForOneMetricType(header,
                                                             iterator,
                                                             blockIteratorType,
                                                             binArray,
                                                             metricType,
                                                             overlapRatioThreshold,
                                                             numberOfLabelsWithUnknown,
                                                             labelNamesWithUnknown,
                                                             catalog,
                                                             threadPool);
        } // end loop over metric types
	
	// wait until all jobs are done
        tpool_wait(threadPool);

        //submit table creation jobs for METRIC_AUN
        SummaryTableList_addCreationJobsForOneMetricType(header,
                                                         iterator,
                                                         blockIteratorType,
                                                         binArray,
                                                         METRIC_AUN,
                                                         overlapRatioThreshold,
                                                         numberOfLabelsWithUnknown,
                                                         labelNamesWithUnknown,
                                                         catalog,
                                                         threadPool);
    }
    // wait until all jobs are done
    tpool_wait(threadPool);
    tpool_destroy(threadPool);


    SummaryTableListFullCatalog_write(catalog, header, labelNamesWithUnknown, outputPath);

    // destruct bin array
    IntBinArray_destruct(binArray);
    // free full catalog
    SummaryTableListFullCatalog_destruct(catalog);
}
