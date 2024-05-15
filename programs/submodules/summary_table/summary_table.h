
#ifndef FLAGGER_SUMMARY_TABLE_H
#define FLAGGER_SUMMARY_TABLE_H

#include <stdint.h>
#include <string.h>
#include <sys/stat.h>   // stat
#include <stdbool.h>    // bool type
#include <assert.h>
#include <stdio.h>
#include <time.h>
#include "sonLib.h"



/*! @typedef
 * @abstract Structure for keeping a summary table
 */
typedef struct SummaryTable {
    double **table;
    double **tablePercentage;
    double *totalPerRow;
    int numberOfRows;
    int numberOfColumns;
} SummaryTable;


SummaryTable *SummaryTable_construct(int numberOfRows, int numberOfColumns);
void SummaryTable_destruct(SummaryTable *summaryTable);
void SummaryTable_increment(SummaryTable *summaryTable, int rowIndex, int columnIndex, double value);
char *SummaryTable_getRowString(SummaryTable *summaryTable, int rowIndex, char delimiter);
char *SummaryTable_getRowStringPercentage(SummaryTable *summaryTable, int rowIndex, char delimiter);


/*! @typedef
 * @abstract Structure for keeping a list of summary tables
 */
typedef struct SummaryTableList {
    stList *summaryTables;
    stList *categoryNames1;
    stList *categoryNames2;
    int numberOfRows;
    int numberOfColumns;
    int numberOfCategories1;
    int numberOfCategories2;
    int totalNumberOfTables;
} SummaryTableList;

int SummaryTableList_getTableIndex(SummaryTableList * summaryTableList, int catIndex1, int catIndex2);
SummaryTable *SummaryTableList_getTable(SummaryTableList *summaryTableList, int catIndex1, int catIndex2);
SummaryTableList *SummaryTableList_construct(stList *categoryNames1,
                                             stList *categoryNames2,
                                             int numberOfRows,
                                             int numberOfColumns);
void SummaryTableList_destruct(SummaryTableList *summaryTableList);

#endif //FLAGGER_SUMMARY_TABLE_H
