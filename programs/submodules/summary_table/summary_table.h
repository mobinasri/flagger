
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
#include "common.h"
#include "chunk.h"
#include "ptBlock.h"

typedef enum BlockIteratorType {
    ITERATOR_BY_CHUNK = 0,
    ITERATOR_BY_COV_BLOCK = 1
} BlockIteratorType;

typedef enum MetricType {
    METRIC_OVERLAP_BASED = 0,
    METRIC_BASE_LEVEL = 1,
    METRIC_AUN = 2
} MetricType;

typedef enum CategoryType {
    CATEGORY_REGION = 0,
    CATEGORY_ANNOTATION = 1
} CategoryType;

typedef enum ComparisonType {
    COMPARISON_TRUTH_VS_PREDICTION = 0,
    COMPARISON_PREDICTION_VS_TRUTH = 1,
    COMPARISON_TRUTH_VS_TRUTH = 2,
    COMPARISON_PREDICTION_VS_PREDICTION = 3,
} ComparisonType;

#define NUMBER_OF_CATEGORY_TYPES 2
#define NUMBER_OF_METRIC_TYPES 3
#define NUMBER_OF_COMPARISON_TYPES 4

static const char *MetricTypeToString[3] = {"overlap_based", "base_level", "truth_based_auN"};
static const char *CategoryTypeToString[2] = {"region", "annotation"};
static const char *ComparisonTypeToString[4] = {"TRUTH_VS_PREDICTION", "PREDICTION_VS_TRUTH", "TRUTH", "PREDICTION"};

/*! @typedef
 * @abstract Structure for keeping a summary table
 */
typedef struct SummaryTable {
    double **table;
    double **tablePercentage;
    double *totalPerRow;
    double *totalPerRowPercentage;
    double totalSum;
    int numberOfRows;
    int numberOfColumns;
    stList *rowNames;
    stList *columnNames;
    char rowString[10000];
    char rowName[1000];
} SummaryTable;


SummaryTable *SummaryTable_construct(int numberOfRows, int numberOfColumns);

SummaryTable *SummaryTable_constructByNames(stList *rowNames, stList *columnNames);

void SummaryTable_destruct(SummaryTable *summaryTable);

void SummaryTable_increment(SummaryTable *summaryTable, int rowIndex, int columnIndex, double value);

double SummaryTable_getValue(SummaryTable *summaryTable, int rowIndex, int columnIndex);

double SummaryTable_getTPCountInRow(SummaryTable *summaryTable, int rowIndex);

double SummaryTable_getAllCountInRow(SummaryTable *summaryTable, int rowIndex);

double SummaryTable_getAccuracyPercentage(SummaryTable *summaryTable, bool excludeLastRowAndColumn);

double SummaryTable_getMicroAverageAcrossRowsPercentage(SummaryTable *summaryTable, bool excludeLastRow);

double SummaryTable_getMacroAverageAcrossRowsPercentage(SummaryTable *summaryTable, bool excludeLastRow);

char *SummaryTable_getRowName(SummaryTable *summaryTable, int rowIndex);

char *SummaryTable_getRowString(SummaryTable *summaryTable, int rowIndex, char delimiter, bool addRowIndex);

char *SummaryTable_getRowStringPercentage(SummaryTable *summaryTable, int rowIndex, char delimiter, bool addRowIndex);

char *SummaryTable_getTotalPerRowString(SummaryTable *summaryTable, char delimiter, const char *prefixEntry);

char *SummaryTable_getTotalPerRowStringPercentage(SummaryTable *summaryTable, char delimiter, const char *prefixEntry);

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

int SummaryTableList_getTableIndex(SummaryTableList *summaryTableList, int catIndex1, int catIndex2);

SummaryTable *SummaryTableList_getTable(SummaryTableList *summaryTableList, int catIndex1, int catIndex2);

SummaryTableList *SummaryTableList_construct(stList *categoryNames1,
                                             stList *categoryNames2,
                                             int numberOfRows,
                                             int numberOfColumns);

SummaryTableList *SummaryTableList_constructByNames(stList *categoryNames1,
                                                    stList *categoryNames2,
                                                    stList *rowNames,
                                                    stList *columnNames);

void SummaryTableList_increment(SummaryTableList *summaryTableList,
                                int catIndex1,
                                int catIndex2,
                                int rowIndex,
                                int columnIndex,
                                double value);

double SummaryTableList_getValue(SummaryTableList *summaryTableList,
                                 int catIndex1,
                                 int catIndex2,
                                 int rowIndex,
                                 int columnIndex);

double SummaryTableList_getValuePercentage(SummaryTableList *summaryTableList,
                                           int catIndex1,
                                           int catIndex2,
                                           int rowIndex,
                                           int columnIndex);

char *SummaryTableList_getRowString(SummaryTableList *summaryTableList,
                                    int catIndex1,
                                    int catIndex2,
                                    int rowIndex,
                                    char delimiter,
                                    bool addRowIndex);

char *SummaryTableList_getRowStringPercentage(SummaryTableList *summaryTableList,
                                              int catIndex1,
                                              int catIndex2,
                                              int rowIndex,
                                              char delimiter,
                                              bool addRowIndex);

char *SummaryTableList_getTotalPerRowString(SummaryTableList *summaryTableList,
                                            int catIndex1,
                                            int catIndex2,
                                            char delimiter,
                                            const char *prefixEntry);

char *SummaryTableList_getTotalPerRowStringPercentage(SummaryTableList *summaryTableList,
                                                      int catIndex1,
                                                      int catIndex2,
                                                      char delimiter,
                                                      const char *prefixEntry);

void SummaryTableList_writeIntoFile(SummaryTableList *summaryTableList, FILE *fp, const char *linePrefix);

void SummaryTableList_writePercentageIntoFile(SummaryTableList *summaryTableList, FILE *fp, const char *linePrefix);

void SummaryTableList_writeTotalPerRowIntoFile(SummaryTableList *summaryTableList, FILE *fp, const char *linePrefix);

void SummaryTableList_writeTotalPerRowPercentageIntoFile(SummaryTableList *summaryTableList, FILE *fp,
                                                         const char *linePrefix);

void SummaryTableList_writeFinalStatisticsIntoFile(SummaryTableList *recallTableList,
                                                   SummaryTableList *precisionTableList,
                                                   FILE *fout,
                                                   const char *linePrefix);

void SummaryTableList_destruct(SummaryTableList *summaryTableList);


typedef struct SummaryTableListFullCatalog {
    SummaryTableList ****array;
    int dimMetricType;
    int dimComparisonType;
    int dimCategoryType;
} SummaryTableListFullCatalog;

SummaryTableListFullCatalog *
SummaryTableListFullCatalog_constructNull(int dimCategoryType, int dimMetricType, int dimComparisonType);

void SummaryTableListFullCatalog_update(SummaryTableListFullCatalog *catalog,
                                        SummaryTableList *summaryTableList,
                                        CategoryType categoryType,
                                        MetricType metricType,
                                        ComparisonType comparisonType);

SummaryTableList *SummaryTableListFullCatalog_get(SummaryTableListFullCatalog *catalog,
                                                  CategoryType categoryType,
                                                  MetricType metricType,
                                                  ComparisonType comparisonType);

void SummaryTableListFullCatalog_write(SummaryTableListFullCatalog *catalog,
                                       CoverageHeader *header,
                                       stList *labelNamesWithUnknown,
                                       char *outputPath);

void SummaryTableListFullCatalog_destruct(SummaryTableListFullCatalog *catalog);


typedef struct SummaryTableUpdaterArgs {
    // a ptBlock iterator
    // can be either ChunkIterator or ptBlockItrPerContig
    void *blockIterator;

    // a function for making a copy of the iterator
    void *(*copyBlockIterator)(void *);

    // a function for resetting the iterator to start from the first block
    void (*resetBlockIterator)(void *);

    // a function for freeing iterator memory
    void (*destructBlockIterator)(void *);

    // a function for fetching the next ptBlock from iterator
    // can be either ChunkIterator_getNextPtBlock or ptBlockItrPerContig_next
    ptBlock *(*getNextBlock)(void *, char *);

    // a struct containing a list of summary tables
    SummaryTableList *summaryTableList;
    // a struct containing the size bin intervals to stratify events based on their sizes
    IntBinArray *sizeBinArray;

    // a function to check if a coverage info has overlap with a category index
    // can be either CoverageInfo_overlapRegionIndex or CoverageInfo_overlapAnnotationIndex
    int (*overlapFuncCategoryIndex1)(CoverageInfo *, int);

    // index of category 1 to update its related sumary table. It can be either region index or annotation index
    int categoryIndex1;

    // a function for getting the reference label (rows in summary table) from the
    // Inference data embedded inside a CoverageInfo struct
    // (it can be either get_inference_prediction_label or get_inference_truth_label)
    // for calculating recall ref should be truth and for precision it should be prediction
    int8_t (*getRefLabelFunction)(Inference *);

    // a function for getting the query label (columns in summary table) from the
    // Inference data embedded inside a CoverageInfo struct
    // (it can be either get_inference_prediction_label or get_inference_truth_label)
    // for calculating recall query should be prediction and for precision it should be truth
    int8_t (*getQueryLabelFunction)(Inference *);

    // There are two types of metric for recall and precision
    // The default metric is base-level which is calculated by just counting the number of bases
    // fall into each entry of the confusion matrix (represented in summary table)
    // For example for entry[0][3] if isMetricOverlapBases=false it counts the number of bases/windows
    // that has ref label of 0 and prediction label of 3
    // However if this attribute is true (overlap-based metric) it counts the number of hits based on overlap threshold.
    // For example if there is a contiguous block with a ref label of 0 and its length is 20 then if at least
    // (overlapThreshold * 20) bases/windows in this contiguous block has the query label of 3
    // then entry[0][3] will be incremented by one.
    MetricType metricType;
    // the overlap ratio threshold for considering an overlap as a hit in calculating overlap-based metrics
    double overlapThreshold;

    // auxiliary summary table list
    // it can be useful when for creating one summary table we need to get some information from another table
    // for example for METRIC_AUN, which is defined only for TRUTH_VS_PREDICTION or PREDICTION_VS_TRUTH,
    // we need the summary table of METRIC_BASE_LEVEL for TRUTH_VS_TRUTH and PREDICTION_VS_PREDICTION respectively
    // to compute what ratio of whole bases with a specific truth label (or prediction label) is covered by
    // each truth label (prediction label).
    SummaryTableList *auxiliarySummaryTableList;
} SummaryTableUpdaterArgs;

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
                                                           SummaryTableList *auxiliarySummaryTableList);

SummaryTableUpdaterArgs *SummaryTableUpdaterArgs_copy(SummaryTableUpdaterArgs *src);

void SummaryTableUpdaterArgs_destruct(SummaryTableUpdaterArgs *args);

void convertBaseLevelToOverlapBased(double *refLabelConfusionRow,
                                    int columnSize,
                                    int refLabelBlockLength,
                                    double overlapThreshold);

void SummaryTableList_addCreationJobsForAllCategory1(SummaryTableUpdaterArgs *argsTemplate, int sizeOfCategory1, tpool_t *threadPool);

void SummaryTableList_updateByUpdaterArgsForThreadPool(void *argWork_);

void SummaryTableList_updateByUpdaterArgs(SummaryTableUpdaterArgs *args);


void SummaryTableList_addCreationJobsForOneMetricType(CoverageHeader *header,
                                                      void *blockIterator,
                                                      BlockIteratorType blockIteratorType,
                                                      IntBinArray *sizeBinArray,
                                                      MetricType metricType,
                                                      double overlapRatioThreshold,
                                                      int numberOfLabelsWithUnknown,
                                                      stList *labelNamesWithUnknown,
                                                      SummaryTableListFullCatalog *catalog,
                                                      tpool_t *threadPool);


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
                                                 tpool_t *threadPool);

void SummaryTableList_createAndWriteAllTables(void *iterator,
                                              BlockIteratorType blockIteratorType,
                                              CoverageHeader *header,
                                              const char *outputPath,
                                              const char *binArrayFilePath,
                                              stList *labelNamesWithUnknown,
                                              double overlapRatioThreshold,
                                              int threads);

#endif //FLAGGER_SUMMARY_TABLE_H
