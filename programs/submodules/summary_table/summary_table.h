
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
    METRIC_BASE_LEVEL = 1
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
#define NUMBER_OF_METRIC_TYPES 2
#define NUMBER_OF_COMPARISON_TYPES 4

static const char* MetricTypeToString[2] = {"overlap_based", "base_level"};
static const char* CategoryTypeToString[2] = {"region", "annotation"};
static const char* ComparisonTypeToString[4] = {"recall", "precision", "truth", "prediction"};

/*! @typedef
 * @abstract Structure for keeping a summary table
 */
typedef struct SummaryTable {
    double **table;
    double **tablePercentage;
    double *totalPerRow;
    int numberOfRows;
    int numberOfColumns;
    stList *rowNames;
    stList *columnNames;
    char rowString[1000];
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

int SummaryTableList_getTableIndex(SummaryTableList *summaryTableList, int catIndex1, int catIndex2);

SummaryTable *SummaryTableList_getTable(SummaryTableList *summaryTableList, int catIndex1, int catIndex2);

SummaryTableList *SummaryTableList_construct(stList *categoryNames1,
                                             stList *categoryNames2,
                                             int numberOfRows,
                                             int numberOfColumns);

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
                                    char delimiter);

char *SummaryTableList_getRowStringPercentage(SummaryTableList *summaryTableList,
                                              int catIndex1,
                                              int catIndex2,
                                              int rowIndex,
                                              char delimiter);

void SummaryTableList_writeIntoFile(SummaryTableList *summaryTableList, FILE *fp, const char* linePrefix);
void SummaryTableList_writePercentageIntoFile(SummaryTableList *summaryTableList, FILE *fp, const char* linePrefix);

void SummaryTableList_destruct(SummaryTableList *summaryTableList);


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
    bool isMetricOverlapBased;
    // the overlap ratio threshold for considering an overlap as a hit in calculating overlap-based metrics
    double overlapThreshold;
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
                                                           bool isMetricOverlapBased,
                                                           double overlapThreshold);

SummaryTableUpdaterArgs *SummaryTableUpdaterArgs_copy(SummaryTableUpdaterArgs *src);

void SummaryTableUpdaterArgs_destruct(SummaryTableUpdaterArgs *args);

void convertBaseLevelToOverlapBased(int *refLabelConfusionRow,
                                    int columnSize,
                                    int refLabelBlockLength,
                                    double overlapThreshold);

void SummaryTableList_updateForAllCategory1(SummaryTableUpdaterArgs *argsTemplate, int sizeOfCategory1, int threads);

void SummaryTableList_updateByUpdaterArgsForThreadPool(void *argWork_);

void SummaryTableList_updateByUpdaterArgs(SummaryTableUpdaterArgs *args);

SummaryTableList *SummaryTableList_constructAndFillByIterator(void *blockIterator,
                                                              BlockIteratorType blockIteratorType,
                                                              stList *categoryNames,
                                                              CategoryType categoryType,
                                                              IntBinArray *sizeBinArray,
                                                              MetricType metricType,
                                                              double overlapRatioThreshold,
                                                              int numberOfLabels,
                                                              ComparisonType comparisonType,
                                                              int numberOfThreads);

#endif //FLAGGER_SUMMARY_TABLE_H
