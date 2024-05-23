#include <getopt.h>
#include <time.h>
#include "bgzf.h"
#include "sonLib.h"
#include <assert.h>
#include <math.h>
#include <float.h>
#include "common.h"
#include <time.h>
#include <string.h>
#include "ptBlock.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include "stdlib.h"
#include "sonLib.h"
#include "chunk.h"
#include "track_reader.h"
#include "summary_table.h"


static struct option long_options[] =
        {
                {"input",                 required_argument, NULL, 'i'},
                {"binArrayFile",          required_argument, NULL, 'b'},
                {"output",                required_argument, NULL, 'o'},
                {"threads",               required_argument, NULL, 't'},
                {"overlapRatioThreshold", required_argument, NULL, 'v'},
                {"labelNames",            required_argument, NULL, 'l'},
                {NULL,                    0,                 NULL, 0}
        };


int main(int argc, char *argv[]) {
    int c;
    double overlapRatioThreshold = 0.4;
    char *inputPath = NULL;
    char *outputPath = NULL;
    char *binArrayFilePath = NULL;
    stList *labelNamesWithUnknown = NULL;
    int threads = 4;
    char *program;
    (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
    while (~(c = getopt_long(argc, argv, "i:o:b:t:v:l:h", long_options, NULL))) {
        switch (c) {
            case 'i':
                inputPath = optarg;
                break;
            case 'o':
                outputPath = optarg;
                break;
            case 'b':
                binArrayFilePath = optarg;
                break;
            case 't':
                threads = atoi(optarg);
                break;
            case 'v':
                overlapRatioThreshold = atof(optarg);
                break;
            case 'l':
                labelNamesWithUnknown = Splitter_getStringList(optarg, ',');
                stList_append(labelNamesWithUnknown, copyString("Unk"));
                break;
            default:
                if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
            help:
                fprintf(stderr, "\nUsage: %s  -i <INPUT_FILE> -b <BIN_ARRAY_TSV> -o <OUTPUT_FILE> \n", program);
                fprintf(stderr, "Options:\n");
                fprintf(stderr,
                        "         -i,--input         input path (can have formats '.cov', '.cov.gz', '.bed', '.bed.gz' or '.bin')\n");
                fprintf(stderr,
                        "         -b,--binArrayFile         a tsv file (tab-delimited) that contains bin arrays for stratifying results by event size. "
                        "It should contain three columns. 1st column is the closed start of the bin and the 2nd column is the open end."
                        "The 3rd column has a name for each bin. For example one row can be '0\t100\t[0-100). If no file is passed it will"
                        "consider one large bin as the default value. (Default = [0,1e9) with the name 'ALL')\n");
                fprintf(stderr,
                        "         -o,--output         output path for saving summary table (should have format '.tsv')\n");
                fprintf(stderr,
                        "         -v, --overlapRatioThreshold                        Minimum overlap ratio in calculating overlap-based metrics for considering a hit between a ref label "
                        "(for example truth label for recall) and query label (for example prediction label for recall) [default: 0.4]\n");
                fprintf(stderr,
                        "         -t, --threads                        Number of threads for parallelizing creating table for each annotation/region [default: 4]\n");
                fprintf(stderr,
                        "         -l, --labelNames                    (Optional) A comma-delimited string of label names (for example 'Err,Dup,Hap,Col'). It should match the number of labels in the header of the input file.[default: none]\n");

                return 1;
        }
    }

    double realtimeStart = System_getRealTimePoint();

    if (inputPath == NULL || outputPath == NULL) {
        fprintf(stderr, "[%s] Error: Input/Output path cannot be NULL.\n", get_timestamp());
        exit(EXIT_FAILURE);
    }

    // check input/output file extensions
    char *inputExtension = extractFileExtension(inputPath);
    char *outputExtension = extractFileExtension(outputPath);
    if (strcmp(outputExtension, "tsv") != 0) {
        fprintf(stderr, "[%s] Error: output file should have tsv extension  !\n", get_timestamp());
        free(inputExtension);
        free(outputExtension);
        exit(EXIT_FAILURE);
    }
    if (strcmp(inputExtension, "cov") != 0 &&
        strcmp(inputExtension, "cov.gz") != 0 &&
        strcmp(inputExtension, "bed") != 0 &&
        strcmp(inputExtension, "bed.gz") != 0 &&
        strcmp(inputExtension, "bin") != 0) {
        fprintf(stderr,
                "[%s] Error: input file should either cov/cov.gz/bed/bed.gz or a binary file made with create_bin_chunks.\n",
                get_timestamp());
        free(inputExtension);
        free(outputExtension);
        exit(EXIT_FAILURE);
    }

    // define object and functions for iterating blocks
    void *iterator;
    CoverageHeader *header = NULL;
    ChunksCreator *chunksCreator = NULL;
    stHash *blocksPerContig = NULL;
    BlockIteratorType blockIteratorType =
            strcmp(inputExtension, "bin") == 0 ? ITERATOR_BY_CHUNK : ITERATOR_BY_COV_BLOCK;

    // create block iterator
    if (blockIteratorType == ITERATOR_BY_CHUNK) { // for binary file
        chunksCreator = ChunksCreator_constructEmpty();
        ChunkCreator_parseChunksFromBinaryFile(chunksCreator, inputPath);
        iterator = (void *) ChunkIterator_construct(chunksCreator);
        header = chunksCreator->header;
    } else if (blockIteratorType == ITERATOR_BY_COV_BLOCK) { // for cov or bed file
        blocksPerContig = ptBlock_parse_coverage_info_blocks(inputPath);
        iterator = (void *) ptBlockItrPerContig_construct(blocksPerContig);
        header = CoverageHeader_construct(inputPath);
    }

    // parse bin intervals
    IntBinArray *binArray = IntBinArray_constructFromFile(binArrayFilePath);

    // open file for writing summary tables
    FILE *fout = fopen(outputPath, "w");
    if (fout == NULL) {
        fprintf(stderr, "[%s] Error: %s cannot be opened.\n", get_timestamp(), outputPath);
        exit(EXIT_FAILURE);
    }

    char outputPathFinalStats[1000];
    char prefix[1000];
    memcpy(prefix, outputPath, strlen(outputPath) - 4); // ".tsv" has 4 characters
    prefix[strlen(outputPath) - 4] = '\0'; // ".tsv" has 4 characters
    sprintf(outputPathFinalStats, "%s.final_stats.tsv", prefix);

    FILE *foutFinalStats = NULL;
    if (header->isTruthAvailable && header->isPredictionAvailable) {
        foutFinalStats = fopen(outputPathFinalStats, "w");
        if (foutFinalStats == NULL) {
            fprintf(stderr, "[%s] Error: %s cannot be opened for making a tsv file with final benchmarking stats.\n",
                    get_timestamp(), outputPathFinalStats);
            exit(EXIT_FAILURE);
        }
        fprintf(foutFinalStats,
                "#Metric_Type\tCategory_Type\tCategory_Name\tSize_Bin_Name\tLabel\tTP_Prediction_Ref\tTP_Truth_Ref\tFP\tFN\tTotal_Prediction_Ref\tTotal_Truth_Ref\tPrecision\tRecall\tF1-Score\tAccuracy_Prediction_Ref\tAccuracy_Truth_Ref\n");
    }

    // The list of label names has an additional name "Unk" for handling labels with the value of -1
    if (labelNamesWithUnknown != NULL && stList_length(labelNamesWithUnknown) - 1 != header->numberOfLabels) {
        fprintf(stderr, "[%s] Error: Number of label names %d  does not match the number of labels in the header %d.\n",
                get_timestamp(),
                stList_length(labelNamesWithUnknown) - 1,
                header->numberOfLabels);
        exit(EXIT_FAILURE);
    }


    int numberOfLabelsWithUnknown = header->numberOfLabels + 1;

    // write column names in the header line
    char linePrefix[1000];
    sprintf(linePrefix, "#Statistic\tMetric_Type\tEntry_Type\tCategory_Type\tCategory_Name\tSize_Bin_Name\tRef_Label");
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


    SummaryTableList **summaryTableListPerComparison = (SummaryTableList **) malloc(
            NUMBER_OF_COMPARISON_TYPES * sizeof(SummaryTableList *));
    for (int comparisonType = 0; comparisonType < NUMBER_OF_COMPARISON_TYPES; comparisonType++) {
        summaryTableListPerComparison[comparisonType] = NULL;
    }

    if (header->isTruthAvailable || header->isPredictionAvailable) {
        // iterating over category types; annotation and region
        for (int categoryType = 0; categoryType < NUMBER_OF_CATEGORY_TYPES; categoryType++) {
            // iterating over metric types; base-level and overlap-based
            for (int metricType = 0; metricType < NUMBER_OF_METRIC_TYPES; metricType++) {
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
                    SummaryTableList *summaryTableList =
                            SummaryTableList_constructAndFillByIterator(iterator,
                                                                        blockIteratorType,
                                                                        categoryNames,
                                                                        categoryType,
                                                                        binArray,
                                                                        metricType,
                                                                        overlapRatioThreshold,
                                                                        numberOfLabelsWithUnknown,
                                                                        labelNamesWithUnknown,
                                                                        comparisonType,
                                                                        threads);
                    summaryTableListPerComparison[comparisonType] = summaryTableList;

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
                } // end comparison type
                // if both labels are present for printing benchmarking stats
                if (header->isTruthAvailable && header->isPredictionAvailable) {
                    SummaryTableList *precisionTables = summaryTableListPerComparison[COMPARISON_PREDICTION_VS_TRUTH];
                    SummaryTableList *recallTables = summaryTableListPerComparison[COMPARISON_TRUTH_VS_PREDICTION];

                    // write percentages
                    sprintf(linePrefix, "%s\t%s",
                            MetricTypeToString[metricType],
                            CategoryTypeToString[categoryType]);
                    SummaryTableList_writeFinalStatisticsIntoFile(recallTables,
                                                                  precisionTables,
                                                                  foutFinalStats,
                                                                  linePrefix);
                }
            }// end metric type
        }// end category type
    }


    // free summary tables
    for (int comparisonType = 0; comparisonType < NUMBER_OF_COMPARISON_TYPES; comparisonType++) {
        if (summaryTableListPerComparison[comparisonType] != NULL) {
            SummaryTableList_destruct(summaryTableListPerComparison[comparisonType]);
        }
    }
    // close file
    fclose(fout);

    // free memory
    free(outputExtension);
    free(inputExtension);
    IntBinArray_destruct(binArray);

    if (blockIteratorType == ITERATOR_BY_CHUNK) {
        // free iterator
        ChunkIterator_destruct((ChunkIterator *) iterator);
        // free chunks
        ChunksCreator_destruct(chunksCreator);
    } else {
        // free iterator
        ptBlockItrPerContig_destruct((ptBlockItrPerContig *) iterator);
        // free blocks
        stHash_destruct(blocksPerContig);
        // free header
        CoverageHeader_destruct(header);
    }
    fprintf(stderr, "[%s] Done!\n", get_timestamp());

    // log used time/resources
    double realtime = System_getRealTimePoint() - realtimeStart;
    double cputime = System_getCpuTime();
    double rssgb = System_getPeakRSSInGB();
    double usage = System_getCpuUsage(cputime, realtime);
    // copied from https://github.com/chhylp123/hifiasm/blob/70fd9a0b1fea45e442eb5f331922ea91ef4f71ae/main.cpp#L73
    fprintf(stderr, "Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB; CPU usage: %.1f\%\n", realtime, cputime,
            rssgb, usage * 100.0);
}
