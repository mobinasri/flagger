#include <getopt.h>
#include "sam.h"
#include "faidx.h"
#include <time.h>
#include "bgzf.h"
#include <regex.h>
#include "sonLib.h"
#include <assert.h>
#include <math.h>
#include <float.h>
#include "cigar_it.h"
#include "common.h"
#include "vcf.h"
#include "ptBlock.h"
#include "ptAlignment.h"
#include <time.h>
#include <string.h>
#include "ptBlock.h"
#include "ptAlignment.h"
#include "tpool.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include "cJSON.h"
#include "stdlib.h"
#include "bias_detector.h"
#include "hmm_utils.h"

BiasDetector *BiasDetector_construct(stList *annotationNames,
                                     const char *baselineAnnotationName,
                                     int minCoverage,
                                     int minTotalCount,
                                     double minCovDiffNormalized) {
    BiasDetector *biasDetector = malloc(sizeof(BiasDetector));
    biasDetector->numberOfAnnotations = stList_length(annotationNames);
    biasDetector->countDataPerAnnotation = CountData_construct1DArray(MAX_COVERAGE_VALUE,
                                                                      biasDetector->numberOfAnnotations);
    biasDetector->minCoverage = minCoverage;
    biasDetector->minTotalCount = minTotalCount;
    biasDetector->baselineAnnotationName = copyString(baselineAnnotationName);
    biasDetector->baselineAnnotationIndex = get_annotation_index(annotationNames, baselineAnnotationName);
    biasDetector->minCovDiffNormalized = minCovDiffNormalized;

    biasDetector->annotationNames = stList_construct3(0, free);
    for (int i = 0; i < stList_length(annotationNames); i++) {
        stList_append(biasDetector->annotationNames, copyString(stList_get(annotationNames, i)));
    }
    biasDetector->mostFrequentCoveragePerAnnotation = Int_construct1DArray(biasDetector->numberOfAnnotations);
    biasDetector->maxCountPerAnnotation = Int_construct1DArray(biasDetector->numberOfAnnotations);
    biasDetector->totalCountPerAnnotation = Int_construct1DArray(biasDetector->numberOfAnnotations);
    biasDetector->annotationToRegionMap = Int_construct1DArray(biasDetector->numberOfAnnotations);
    biasDetector->numberOfRegions = 0;
    biasDetector->coveragePerRegion = NULL; // we don't know the number of regions yet
    biasDetector->statsAreUpdated = false;
    return biasDetector;
}

void BiasDetector_destruct(BiasDetector *biasDetector) {
    CountData_destruct1DArray(biasDetector->countDataPerAnnotation, biasDetector->numberOfAnnotations);
    stList_destruct(biasDetector->annotationNames);
    free(biasDetector->mostFrequentCoveragePerAnnotation);
    free(biasDetector->maxCountPerAnnotation);
    free(biasDetector->totalCountPerAnnotation);
    free(biasDetector->baselineAnnotationName);
    free(biasDetector->coveragePerRegion);
    free(biasDetector->annotationToRegionMap);
    free(biasDetector);
}


void BiasDetector_setStatisticsPerAnnotation(BiasDetector *biasDetector, stHash *blockTable) {
    BiasDetector_setCountDataPerAnnotation(biasDetector, blockTable);
    BiasDetector_setMostFrequentCoveragePerAnnotation(biasDetector);
    BiasDetector_setMaxCountPerAnnotation(biasDetector);
    BiasDetector_setTotalCountPerAnnotation(biasDetector);
    biasDetector->statsAreUpdated = true;
}

void BiasDetector_setCountDataPerAnnotation(BiasDetector *biasDetector, stHash *blockTable) {
    stHashIterator *it = stHash_getIterator(blockTable);
    char *contigName;
    while ((contigName = stHash_getNext(it)) != NULL) {
        stList *blocks = stHash_search(blockTable, contigName);
        for (int i = 0; i < stList_length(blocks); i++) {
            ptBlock *block = stList_get(blocks, i);
            CoverageInfo *covInfo = (CoverageInfo *) block->data;
            for (int annotationIndex = 0; annotationIndex < biasDetector->numberOfAnnotations; annotationIndex++) {
                if (CoverageInfo_overlapAnnotationIndex(covInfo,
                                                        annotationIndex)) { // checking overlap is fast (a few bit-wise operations)
                    int count = block->rfe - block->rfs + 1;
                    int value = covInfo->coverage;
                    CountData_increment(biasDetector->countDataPerAnnotation[annotationIndex], value, (double) count);
                }
            }
        }
    }
}

void BiasDetector_setMostFrequentCoveragePerAnnotation(BiasDetector *biasDetector) {
    for (int annotationIndex = 0; annotationIndex < biasDetector->numberOfAnnotations; annotationIndex++) {
        int mostFreqCoverage = CountData_getMostFrequentValue(biasDetector->countDataPerAnnotation[annotationIndex],
                                                              biasDetector->minCoverage, MAX_COVERAGE_VALUE);
        biasDetector->mostFrequentCoveragePerAnnotation[annotationIndex] = mostFreqCoverage;
    }
}

void BiasDetector_setMaxCountPerAnnotation(BiasDetector *biasDetector) {
    for (int annotationIndex = 0; annotationIndex < biasDetector->numberOfAnnotations; annotationIndex++) {
        int maxCount = (int) CountData_getMaxCount(biasDetector->countDataPerAnnotation[annotationIndex],
                                                   biasDetector->minCoverage, MAX_COVERAGE_VALUE);
        biasDetector->maxCountPerAnnotation[annotationIndex] = maxCount;
    }
}

void BiasDetector_setTotalCountPerAnnotation(BiasDetector *biasDetector) {
    for (int annotationIndex = 0; annotationIndex < biasDetector->numberOfAnnotations; annotationIndex++) {
        int totalCount = (int) CountData_getTotalCount(biasDetector->countDataPerAnnotation[annotationIndex], 0,
                                                       MAX_COVERAGE_VALUE);
        biasDetector->totalCountPerAnnotation[annotationIndex] = totalCount;
    }
}


void BiasDetector_runBiasDetection(BiasDetector *biasDetector,
                                   stList *annotationNamesToCheck,
                                   const char *tsvPathToWriteTable) {

    if (biasDetector->statsAreUpdated == false) {
        fprintf(stderr,
                "[%s] Warning: Bias detection cannot be executed. Run the function BiasDetector_setStatisticsPerAnnotation() for updating stats.\n",
                get_timestamp());
        return;
    }

    int numberOfAnnotations = biasDetector->numberOfAnnotations;
    int *annotationToRegionMap = biasDetector->annotationToRegionMap;

    int baselineIndex = biasDetector->baselineAnnotationIndex;
    int baselineCoverage = biasDetector->mostFrequentCoveragePerAnnotation[baselineIndex];
    if (baselineCoverage < 1) {
        fprintf(stderr, "[%s] Error: baseline coverage is too low (< 1) %d!\n", get_timestamp(), baselineCoverage);
        exit(EXIT_FAILURE);
    }

    if (annotationNamesToCheck == NULL || stList_length(annotationNamesToCheck) == 0) {
        fprintf(stderr,
                "[%s] Warning: annotationNamesToCheck is empty for bias detection. All annotations are mapped to region index 0!\n",
                get_timestamp());
        biasDetector->numberOfRegions = 1;
        biasDetector->coveragePerRegion = malloc(1 * sizeof(int));
        biasDetector->coveragePerRegion[0] = baselineCoverage;
    }

    FILE *tsvFile = NULL;
    if (tsvPathToWriteTable == NULL) {
        fprintf(stderr,
                "[%s] Warning: Summary table for coverage bias detection will not be written to a file since no path is provided!\n",
                get_timestamp());
    } else {
        tsvFile = fopen(tsvPathToWriteTable, "w");
        if (tsvFile == NULL) {
            fprintf(stderr, "[%s] Warning: %s cannot be opened for writing coverage bias detection table!\n",
                    get_timestamp(), tsvPathToWriteTable);
        }
    }

    if (tsvFile != NULL) {
        fprintf(stderr, "[%s] Writing coverage bias table into this file %s \n", get_timestamp(), tsvPathToWriteTable);
        fprintf(tsvFile,
                "annotation\tstatus\tmost_freq_cov\tmax_count\ttotal_count\tcov_diff_normalized\tregion_index\n");
    }

    int regionIndex = 1; // starts from 1 since 0 is reserved for regions with no bias
    biasDetector->coveragePerRegion = realloc(biasDetector->coveragePerRegion, regionIndex * sizeof(int));
    biasDetector->coveragePerRegion[0] = baselineCoverage;

    for (int annotIndex = 0; annotIndex < numberOfAnnotations; annotIndex++) {
        char *annotationName = stList_get(biasDetector->annotationNames, annotIndex);
        // if this annotation is not among the ones to check go to the next one
        if (get_annotation_index(annotationNamesToCheck, annotationName) == -1) continue;

        // get stats for detecting coverage bias
        int mostFrequentCoverage = biasDetector->mostFrequentCoveragePerAnnotation[annotIndex];
        int maxCount = biasDetector->maxCountPerAnnotation[annotIndex];
        int totalCount = biasDetector->totalCountPerAnnotation[annotIndex];

        double covDiffNormalized = ((double) mostFrequentCoverage - baselineCoverage) / baselineCoverage;
        if (mostFrequentCoverage < 1) {
            covDiffNormalized = 0.0;
        }
        double covDiffNormalizedAbs = covDiffNormalized < 0 ? -1 * covDiffNormalized : covDiffNormalized;
        bool isBiased = (biasDetector->minCovDiffNormalized < covDiffNormalizedAbs) &&
                        (biasDetector->minTotalCount < totalCount) && mostFrequentCoverage != -1;
        // if annotation is biased add it as a new region index
        if (isBiased == true) {
            annotationToRegionMap[annotIndex] = regionIndex;
            biasDetector->coveragePerRegion = realloc(biasDetector->coveragePerRegion, (regionIndex + 1) * sizeof(int));
            biasDetector->coveragePerRegion[regionIndex] = mostFrequentCoverage;
            regionIndex++;
        }
        if (tsvFile != NULL) {
            fprintf(tsvFile, "%s\t%s\t%d\t%d\t%d\t%+.3f\t%d\n",
                    annotationName,
                    isBiased ? "biased" : "not_biased",
                    mostFrequentCoverage,
                    maxCount,
                    totalCount,
                    covDiffNormalized,
                    annotationToRegionMap[annotIndex]);
        }

    }

    if (tsvFile != NULL) {
        fclose(tsvFile);
    }
    biasDetector->numberOfRegions = regionIndex;
}


