#ifndef BIAS_DETECTOR_H
#define BIAS_DETECTOR_H


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
#include "count_data.h"
#include "hmm_utils.h"

/*! @typedef
 * @abstract Structure for finding annotations with coverage biases
 * @field countDataPerAnnotation                An array of CountData; one per annotation
 * @field numberOfAnnotations                   Total number of annotations
 * @field minCoverage                           Minimum coverage value (coverages lower than this will be ignored for bias detection)
 * @field minTotalCount                         An annotation to be detected as biased should have a total length higher than this
 * @field baselineAnnotationName                Name of the annotation whose most frequent coverage will used as a baseline for bias detection
 * @field minCovDiffNormalized                  The minimum value of the normalized coverage deviation from baseline to report an annotation as biased
 * @field annotationNames                       Annotation names
 */
typedef struct {
    CountData **countDataPerAnnotation;
    int numberOfAnnotations;
    char *baselineAnnotationName;
    int baselineAnnotationIndex;
    // thresholds
    double minCovDiffNormalized;
    int minCoverage;
    int minTotalCount;
    stList *annotationNames;
    // statistics
    bool statsAreUpdated;
    int *mostFrequentCoveragePerAnnotation;
    int *maxCountPerAnnotation;
    double *totalCountPerAnnotation;
    // result variables contating the mapping from annotations to region indices
    int *annotationToRegionMap;
    int *coveragePerRegion; // region 0 is reserved for baseline annotation
    int numberOfRegions; // region 0 is not biased but other regions are
    // for start_only mode
    bool startOnlyMode;
    int averageAlignmentLength;
    stHash *contigLengthTable;
    stList *chunks;
} BiasDetector;


BiasDetector *BiasDetector_construct(stList *annotationNames,
                                     const char *baselineAnnotationName,
                                     int minCoverage,
                                     int minTotalCount,
                                     double minCovDiffNormalized,
                                     int averageAlignmentLength,
                                     bool startOnlyMode,
                                     stHash *contigLengthTable);

void BiasDetector_destruct(BiasDetector *biasDetector);

void BiasDetector_setStatisticsPerAnnotation(BiasDetector *biasDetector, stHash *blockTable);

void BiasDetector_setMostFrequentCoveragePerAnnotation(BiasDetector *biasDetector);

void BiasDetector_setMaxCountPerAnnotation(BiasDetector *biasDetector);

void BiasDetector_setTotalCountPerAnnotation(BiasDetector *biasDetector);

void BiasDetector_runBiasDetection(BiasDetector *biasDetector,
                                   stList *annotationNamesToCheck,
                                   const char *tsvPathToWriteTable);

void BiasDetector_setCountDataPerAnnotationForStartOnlyMode(BiasDetector *biasDetector, stHash *blockTable);

void BiasDetector_setCountDataPerAnnotation(BiasDetector *biasDetector, stHash *blockTable);

#endif //BIAS_DETECTOR_H
