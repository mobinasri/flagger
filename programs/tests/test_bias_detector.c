#include "bias_detector.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "ptBlock.h"
#include "sonLib.h"

bool test_BiasDetector_getAnnotationToRegionMap(char *bamPath, char *jsonPath) {
    int threads = 2;
    int minMapq = 10;
    int minClip = 0.1;
    double downsampleRate = 1.0;
    int minAlignmentLength = 0;
    bool startOnlyMode = false;
    int averageAlignmentLength = 0;
    stHash *blockTable = ptBlock_multi_threaded_coverage_extraction_with_zero_coverage_and_annotation(bamPath, NULL,
                                                                                                      downsampleRate, jsonPath,
                                                                                                      threads, minMapq,
                                                                                                      minClip,
                                                                                                      minAlignmentLength,
                                                                                                      startOnlyMode,
                                                                                                      &averageAlignmentLength);

    const char *annotationZeroName = "no_annotation";
    stList *annotationNames = parse_annotation_names_and_save_in_stList(jsonPath, annotationZeroName);
    if (stList_length(annotationNames) != 7) return false;

    const char *baselineAnnotationName = "annot_1";
    int minCoverage = 1;
    int minTotalCount = 1;
    double minCovDiffNormalized = 0.2;
    BiasDetector *biasDetector = BiasDetector_construct(annotationNames, 
		                                        baselineAnnotationName, 
							minCoverage,
                                                        minTotalCount, 
							minCovDiffNormalized,
							averageAlignmentLength,
							startOnlyMode,
							NULL);

    BiasDetector_setStatisticsPerAnnotation(biasDetector, blockTable);
    /*for(int i=0; i < 7; i++){
        fprintf(stderr, "%s %d mod=%d tot=%d max=%d\n", stList_get(annotationNames,i), i, biasDetector->mostFrequentCoveragePerAnnotation[i], biasDetector->totalCountPerAnnotation[i], biasDetector->maxCountPerAnnotation[i]);
    }*/

    stList *annotationNamesToCheck = stList_construct3(0, free);
    stList_append(annotationNamesToCheck, copyString("annot_1"));
    stList_append(annotationNamesToCheck, copyString("annot_2"));
    stList_append(annotationNamesToCheck, copyString("annot_3"));
    stList_append(annotationNamesToCheck, copyString("annot_4"));
    stList_append(annotationNamesToCheck, copyString("annot_5"));

    const char *tsvPathToWriteTable = "tests/test_files/test_bias_detector.tsv";
    BiasDetector_runBiasDetection(biasDetector, annotationNamesToCheck, tsvPathToWriteTable);
    int *annotationToRegionMap = biasDetector->annotationToRegionMap;
    int *coveragePerRegion = biasDetector->coveragePerRegion;

    bool correct = true;
    if (biasDetector->numberOfRegions != 2) return false;

    // check mapping between annotation index and region index
    correct &= (annotationToRegionMap[0] == 0); // no_annotation -> 0
    correct &= (annotationToRegionMap[1] == 0); // annot_1 -> 0
    correct &= (annotationToRegionMap[2] == 0); // annot_2 -> 0
    correct &= (annotationToRegionMap[3] == 0); // annot_3 -> 0
    correct &= (annotationToRegionMap[4] == 1); // annot_4 -> 1
    correct &= (annotationToRegionMap[5] == 0); // annot_5 -> 0
    correct &= (annotationToRegionMap[6] == 0); // annot_5 -> 0


    // check median coverage per region
    correct &= (coveragePerRegion[0] == 1);
    correct &= (coveragePerRegion[1] == 2);

    stHash_destruct(blockTable);
    stList_destruct(annotationNames);
    stList_destruct(annotationNamesToCheck);
    BiasDetector_destruct(biasDetector);
    return correct;
}


int main(int argc, char *argv[]) {
    char bamPath[1000] = "tests/test_files/bias_detector/test.bam";
    char jsonPath[1000] = "tests/test_files/bias_detector/test.json";

    bool allTestsPassed = true;

    // test 1
    bool test1Passed = test_BiasDetector_getAnnotationToRegionMap(bamPath, jsonPath);
    printf("[bias_detector] Test bias detection and getting annotation to region map:");
    printf(test1Passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");
    allTestsPassed &= test1Passed;

    if (allTestsPassed)
        return 0;
    else
        return 1;

}
