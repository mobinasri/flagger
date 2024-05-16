
#include "chunk.h"
#include "ptBlock.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>


bool testChunkIterator(char *covPath) {
    bool correct = true;
    // start (0-based), end (0-based), coverage, high_mapq, high_clip, region_index
    int truthValuesCtg1[6][6] = {{0,   19,  5,  2,  2,  0},
                                 {20,  39,  10, 10, 10, 1},
                                 {40,  59,  10, 10, 10, 1},
                                 {60,  79,  15, 15, 13, 1},
                                 {80,  99,  16, 16, 16, 1},
                                 {100, 109, 16, 16, 16, 1}};
    // annotation indices
    int truthValuesCtg2[1][6] = {{0, 9, 7, 7, 1, 0}};
    int truthAnnotationsCtg1[6][2] = {{1, -1},
                                      {1, -1},
                                      {1, -1},
                                      {1, 2},
                                      {2, -1},
                                      {2, -1}};
    int truthAnnotationsCtg2[1][2] = {{1, 2}};

    int truthAnnotationsLenCtg1[6] = {1,
                                      1,
                                      1,
                                      2,
                                      1,
                                      1};
    int truthAnnotationsLenCtg2[1] = {2};

    int truthLabelsCtg1[6][2] = {{1, 1},
                                 {2, 3},
                                 {2, 3},
                                 {3, 3},
                                 {3, 3},
                                 {3, 3}};

    int truthLabelsCtg2[1][2] = {{2, 2}};

    int windowLen = 20;
    int chunkCanonicalLen = 40;
    int nThreads = 2;
    int trackIndexCtg1 = -1;
    int trackIndexCtg2 = -1;
    ChunksCreator *chunksCreator = ChunksCreator_constructFromCov(covPath, NULL, chunkCanonicalLen, nThreads,
                                                                  windowLen);
    if (ChunksCreator_parseChunks(chunksCreator) != 0) {
        return false;
    }

    ChunkIterator *chunkIterator = ChunkIterator_construct(chunksCreator);

    char ctg[200];
    ptBlock *block = NULL;

    while ((block = ChunkIterator_getNextPtBlock(chunkIterator, ctg)) != NULL) {
        if (strcmp(ctg, "ctg1") == 0) {
            trackIndexCtg1 += 1;
            CoverageInfo *coverageInfo = block->data;
            int *truthValues = truthValuesCtg1[trackIndexCtg1];
            int s = block->rfs;
            int e = block->rfe;
            correct &= (s == truthValues[0]);
            correct &= (e == truthValues[1]);
            correct &= (coverageInfo->coverage == truthValues[2]);
            correct &= (coverageInfo->coverage_high_mapq == truthValues[3]);
            correct &= (coverageInfo->coverage_high_clip == truthValues[4]);
            correct &= (CoverageInfo_getRegionIndex(coverageInfo) == truthValues[5]);

            // check annotations
            int *truthAnnotations = truthAnnotationsCtg1[trackIndexCtg1];
            int truthAnnnotationsLen = truthAnnotationsLenCtg1[trackIndexCtg1];
            int len;
            int *annotationIndices = CoverageInfo_getAnnotationIndices(coverageInfo, &len);
            if (len != truthAnnnotationsLen) return false;
            for (int i = 0; i < len; i++) {
                correct &= (annotationIndices[i] == truthAnnotations[i]);
            }
            free(annotationIndices);

            // check labels
            int *labels = truthLabelsCtg1[trackIndexCtg1];
            if (coverageInfo->data == NULL) return false;
            Inference *inference = coverageInfo->data;
            correct &= (inference->truth == labels[0]);
            correct &= (inference->prediction == labels[1]);
        }
        if (strcmp(ctg, "ctg2") == 0) {
            trackIndexCtg2 += 1;
            CoverageInfo *coverageInfo = block->data;
            int *truthValues = truthValuesCtg2[trackIndexCtg2];
            int s = block->rfs;
            int e = block->rfe;
            correct &= (s == truthValues[0]);
            correct &= (e == truthValues[1]);
            correct &= (coverageInfo->coverage == truthValues[2]);
            correct &= (coverageInfo->coverage_high_mapq == truthValues[3]);
            correct &= (coverageInfo->coverage_high_clip == truthValues[4]);
            correct &= (CoverageInfo_getRegionIndex(coverageInfo) == truthValues[5]);


            // check annotations
            int *truthAnnotations = truthAnnotationsCtg2[trackIndexCtg2];
            int truthAnnnotationsLen = truthAnnotationsLenCtg2[trackIndexCtg2];
            int len;
            int *annotationIndices = CoverageInfo_getAnnotationIndices(coverageInfo, &len);
            if (len != truthAnnnotationsLen) return false;
            for (int i = 0; i < len; i++) {
                correct &= (annotationIndices[i] == truthAnnotations[i]);
            }
            free(annotationIndices);

            // check labels
            int *labels = truthLabelsCtg2[trackIndexCtg2];
            if (coverageInfo->data == NULL) return false;
            Inference *inference = coverageInfo->data;
            correct &= (inference->truth == labels[0]);
            correct &= (inference->prediction == labels[1]);
        }
    }
    ChunksCreator_destruct(chunksCreator);
    ChunkIterator_destruct(chunkIterator);
    return correct;
}




int main(int argc, char *argv[]) {

    bool allTestsPassed = true;

    // test 1
    bool test1Passed = testChunkIterator("tests/test_files/chunk_iterator/test_1_with_labels.cov");
    printf("[chunk_iterator] Test iterating over chunks with ChunkIterator:");
    printf(test1Passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");
    allTestsPassed &= test1Passed;


    if (allTestsPassed)
        return 0;
    else
        return 1;
}
