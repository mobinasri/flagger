#include "track_reader.h"
#include "chunk.h"
#include "sonLib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

bool testCreatingChunks(char *covPath) {
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

    int truthAnnotationsLenCtg1[6] = {1, 1, 1, 2, 1, 1};
    int truthAnnotationsLenCtg2[1] = {2};

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
    stList *chunks = chunksCreator->chunks;
    int numberOfChunks = stList_length(chunks);
    for (int chunkIndex = 0; chunkIndex < numberOfChunks; chunkIndex++) {
        Chunk *chunk = stList_get(chunks, chunkIndex);
        if (strcmp(chunk->ctg, "ctg1") == 0) {
            for (int windowIndex = 0; windowIndex < chunk->coverageInfoSeqLen; windowIndex++) {
                trackIndexCtg1 += 1;
                CoverageInfo *coverageInfo = chunk->coverageInfoSeq[windowIndex];
                int *truthValues = truthValuesCtg1[trackIndexCtg1];
                int *truthAnnotations = truthAnnotationsCtg1[trackIndexCtg1];
                int truthAnnnotationsLen = truthAnnotationsLenCtg1[trackIndexCtg1];
                int s = chunk->s + windowLen * windowIndex;
                int e = min(chunk->e, s + windowLen - 1);
                correct &= (s == truthValues[0]);
                correct &= (e == truthValues[1]);
                correct &= (coverageInfo->coverage == truthValues[2]);
                correct &= (coverageInfo->coverage_high_mapq == truthValues[3]);
                correct &= (coverageInfo->coverage_high_clip == truthValues[4]);
                correct &= (CoverageInfo_getRegionIndex(coverageInfo) == truthValues[5]);
                int len;
                int *annotationIndices = CoverageInfo_getAnnotationIndices(coverageInfo, &len);
                if (len != truthAnnnotationsLen) return false;
                for (int i = 0; i < len; i++) {
                    correct &= (annotationIndices[i] == truthAnnotations[i]);
                }
                free(annotationIndices);
            }
        }
        if (strcmp(chunk->ctg, "ctg2") == 0) {
            for (int windowIndex = 0; windowIndex < chunk->coverageInfoSeqLen; windowIndex++) {
                trackIndexCtg2 += 1;
                CoverageInfo *coverageInfo = chunk->coverageInfoSeq[windowIndex];
                int *truthValues = truthValuesCtg2[trackIndexCtg2];
                int *truthAnnotations = truthAnnotationsCtg2[trackIndexCtg2];
                int truthAnnnotationsLen = truthAnnotationsLenCtg2[trackIndexCtg2];
                int s = chunk->s + windowLen * windowIndex;
                int e = min(chunk->e, s + windowLen - 1);
                correct &= (s == truthValues[0]);
                correct &= (e == truthValues[1]);
                correct &= (coverageInfo->coverage == truthValues[2]);
                correct &= (coverageInfo->coverage_high_mapq == truthValues[3]);
                correct &= (coverageInfo->coverage_high_clip == truthValues[4]);
                correct &= (CoverageInfo_getRegionIndex(coverageInfo) == truthValues[5]);
                int len;
                int *annotationIndices = CoverageInfo_getAnnotationIndices(coverageInfo, &len);
                if (len != truthAnnnotationsLen) return false;
                for (int i = 0; i < len; i++) {
                    correct &= (annotationIndices[i] == truthAnnotations[i]);
                }
                free(annotationIndices);
            }
        }
    }
    ChunksCreator_destruct(chunksCreator);
    return correct;
}


bool testCreatingChunksWithLabels(char *covPath, bool doWriteAndRead) {
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

    if(doWriteAndRead){
        ChunksCreator_writeChunksIntoBinaryFile(chunksCreator, "tests/test_files/chunks_creator/tmp.bin");

        // free previous ChunksCreator and make an empty for parsing the saved binary file
        ChunksCreator_destruct(chunksCreator);
        chunksCreator = ChunksCreator_constructEmpty();
        ChunksCreator_parseChunksFromBinaryFile(chunksCreator, "tests/test_files/chunks_creator/tmp.bin");
    }


    stList *chunks = chunksCreator->chunks;
    int numberOfChunks = stList_length(chunks);
    for (int chunkIndex = 0; chunkIndex < numberOfChunks; chunkIndex++) {
        Chunk *chunk = stList_get(chunks, chunkIndex);
        if (strcmp(chunk->ctg, "ctg1") == 0) {
            for (int windowIndex = 0; windowIndex < chunk->coverageInfoSeqLen; windowIndex++) {
                trackIndexCtg1 += 1;
                CoverageInfo *coverageInfo = chunk->coverageInfoSeq[windowIndex];
                int *truthValues = truthValuesCtg1[trackIndexCtg1];
                int s = chunk->s + windowLen * windowIndex;
                int e = min(chunk->e, s + windowLen - 1);
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
        }
        if (strcmp(chunk->ctg, "ctg2") == 0) {
            for (int windowIndex = 0; windowIndex < chunk->coverageInfoSeqLen; windowIndex++) {
                trackIndexCtg2 += 1;
                CoverageInfo *coverageInfo = chunk->coverageInfoSeq[windowIndex];
                int *truthValues = truthValuesCtg2[trackIndexCtg2];
                int s = chunk->s + windowLen * windowIndex;
                int e = min(chunk->e, s + windowLen - 1);
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
    }
    ChunksCreator_destruct(chunksCreator);
    return correct;
}


int main(int argc, char *argv[]) {

    bool allTestsPassed = true;

    // test 1
    bool test1Passed = testCreatingChunks("tests/test_files/chunks_creator/test_1.cov");
    printf("[chunks_creator] Test creating chunks from uncompressed file:");
    printf(test1Passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");
    allTestsPassed &= test1Passed;

    // test 2
    bool test2Passed = testCreatingChunks("tests/test_files/chunks_creator/test_1.cov.gz");
    printf("[chunks_creator] Test creating chunks from compressed file:");
    printf(test2Passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");
    allTestsPassed &= test2Passed;

    // test 3
    bool test3Passed = testCreatingChunksWithLabels("tests/test_files/chunks_creator/test_1_with_labels.cov", false);
    printf("[chunks_creator] Test creating chunks from uncompressed file with labels:");
    printf(test3Passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");
    allTestsPassed &= test3Passed;

    // test 3
    bool test4Passed = testCreatingChunksWithLabels("tests/test_files/chunks_creator/test_1_with_labels.cov", true);
    printf("[chunks_creator] Test creating chunks from uncompressed file with labels (also checked writing and reading):");
    printf(test4Passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");
    allTestsPassed &= test4Passed;


    if (allTestsPassed)
        return 0;
    else
        return 1;
}
