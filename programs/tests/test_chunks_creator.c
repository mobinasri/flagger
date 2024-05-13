#include "track_reader.h"
#include "chunk.h"
#include "sonLib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

bool testCreatingChunks(char *covPath) {
    bool correct = true;
    int truthValuesCtg1[6][6] = {{1,   20,  5,  2,  2,  0},
                                 {21,  40,  10, 10, 10, 1},
                                 {41,  60,  10, 10, 10, 1},
                                 {61,  80,  15, 15, 13, 1},
                                 {81,  100, 16, 16, 16, 1},
                                 {101, 110, 16, 16, 16, 1}};
    int truthValuesCtg2[1][6] = {{1, 10, 7, 7, 1, 0}};
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
        //fprintf(stderr,"### %s:%d-%d\n", chunk->ctg, chunk->s, chunk->e);
        if (strcmp(chunk->ctg, "ctg1") == 0) {
            for (int windowIndex = 0; windowIndex < chunk->coverageInfoSeqLen; windowIndex++) {
                trackIndexCtg1 += 1;
                CoverageInfo *coverageInfo = chunk->coverageInfoSeq[windowIndex];
                int *truthValues = truthValuesCtg1[trackIndexCtg1];
                //correct &= (chunk->s == (truthValues[0] - 1));
                //correct &= (chunk->e == (truthValues[1] - 1));
                fprintf(stderr, "**** %d %d %d %d\n", coverageInfo->coverage, coverageInfo->coverage_high_mapq,
                        coverageInfo->coverage_high_clip, coverageInfo->annotation_flag);
                correct &= (coverageInfo->coverage == truthValues[2]);
                correct &= (coverageInfo->coverage_high_mapq == truthValues[3]);
                correct &= (coverageInfo->coverage_high_clip == truthValues[4]);
                correct &= (CoverageInfo_getRegionIndex(coverageInfo) == truthValues[5]);
            }
        }
        if (strcmp(chunk->ctg, "ctg2") == 0) {
            for (int windowIndex = 0; windowIndex < chunk->coverageInfoSeqLen; windowIndex++) {
                trackIndexCtg2 += 1;
                CoverageInfo *coverageInfo = chunk->coverageInfoSeq[windowIndex];
                int *truthValues = truthValuesCtg2[trackIndexCtg2];
                //correct &= (chunk->s == (truthValues[0] - 1));
                //correct &= (chunk->e == (truthValues[1] - 1));
                fprintf(stderr, "**** %d %d %d %d\n", coverageInfo->coverage, coverageInfo->coverage_high_mapq,
                        coverageInfo->coverage_high_clip, coverageInfo->annotation_flag);
                correct &= (coverageInfo->coverage == truthValues[2]);
                correct &= (coverageInfo->coverage_high_mapq == truthValues[3]);
                correct &= (coverageInfo->coverage_high_clip == truthValues[4]);
                correct &= (CoverageInfo_getRegionIndex(coverageInfo) == truthValues[5]);
            }
        }
    }
    ChunksCreator_destruct(chunksCreator);
    return correct;
}



int main(int argc, char *argv[]) {
    // test 1
    if (testCreatingChunks("tests/test_files/test_1.cov") == true) {
        fprintf(stderr, "Test CreatingChunks for uncompressed file passed!\n");
    } else {
        fprintf(stderr, "Test CreatingChunks for uncompressed file failed!\n");
        return 1;
    }

    // test 2
    if (testCreatingChunks("tests/test_files/test_1.cov.gz") == true) {
        fprintf(stderr, "Test CreatingChunks for compressed file passed!\n");
    } else {
        fprintf(stderr, "Test CreatingChunks for compressed file failed!\n");
        return 1;
    }

    return 0;

}
