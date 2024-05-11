#include "track_reader.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

bool testParsingCov(char *covPath) {
    bool correct = true;
    int truthValuesCtg1[7][6] = {{1,  10,  4,  4,  4,  0},
                                 {11, 15,  6,  0,  0,  0},
                                 {16, 20,  6,  0,  0,  1},
                                 {21, 60,  10, 10, 10, 1},
                                 {61, 64,  14, 14, 10, 1},
                                 {65, 70,  14, 14, 10, 0},
                                 {71, 110, 16, 16, 16, 1}};
    int truthValuesCtg2[2][6] = {{1, 2,  4, 4, 4, 1},
                                 {3, 10, 8, 8, 0, 0}};
    int trackIndexCtg1 = -1;
    int trackIndexCtg2 = -1;
    bool zeroBasedCoors = false;
    TrackReader *trackReader = TrackReader_construct(covPath, NULL, zeroBasedCoors);
    fprintf(stderr, "1\n");
    while (0 < TrackReader_next(trackReader)) {
        fprintf(stderr, "2\n");
        fprintf(stderr, "%s:%d-%d\n", trackReader->ctg, trackReader->s, trackReader->e);
        if (strcmp(trackReader->ctg, "ctg1") == 0) {
            trackIndexCtg1 += 1;
            int *truthValues = truthValuesCtg1[trackIndexCtg1];
            if (trackReader->attrbsLen != 5) {
                fprintf(stderr, "Number of parsed attributes per track %d does not match truth (5)\n",
                        trackReader->attrbsLen);
                TrackReader_destruct(trackReader);
                return false;
            }
            correct &= (trackReader->s == truthValues[0]);
            correct &= (trackReader->e == truthValues[1]);
            correct &= (atoi(trackReader->attrbs[0]) == truthValues[2]);
            correct &= (atoi(trackReader->attrbs[1]) == truthValues[3]);
            correct &= (atoi(trackReader->attrbs[2]) == truthValues[4]); // skip annotation column for now
            correct &= (atoi(trackReader->attrbs[4]) == truthValues[5]);
        }
        if (strcmp(trackReader->ctg, "ctg2") == 0) {
            trackIndexCtg2 += 1;
            int *truthValues = truthValuesCtg2[trackIndexCtg2];
            if (trackReader->attrbsLen != 5) {
                fprintf(stderr, "Number of parsed attributes per track %d does not match truth (5)\n",
                        trackReader->attrbsLen);
                TrackReader_destruct(trackReader);
                return false;
            }
            correct &= (trackReader->s == truthValues[0]);
            correct &= (trackReader->e == truthValues[1]);
            correct &= (atoi(trackReader->attrbs[0]) == truthValues[2]);
            correct &= (atoi(trackReader->attrbs[1]) == truthValues[3]);
            correct &= (atoi(trackReader->attrbs[2]) == truthValues[4]);
            correct &= (atoi(trackReader->attrbs[4]) == truthValues[5]);
        }
        fprintf(stderr, "**2\n");
    }
    fprintf(stderr, "***2\n");
    TrackReader_destruct(trackReader);
    return correct;
}


int main(int argc, char *argv[]) {
    char covPath[1000] = "tests/test_files/test_1.cov";
    char covGzPath[1000] = "tests/test_files/test_1.cov.gz";
    if (testParsingCov(covPath) == true) {
        fprintf(stderr, "Test for tests/test_files/test_1.cov passed!\n");
    } else {
        fprintf(stderr, "Test for tests/test_files/test_1.cov failed!\n");
        return 1;
    }

    if (testParsingCov(covGzPath) == true) {
        fprintf(stderr, "Test for tests/test_files/test_1.cov.gz passed!\n");
    } else {
        fprintf(stderr, "Test for tests/test_files/test_1.cov.gz failed!\n");
        return 1;
    }
    return 0;

}
