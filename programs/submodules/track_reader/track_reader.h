#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include "stdio.h"
#include "stdlib.h"
#include "stdbool.h"
#include "ptBlock.h"

#ifndef PT_TRACK_H
#define PT_TRACK_H

typedef enum TrackFileFormat {
    TRACK_FILE_FORMAT_COV,
    TRACK_FILE_FORMAT_COV_GZ,
    TRACK_FILE_FORMAT_BED,
    TRACK_FILE_FORMAT_BED_GZ,
    TRACK_MEMORY_COV,
    TRACK_FILE_FORMAT_UNDEFINED
} TrackFileFormat;

typedef struct TrackReader {
    TrackFileFormat trackFileFormat;
    stHash *contigLengthTable;
    void *fileReaderPtr;
    char ctg[1000];
    int ctgLen;
    int s; // 1-based (0-based if zeroBasedCoors is true)
    int e; // 1-based (0-based if zeroBasedCoors in true)
    char **attrbs;
    int attrbsLen;
    bool zeroBasedCoors;
    // attributes for iterating over coverage blocks in memory
    stList *contigList;
    stHash *coverageBlockTable; // would be NULL if reading from file
    stList *coverageBlockListBeingIterated; // would be NULL if reading from file
    int nextBlockIndexToRead; // would be -1 if reading from file
    int nextContigIndexToRead;
} TrackReader;


typedef struct CoverageHeader {
    stList *headerLines;
    int numberOfAnnotations;
    stList *annotationNames;
    int numberOfRegions;
    int *regionCoverages;
    stList *regionNames;
    int numberOfLabels;
    bool isTruthAvailable;
    bool isPredictionAvailable;
    bool startOnlyMode;
    int averageAlignmentLength;
} CoverageHeader;

CoverageHeader *CoverageHeader_construct(char *filePath);

void CoverageHeader_writeIntoFile(CoverageHeader *header,
                                  void *filePtr,
                                  bool isCompressed);

CoverageHeader *CoverageHeader_constructByAttributes(stList *annotationNames,
                                                     int *regionCoverages,
                                                     int numberOfRegions,
                                                     int numberOfLabels,
                                                     bool isTruthAvailable,
                                                     bool isPredictionAvailable,
                                                     bool startOnlyMode,
                                                     int averageAlignmentLength);

void CoverageHeader_destruct(CoverageHeader *header);

void CoverageHeader_updateNumberOfAnnotations(CoverageHeader *header);

void CoverageHeader_updateAnnotationNames(CoverageHeader *header);

void CoverageHeader_updateNumberOfRegions(CoverageHeader *header);

void CoverageHeader_updateRegionCoverages(CoverageHeader *header);

void CoverageHeader_updateRegionNames(CoverageHeader *header);

void CoverageHeader_updateNumberOfLabels(CoverageHeader *header);

void CoverageHeader_updateTruthAvailability(CoverageHeader *header);

void CoverageHeader_updatePredictionAvailability(CoverageHeader *header);

void CoverageHeader_updateAverageAlignmentLength(CoverageHeader *header);

void CoverageHeader_updateStartOnlyMode(CoverageHeader *header);

void *TrackReader_openFile(char *filePath, TrackFileFormat format);

stList *TrackReader_parseHeaderLines(TrackReader *trackReader);

TrackFileFormat TrackReader_getTrackFileFormat(char *filePath);

int TrackReader_readLine(TrackReader *trackReader, char **line, int maxSize);

int64_t TrackReader_getFilePosition(TrackReader *trackReader);

void TrackReader_setFilePosition(TrackReader *trackReader, int64_t filePosition);

TrackReader *TrackReader_construct(char *filePath, char *faiPath, bool zeroBasedCoors);

TrackReader *TrackReader_constructFromTableInMemory(stHash *blockTable, bool zeroBasedCoors);

void TrackReader_destruct(TrackReader *trackReader);

// Read next trackReader in BED or COV
int TrackReader_next(TrackReader *trackReader);

int TrackReader_readNextTrackBed(TrackReader *trackReader);

int TrackReader_readNextTrackCov(TrackReader *trackReader);

#endif
