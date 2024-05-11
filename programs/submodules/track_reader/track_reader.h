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

typedef enum TrackFileFormat {TRACK_FILE_FORMAT_COV, TRACK_FILE_FORMAT_COV_GZ, TRACK_FILE_FORMAT_BED, TRACK_FILE_FORMAT_BED_GZ, TRACK_FILE_FORMAT_UNDEFINED} TrackFileFormat;

typedef struct TrackReader{
	TrackFileFormat trackFileFormat;
	stList *headerLines;
	stHash *contigLengthTable;
	void *fileReaderPtr;
        char ctg[1000];
        int ctgLen;
        int s; // 1-based (0-based if zeroBasedCoors is true)
        int e; // 1-based (0-based if zeroBasedCoors in true)
        char** attrbs;
        int attrbsLen;
	bool zeroBasedCoors;
} TrackReader;

void *TrackReader_openFile(char *filePath, TrackFileFormat format);

TrackFileFormat TrackReader_getTrackFileFormat(char *filePath);

int TrackReader_readLine(TrackReader *trackReader, char **line, int maxSize);

void TrackReader_updateHeaderLines(TrackReader *trackReader);

int64_t TrackReader_getFilePosition(TrackReader *trackReader);

void TrackReader_setFilePosition(TrackReader *trackReader, int64_t filePosition);

TrackReader* TrackReader_construct(char* filePath, char *faiPath, bool zeroBasedCoors);

void TrackReader_destruct(TrackReader* trackReader);

// Read next trackReader in BED or COV
int TrackReader_next(TrackReader* trackReader);

int TrackReader_readNextTrackBed(TrackReader* trackReader);

int TrackReader_readNextTrackCov(TrackReader* trackReader);

#endif
