//
// Created by mobin on 5/30/24.
//

#ifndef FLAGGER_COV_FAST_READER_H
#define FLAGGER_COV_FAST_READER_H

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

typedef struct CovFastReader {
    ChunksCreator *chunksCreator;
    stHash *blockTablePerContig;
    int threads;
    pthread_mutex_t *mutex;
} CovFastReader;

CovFastReader *CovFastReader_construct(char *covPath, int chunkLen, int threads);

stHash *CovFastReader_getBlockTablePerContig(CovFastReader *covFastReader);

void CovFastReader_parseBlocks(CovFastReader *covFastReader);

void CovFastReaderPerThread_parseBlocks(void *arg_);

void CovFastReader_destruct(CovFastReader *covFastReader);

typedef struct CovFastReaderPerThread {
    CovFastReader *covFastReader;
    int chunkIndexToParse;
} CovFastReaderPerThread;

CovFastReaderPerThread *CovFastReaderPerThread_construct(CovFastReader *covFastReader, int chunkIndexToParse);

void CovFastReaderPerThread_destruct(CovFastReaderPerThread *covFastReaderPerThread);

#endif //FLAGGER_COV_FAST_READER_H
