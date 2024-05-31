//
// Created by mobin on 5/30/24.
//

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
#include "cov_fast_reader.h"

CovFastReaderPerThread *CovFastReaderPerThread_construct(CovFastReader *covFastReader, int chunkIndexToParse) {
    CovFastReaderPerThread *covFastReaderPerThread = malloc(sizeof(CovFastReaderPerThread));
    covFastReaderPerThread->covFastReader = covFastReader;
    covFastReaderPerThread->chunkIndexToParse = chunkIndexToParse;
    return covFastReaderPerThread;
}

void CovFastReaderPerThread_destruct(CovFastReaderPerThread *covFastReaderPerThread) {
    free(covFastReaderPerThread);
}

CovFastReader *CovFastReader_construct(char *covPath, int chunkLen, int threads) {
    CovFastReader *covFastReader = malloc(sizeof(CovFastReader));
    covFastReader->chunksCreator = ChunksCreator_constructFromCov(covPath, NULL, chunkLen, threads, chunkLen);
    covFastReader->blockTablePerContig = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL,
                                                           (void (*)(void *)) stList_destruct);
    covFastReader->threads = threads;
    covFastReader->mutex = malloc(sizeof(pthread_mutex_t));
    pthread_mutex_init(covFastReader->mutex, NULL);

    CovFastReader_parseBlocks(covFastReader);
    return covFastReader;
}

stHash *CovFastReader_getBlockTablePerContig(CovFastReader *covFastReader) {
    return covFastReader->blockTablePerContig;
}

void CovFastReader_parseBlocks(CovFastReader *covFastReader) {
    int numberOfChunks = stList_length(covFastReader->chunksCreator->templateChunks);
    stList *covFastReaderPerThreads = stList_construct3(0, CovFastReaderPerThread_destruct);
    tpool_t *tm = tpool_create(covFastReader->threads);
    for (int chunkIndex = 0; chunkIndex < numberOfChunks; chunkIndex++) {
        CovFastReaderPerThread *covFastReaderPerThread =
                CovFastReaderPerThread_construct(covFastReader, chunkIndex);
        stList_append(covFastReaderPerThreads, covFastReaderPerThread);
        work_arg_t *argWork = malloc(sizeof(work_arg_t));
        argWork->data = (void *) covFastReaderPerThread;
        // queue job
        tpool_add_work(tm, CovFastReaderPerThread_parseBlocks, argWork);
	ChunksCreator *chunksCreator = covFastReader->chunksCreator;
	Chunk *templateChunk = stList_get(covFastReader->chunksCreator->templateChunks, chunkIndex);
	fprintf(stderr, "[%s] Parsing cov: submitted job for parsing chunk index:%d/%d %s:%d-%d\n", get_timestamp(), chunkIndex , numberOfChunks, templateChunk->ctg, templateChunk->s, templateChunk->e);
    }
    // wait until all jobs are finished
    tpool_wait(tm);
    // destroy thread pool
    tpool_destroy(tm);

    stList_destruct(covFastReaderPerThreads);

    //sort
    ptBlock_sort_stHash_by_rfs(covFastReader->blockTablePerContig);
}

void CovFastReaderPerThread_parseBlocks(void *arg_) {
    work_arg_t *arg = arg_;
    CovFastReaderPerThread *covFastReaderPerThread = arg->data;
    CovFastReader *covFastReader = covFastReaderPerThread->covFastReader;
    ChunksCreator *chunksCreator = covFastReader->chunksCreator;
    CoverageHeader *header = chunksCreator->header;
    int chunkIndex = covFastReaderPerThread->chunkIndexToParse;
    Chunk *templateChunk = stList_get(covFastReader->chunksCreator->templateChunks, chunkIndex);

    // open track reader
    bool zeroBasedCoors = true;
    TrackReader *trackReader = TrackReader_construct(chunksCreator->covPath, NULL, zeroBasedCoors);
    // Open cov file and jump to the first pos of chunk
    TrackReader_setFilePosition(trackReader, templateChunk->fileOffset);

    //set contig name
    strcpy(trackReader->ctg, templateChunk->ctg);
    trackReader->ctgLen = templateChunk->ctgLen;

    stList *blocks = NULL;
    while (0 < TrackReader_next(trackReader)) {
	if(templateChunk->e < trackReader->s) break;
	if(strcmp(trackReader->ctg, templateChunk->ctg) != 0) break;
	if(trackReader->s < templateChunk->s) continue;
        // create a ptBlock based on the parsed track
        ptBlock *block = ptBlock_constructFromTrackReader(trackReader, header);

        pthread_mutex_lock(covFastReader->mutex);
        // add block to the block stHash table
        blocks = stHash_search(covFastReader->blockTablePerContig, trackReader->ctg);
        // add new contig key to the table if it does not exist
        if (blocks == NULL) {
            blocks = stList_construct3(0, ptBlock_destruct);
            stHash_insert(covFastReader->blockTablePerContig, copyString(trackReader->ctg), blocks);
        }
        stList_append(blocks, block);
        pthread_mutex_unlock(covFastReader->mutex);
    }
    TrackReader_destruct(trackReader);
}

void CovFastReader_destruct(CovFastReader *covFastReader) {
    if (covFastReader->chunksCreator != NULL) {
        ChunksCreator_destruct(covFastReader->chunksCreator);
    }
    if (covFastReader->blockTablePerContig != NULL) {
        stHash_destruct(covFastReader->blockTablePerContig);
    }
    pthread_mutex_destroy(covFastReader->mutex);
    free(covFastReader);
}



