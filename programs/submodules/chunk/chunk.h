#include <stdint.h>
#include "data_types.h"
#include "ptBlock.h"
#include "common.h"
#include "track_reader.h"


typedef struct Chunk {
    // 2 * chunkCanonicalLen is the maximum size for a chunk
    // the last chunk of a contig is most of the times
    // longer than chunkCanonicalLen and shorter than 2 * chunkCanonicalLen
    CoverageInfo **coverageInfoSeq; // [2 * chunkCanonicalLen] the emitted sequence
    char ctg[200]; // the name of the contig where this chunk is located in
    int ctgLen; // the length of the contig
    int s; // start location of the chunk on the contig 0-based
    int e; // end location of the chunk on the contig 0-based
    int coverageInfoSeqLen; // the length of the sequence that has been added to this chunk so far
    int chunkCanonicalLen; // size of actual bases
    int coverageInfoMaxSeqSize; // maximum size of array int(chunkCanonicalLen * 2 / windowLen) + 1
    int windowLen;
    int windowItr;
    CoverageInfo *windowSumCoverageInfo;
    uint64_t windowAnnotationFlag;
    int *windowRegionArray;
    int *windowTruthArray;
    int *windowPredictionArray;
    uint64_t fileOffset; // the number of offset bytes to reach the first trackReader of this chunk
} Chunk;

typedef struct ChunksCreator {
    char *covPath;
    // attributes that should be parsed from header lines
    int numberOfAnnotations;
    stList* annotationNames;
    int numberOfRegions;
    int *regionCoverages;
    int numberOfLabels;
    bool isTruthAvailable;
    // Constant attributes
    int nThreads;
    int chunkCanonicalLen;
    int windowLen;
    // Attributes that may change
    stList *chunks;
    stList *templateChunks; // The list of all chunks (no seq) in the cov file
    int nextChunkIndexToRead;
    pthread_mutex_t *mutex;
} ChunksCreator;


int Chunk_cmp(const void *chunk_1, const void *chunk2);
Chunk *Chunk_construct(int chunkCanonicalLen);
Chunk *Chunk_constructWithAllocatedSeq(int chunkCanonicalLen, int windowLen, int maxSeqSize);
void Chunk_destruct(Chunk *chunk);
stList *Chunk_constructListWithAllocatedSeq(stList *templateChunks, int windowLen);
int32_t Chunk_getWindowRegion(Chunk *chunk);
int Chunk_addWindow(Chunk *chunk);
int Chunk_addTrack(Chunk *chunk, TrackReader *trackReader);


ChunksCreator *ChunksCreator_constructEmpty();
ChunksCreator *ChunksCreator_constructFromCov(char *covPath, char *faiPath, int chunkCanonicalLen, int nThreads, int windowLen);
stList *ChunksCreator_createCovIndex(char *filePath, char *faiPath, int chunkCanonicalLen);
void ChunksCreator_writeCovIndex(stList *templateChunks, char *indexPath);
stList *ChunksCreator_parseCovIndex(char *covIndexPath);
void ChunksCreator_destruct(ChunksCreator *chunksCreator);
int ChunksCreator_parseChunks(ChunksCreator *chunksCreator);
void ChunksCreator_parseOneChunk(void *chunksCreator_);
void ChunksCreator_sortChunks(ChunksCreator *chunksCreator);

void ChunksCreator_parseAnnotationNames(ChunksCreator *chunksCreator);
void ChunksCreator_parseRegionCoverages(ChunksCreator *chunksCreator);
void ChunksCreator_parseNumberOfLabels(ChunksCreator *chunksCreator);
void ChunksCreator_parseTruthAvailability(ChunksCreator *chunksCreator);

void ChunksCreator_writeChunksIntoBedGraph(ChunksCreator *chunksCreator,
                                    const char *outputPath,
                                    const char *trackName,
                                    u_int16_t (*getCoverageInfoAttribute)(CoverageInfo *),
                                    const char *writingMode,
				    const char *color);

void ChunkCreator_parseChunksFromBinaryFile(ChunksCreator *chunksCreator, char *binPath);
void ChunksCreator_writeChunksIntoBinaryFile(ChunksCreator *chunksCreator, char *binPath);
