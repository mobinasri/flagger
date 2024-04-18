#ifndef HMM_H
#define HMM_H

#include "track_reader.h"
#include "common.h"
#include "data_types.h"
#include <stdio.h>
#include "sonLib.h"
#include "math.h"
#include "digamma.h"
#include "hmm_utils.h"



typedef struct HMM {
    EmissionDistSeries **emissionDistSeriesPerRegion;
    Transition **transitionPerRegion;
    ModelType modelType;
    MatrixDouble *alpha; // (#nComps + 1) x (#nComps + 1)
    // Dimension attributes
    int numberOfRegions; // Number of classes like non-HSAT, HSAT1, ...
    int numberOfStates; // Number of states like erroneous, haploid, ...
    int maxNumberOfComps;
    bool excludeMisjoin;
} HMM;

HMM *HMM_construct(int numberOfStates,
                   int numberOfRegions,
                   int *numberOfCompsPerState,
                   double **means,
                   double *meanScalePerRegion,
                   double maxHighMapqRatio,
                   double minHighlyClippedRatio,
                   char* pathToTransitionCounts,
                   ModelType modelType,
                   MatrixDouble *alpha,
		   bool excludeMisjoin);

void HMM_destruct(HMM *model);

int HMM_getStartStateIndex(HMM *model);
int HMM_getEndStateIndex(HMM *model);

void HMM_printTransitionMatrixInTsvFormat(HMM* model, FILE* fout);
void HMM_printEmissionParametersInTsvFormat(HMM* model, FILE* fout);

typedef struct EM {
    CoverageInfo **coverageInfoSeq; // the sequence of emissions
    int seqLen;
    double **f; // Forward matrix: #SEQ_LEN x #nComps
    double **b; // Backward matrix: #SEQ_LEN x #nComps
    double px; // P(x)
    double *scales;
    HMM *model;
} EM;

/*
typedef struct Chunk {
    // 2 * chunkLen is the maximum size for a chunk
    // the last chunk of a contig is most of the times
    // longer than chunkLen and shorter than 2 * chunkLen
    CoverageInfo **coverageInfoSeq; // [2 * chunkLen] the emitted sequence
    uint8_t *regionClassSeq; // [2 * chunkLen] the sequence of class of regions
    char ctg[200]; // the name of the contig this chunk is located
    int ctgLen; // the length of the contig
    int s; // start location of the chunk on the contig 0-based
    int e; // end location of the chunk on the contig 0-based
    int seqLen; // the length of the sequence that has been added to this chunk so far
    int chunkLen; // size of actual bases
    int maxSeqSize; // maximum size of array int(chunkLen * 2 / windowLen) + 1
    int windowLen;
    int windowItr;
    CoverageInfo *coverageInfoWindowSum;
    uint8_t *windowClassSeq;
    uint64_t fileOffset; // the number of offset bytes to reach the first block of this chunk
} Chunk;


typedef struct Batch {
    // Constant attributes
    int nThreads;
    int chunkLen;
    int windowLen;
    // Attributes that may change
    char covPath[200];
    Chunk **threadChunks; // There are at most nThreads chunks in the batch
    stList *templateChunks; // The list of all chunks (no seq) in the cov file
    int templateChunkIdx;
    int nThreadChunks; // The number of chunks added to this batch so far
    pthread_mutex_t *mutex;
} Batch;

stList* Chunk_readAllChunksFromBin(char* covPath, int chunkLen, int windowLen, int nEmit);
*/

EM *EM_construct(CoverageInfo **coverageInfoSeq, int seqLen, HMM *model);

void EM_destruct(EM *em);

void EM_runForward(EM *em);

void EM_runBackward(EM *em);

void EM_updateEstimators(EM *em);

bool EM_estimateParameters(EM *em, double convergenceTol);

void EM_resetEstimators(EM *em);

double *EM_getPosterior(EM *em, int pos);

int EM_getMostProbableState(EM *em, int pos);

void EM_printPosteriorInTsvFormat(EM* em, FILE* fout);

/*
double *getForward(EM *em, int pos);

double *getBackward(EM *em, int pos);

double *getPosterior(EM *em, int pos);

Chunk *Chunk_construct1(int chunkLen);

Chunk *Chunk_construct3(int chunkLen, int emissionDim, int windowLen);

stList *createCovIndex(char *covPath, int chunkLen);

stList *parseCovIndex(char *covIndexPath);

void writeCovIndex(stList *chunks, char *indexPath);

int Chunk_addBlock(Chunk *chunk, Block_t *block);

Batch *Batch_construct(char *covPath, int chunkLen, int nThreads, int emissionDim, int windowLen);

void Batch_destruct(Batch *batch);

int Batch_readThreadChunks(Batch *batch);

int Batch_readNextChunk(void *batch_);
*/

#endif
