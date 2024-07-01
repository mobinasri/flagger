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
    double loglikelihood;
} HMM;

HMM *HMM_construct(int numberOfStates,
                   int numberOfRegions,
                   int *numberOfCompsPerState,
                   double **means,
                   double *meanScalePerRegion,
                   double maxHighMapqRatio,
                   double minHighlyClippedRatio,
                   char *pathToTransitionCounts,
                   ModelType modelType,
                   MatrixDouble *alpha,
                   bool excludeMisjoin);

bool HMM_isFeasible(HMM *model);

void HMM_normalizeWeightsAndTransitionRows(HMM *model);

HMM *HMM_copy(HMM *src);

void HMM_destruct(HMM *model);

int HMM_getStartStateIndex(HMM *model);

int HMM_getEndStateIndex(HMM *model);

bool HMM_estimateParameters(HMM *model, double convergenceTol);

void HMM_resetEstimators(HMM *model);

void HMM_printTransitionMatrixInTsvFormat(HMM *model, FILE *fout);

void HMM_printEmissionParametersInTsvFormat(HMM *model, FILE *fout);

typedef struct EM {
    CoverageInfo **coverageInfoSeq; // the sequence of emissions
    int seqLen;
    double **f; // Forward matrix: #SEQ_LEN x #nComps
    double **b; // Backward matrix: #SEQ_LEN x #nComps
    double px; // P(x)
    double *scales;
    HMM *model;
    EmissionDistSeries **emissionDistSeriesPerRegion;
    Transition **transitionPerRegion;
    int numberOfRegions;
    double loglikelihood;
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

void EM_renewParametersAndEstimatorsFromModel(EM *em, HMM *model);

void EM_destruct(EM *em);

void EM_runForward(EM *em);

void EM_runBackward(EM *em);

void EM_updateModelEstimators(EM *em);

void EM_updateEstimators(EM *em);

bool EM_estimateParameters(EM *em, double convergenceTol);

void EM_resetEstimators(EM *em);

double *EM_getPosterior(EM *em, int pos);

int EM_getMostProbableState(EM *em, int pos);

void EM_printPosteriorInTsvFormat(EM *em, FILE *fout);

void EM_runOneIterationAndUpdateEstimatorsForThreadPool(void *arg_);

void EM_runOneIterationForList(stList *emList, HMM *model, int threads);

void EM_runForwardForThreadPool(void *arg_);

void EM_runForwardForList(stList *emList, HMM *model, int threads);

typedef struct SquareAccelerator {
    HMM *model0;
    HMM *model1;
    HMM *model2;
    HMM *modelPrime;
    HMM *modelRatesV;
    HMM *modelRatesR;
    double alphaRate;
} SquareAccelerator;

// steps :
// 1. run SquareAccelerator_construct to make an accelerator object
// 2. set model 0 with SquareAccelerator_setModel0 (this function makes a copy of model 0)
// 3. run EM_runOneIterationForList on model 0 to get model 1
// 4. set model 1 with SquareAccelerator_setModel1 (this function makes a copy of model 1)
// 5. run EM_runOneIterationForList on model 1 to get model 2
// 6. set model 2 with SquareAccelerator_setModel2 (this function makes a copy of model 2)
// 7. call SquareAccelerator_getModelPrime to get model prime
// 8. run EM_runOneIterationForList on model prime to get the final model after acceleration

SquareAccelerator *SquareAccelerator_construct();
void SquareAccelerator_destruct(SquareAccelerator *accelerator);
void SquareAccelerator_setModel0(SquareAccelerator *accelerator, HMM *model0);
void SquareAccelerator_setModel1(SquareAccelerator *accelerator, HMM *model1);
void SquareAccelerator_setModel2(SquareAccelerator *accelerator, HMM *model2);
HMM *SquareAccelerator_getModelPrime(SquareAccelerator *accelerator, stList *emList, int threads);

// run SquareAccelerator_computeRates before running this function
HMM *SquareAccelerator_computeValuesForModelPrime(SquareAccelerator *accelerator);
void SquareAccelerator_computeRates(SquareAccelerator *accelerator);

#endif
