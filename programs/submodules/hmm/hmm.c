#include <math.h>
#include "float.h"
#include <assert.h>
#include "../data_types/data_types.h"
#include "hmm.h"
#include <stdio.h>
#include "../block_it/block_it.h"
#include "../common/common.h"
#include "../ptBlock/ptBlock.h"
#include "../hmm_utils/hmm_utils.h"
#include "sonLib.h"
#include "stdlib.h"


#define TRANSITION_INITIAL_DIAG_PROB 0.9
static pthread_mutex_t chunkMutex;

///////////////////
// HMM Functions //
///////////////////
HMM *HMM_construct(int numberOfStates,
                   int numberOfRegions,
                   int *numberOfCompsPerState,
                   double **means,
                   double *meanScalePerRegion,
                   double maxHighMapqRatio,
                   double minHighlyClippedRatio,
                   char* pathToTransitionCounts,
                   ModelType modelType,
                   double* alpha){
    HMM *model = malloc(1 * sizeof(HMM));
    model->emissionDistSeriesPerRegion = malloc(numberOfRegions * sizeof(EmissionDistSeries *));
    model->transitionPerRegion =  malloc(numberOfRegions * sizeof(Transition *));
    int maxNumberOfComps = Double_getMaxValue1DArray(numberOfCompsPerState, numberOfStates);
    for(int region=0; region < numberOfRegions; region++){
        // construct emissionDistSeries
        double **meansForRegion = Double_copy2DArray(means, numberOfStates, maxNumberOfComps);
        Double_multiply2DArray(meansForRegion,
                               numberOfStates,
                               maxNumberOfComps,
                               meanScalePerRegion[region]);
        model->emissionDistSeriesPerRegion[region] = EmissionDistSeries_constructForModel(modelType,
                                                                                          means,
                                                                                          numberOfCompsPerState,
                                                                                          numberOfStates);
        Double_destruct2DArray(meansForRegion);

        // construct transition
        model->transitionPerRegion[region] = Transition_constructSymmetricBiased(numberOfStates,
                                                                                 TRANSITION_INITIAL_DIAG_PROB);
        Transition_addRequirements(model->transitionPerRegion[region],
                                   TransitionRequirements_construct(minHighlyClippedRatio, maxHighMapqRatio));
        Transition_addValidityFunction(model->transitionPerRegion[region],
                                       ValidityFunction_checkDupByMapq);
        Transition_addValidityFunction(model->transitionPerRegion[region],
                                       ValidityFunction_checkMsjByClipping);
    }
    model->alpha = MatrixDouble_copy(alpha);
    model->numberOfStates = numberOfStates;
    model->numberOfRegions = numberOfRegions;
    model->maxNumberOfComps = maxNumberOfComps;
    model->modelType = modelType;
    return model;
}

void HMM_destruct(HMM *model){
    for(int region=0; region < model->numberOfRegions; region++) {
        EmissionDistSeries_destruct(model->emissionDistSeriesPerRegion[region]);
        Transition_destruct(model->transitionPerRegion[region]);
    }
    free(model->emissionDistSeriesPerRegion);
    free(model->transitionPerRegion);
    MatrixDouble_destruct(model->alpha);

}

///////////////////////////////////////
// Functions for forward algorithm   //
//////////////////////////////////////

void EM_resetAllColumnsForward(EM *em){
    HMM *model = em->model;
    // Initialize to zero
    for (int i = 0; i < em->seqLen; i++) {
        for (int s = 0; s < model->numberOfStates; s++) {
            em->f[i][s] = 0.0;
        }
    }
}

void EM_fillFirstColumnForward(EM *em){
    HMM *model = em->model;
    double scale = 0.0; // For scaling the forward probabilities at each location to avoid underflow
    double eProb;
    double tProb;
    uint8_t preX = 0;
    double alpha = 0.0; // for the first column alpha can not be greater than 0
    // Set the 0-th block of the forward matrix
    scale = 0.0;
    for (int state = 0; state < model->numberOfStates; state++) {
        uint8_t region = em->coverageInfoSeq[0]->annotation_flag;
        uint8_t x = em->coverageInfoSeq[0]->coverage;
        // Emission probability
        eProb = EmissionDistSeries_getProb(model->emissionDistSeries[region],
                                           state,
                                           x,
                                           preX,
                                           alpha);
        // Transition probability
        tProb = Transition_getProbConditional(model->transitionPerRegion[region],
                                              START_STATE_INDEX,
                                              state,
                                              em->coverageInfoSeq[0]);
        // Update forward
        em->f[0][state] = eProb * tProb;
        scale += em->f[0][state];
    }
    em->scales[0] = scale;
    // Scale f 0-th column
    for (int state = 0; state < model->numberOfStates; state++) {
        em->f[0][state] /= scale;
    }
}

void EM_fillOneColumnForward(EM* em, int columnIndex){
    if (columnIndex == 0){
        EM_fillFirstColumnForward(em);
        return;
    }
    HMM *model = em->model;
    int i = columnIndex;
    double eProb;
    double tProb;
    uint8_t region;
    uint8_t preRegion;
    uint8_t x;
    uint8_t preX;
    double alpha;
    for (int state = 0; state < model->numberOfStates; state++) {
        for (int preState = 0; preState < model->numberOfSates; preState++) { // Transition from c1 comp to c2 comp
            region = em->coverageInfoSeq[i]->annotation_flag;
            preRegion = em->coverageInfoSeq[i-1]->annotation_flag;
            x = em->coverageInfoSeq[i]->coverage;
            preX = em->coverageInfoSeq[i-1]->coverage;
            alpha = model->alpha->data[preState][state];
            // Emission probability
            // Not that alpha can be zero and in that case emission probability is not
            // dependent on the previous observation
            eProb = EmissionDistSeries_getProb(model->emissionDistSeriesPerRegion[region],
                                               state,
                                               x,
                                               preX,
                                               model->alpha[preState][state]);
            if (region != preRegion) { // if the region class has changed
                // Make the transition prob uniform
                tProb = 1.0 / (model->numberOfStates + 1);
            } else {
                tProb = Transition_getProbConditional(model->transitionPerRegion[region],
                                                      preState,
                                                      state,
                                                      em->coverageInfoSeq[i]);
            }
            em->f[i][state] += (em->f[i - 1][preState] * tProb * eProb);
        }
        em->scales[i] += em->f[i][state];
    }
    if (em->scales[i] < 1e-50) {
        fprintf(stderr, "scale (= %.2e) is very low!\n", scale);
        exit(EXIT_FAILURE);
    }
    // Scale f
    for (int state = 0; state < model->numberOfStates; state++) {
        em->f[i][state] /= em->scales[i];
    }
}


void EM_runForward(EM *em) {
    EM_resetAllColumnsForward(em);
    // Fill columns of the forward matrix
    for (int columnIndex = 0; columnIndex < em->seqLen; columnIndex++) {
        EM_fillOneColumnForward(em, columnIndex);
    }
    // Update P(x)
    // TODO: P(x) is not calculated here p(x) = s(1) x s(2) x ... x s(L) check Durbin's
    // It shouldn't make any problem in the training process since px is not needed in the scaled calculations
    em->px = 0;
}

///////////////////////////////////////
// Functions for backward algorithm //
//////////////////////////////////////


void EM_resetAllColumnsBackward(EM *em){
    HMM *model = em->model;
    // Initialize to zero
    for (int i = 0; i < em->seqLen; i++) {
        for (int s = 0; s < model->numberOfStates; s++) {
            em->b[i][s] = 0.0;
        }
    }
}


void EM_fillLastColumnBackward(EM *em){
    HMM *model = em->model;
    double tProb;
    // Set the last block of the backward matrix
    for (int state = 0; state < model->numberOfStates; state++) {
        uint8_t region = em->coverageInfoSeq[em->seqLen - 1]->annotation_flag;
        // Transition probability
        tProb = Transition_getProbConditional(model->transitionPerRegion[region],
                                              state,
                                              END_STATE_INDEX,
                                              em->coverageInfoSeq[em->seqLen - 1]);
        // Update backward
        em->b[em->seqLen - 1][state] = tProb;
    }
    // Scale last column of the backward matrix
    for (int state = 0; state < model->numberOfStates; state++) {
        em->b[em->seqLen - 1][state] /= em->scales[em->seqLen - 1];
    }
}


void EM_fillOneColumnBackward(EM* em, int columnIndex){
    if (columnIndex == em->seqLen - 1){
        EM_fillLastColumnBackward(em);
        return;
    }
    HMM *model = em->model;
    int i = columnIndex;
    double eProb;
    double tProb;
    uint8_t region;
    uint8_t preRegion;
    uint8_t x;
    uint8_t preX;
    CoverageInfo *covInfo;
    CoverageInfo *preCovInfo;
    double alpha;
    for (int state = 0; state < model->numberOfStates; state++) {
        for (int preState = 0; preState < model->numberOfSates; preState++) { // Transition from c1 comp to c2 comp
            region = em->coverageInfoSeq[i + 1]->annotation_flag;
            preRegion = em->coverageInfoSeq[i]->annotation_flag;
            covInfo = em->coverageInfoSeq[i + 1];
            preCovInfo = em->coverageInfoSeq[i];
            preX = covInfo->coverage;
            alpha = model->alpha[preState][state];
            // Emission probability
            // Not that alpha can be zero and in that case emission probability is not
            // dependent on the previous observation
            eProb = EmissionDistSeries_getProb(model->emissionDistSeriesPerRegion[region],
                                               state,
                                               x,
                                               preX,
                                               alpha);
            if (region != preRegion) { // if the region class has changed
                // Make the transition prob uniform
                tProb = 1.0 / (model->numberOfStates + 1);
            } else {
                tProb = Transition_getProbConditional(model->transitionPerRegion[region],
                                                      preState,
                                                      state,
                                                      covInfo);
            }
            em->b[i][preState] += tProb * eProb * em->b[i + 1][state];
        }
    }
    if (em->scales[i] < 1e-50) {
        fprintf(stderr, "scale (= %.2e) is very low!\n", scale);
        exit(EXIT_FAILURE);
    }
    // Scale f
    for (int state = 0; state < model->numberOfStates; state++) {
        em->b[i][state] /= em->scales[i];
    }
}


// EM_runForward should be run before EM_runBackward because of scales
// after running EM_runForward scales will be saved in em->scales
// so they can be used in the backward algorithm
void EM_runBackward(EM *em) {
    EM_resetAllColumnsBackward(em);
    // Fill columns of the backward matrix
    for (int columnIndex = 0; columnIndex < em->seqLen; columnIndex++) {
        EM_fillOneColumnBackward(em, columnIndex);
    }
    // Update P(x)
    // TODO: P(x) is not calculated here p(x) = s(1) x s(2) x ... x s(L) check Durbin's
    // It shouldn't make any problem in the training process since px is not needed in the scaled calculations
    em->px = 0;
}



void EM_updateEstimatorsUsingOneColumn(EM* em, int columnIndex){
    if (columnIndex == em->seqLen - 1){
        return;
    }
    HMM *model = em->model;
    int i = columnIndex;
    double eProb;
    double tProb;
    uint8_t region;
    uint8_t preRegion;
    uint8_t x;
    uint8_t preX;
    double alpha;
    EmissionDistSeries *emissionDistSeries;
    Transition *transition;
    for (int state = 0; state < model->numberOfStates; state++) {
        for (int preState = 0; preState < model->numberOfSates; preState++) { // Transition from c1 comp to c2 comp
            // get observations
            region = em->coverageInfoSeq[i + 1]->annotation_flag;
            preRegion = em->coverageInfoSeq[i]->annotation_flag;
            x = em->coverageInfoSeq[i + 1]->coverage;
            preX = em->coverageInfoSeq[i]->coverage;
            // get model attributes
            alpha = model->alpha->data[preState][state];
            emissionDistSeries = model->emissionDistSeriesPerRegion[region];
            transition = model->transitionPerRegion[region];
            // Emission probability
            // Not that alpha can be zero and in that case emission probability is not
            // dependent on the previous observation
            eProb = EmissionDistSeries_getProb(emissionDistSeries,
                                               state,
                                               x,
                                               preX,
                                               alpha);
            if (region != preRegion) { // if the region class has changed
                // Make the transition prob uniform
                tProb = 1.0 / (model->numberOfStates + 1);
            } else {
                tProb = Transition_getProbConditional(transition,
                                                      preState,
                                                      state,
                                                      em->coverageInfoSeq[i + 1]);
            }
            // P(s_i = preState, s_(i+1) = state|x)
            double count = em->f[i][preState] * tProb * eProb * em->b[i + 1][state];
            EmissionDistSeries_updateEstimators(emissionDistSeries,
                                                x,
                                                preX,
                                                alpha,
                                                count);
            TransitionCountData_increment(transition->transitionCountData,
                                          count,
                                          preState,
                                          state);
        }
    }
}

void EM_updateEstimators(EM *em) {
    // skip first column since alpha might be > 0
    for(int columnIndex=1; columnIndex < em->seqLen; columnIndex++){
        EM_updateEstimatorsUsingOneColumn(em, columnIndex);
    }
}

void EM_estimateParameters(EM* em) {
    HMM* model = em->model;
}
