#include <math.h>
#include "float.h"
#include <assert.h>
#include "data_types.h"
#include "hmm.h"
#include <stdio.h>
#include "block_it.h"
#include "common.h"
#include "ptBlock.h"
#include "hmm_utils.h"
#include "sonLib.h"
#include "stdlib.h"


#define TRANSITION_INITIAL_DIAG_PROB 0.9
#define TRANSITION_PSEUDO_COUNT_VALUE 10
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
                   MatrixDouble* alpha){
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
                                                                                          meansForRegion,
                                                                                          numberOfCompsPerState,
                                                                                          numberOfStates);
        Double_destruct2DArray(meansForRegion);

        // construct transition
        model->transitionPerRegion[region] = Transition_constructSymmetricBiased(numberOfStates,
                                                                                 TRANSITION_INITIAL_DIAG_PROB);
        TransitionCountData_setPseudoCountMatrix(model->transitionPerRegion[region]->transitionCountData,
                                                 TRANSITION_PSEUDO_COUNT_VALUE);
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

void HMM_estimateParameters(HMM* model) {
    for(int region = 0; region < model->numberOfRegions; region++){
        EmissionDistSeries_estimateParameters(model->emissionDistSeriesPerRegion[region]);
        Transition_estimateTransitionMatrix(model->transitionPerRegion[region]);
    }
}

void HMM_printTransitionMatrixInTsvFormat(HMM* model, FILE* fout){
    // 100 = maximum number of rows
    // 10 = maximum number of columns
    // 20 = maximum size of each entry/word
    char table[100][10][20];
    int numberOfColumns = model->numberOfStates + 3; // +3; because "Region" + "PreState/State" and end state
    // fill first row (header)
    sprintf(table[0][0], "#Region");
    sprintf(table[0][1], "State");
    for(int state=0; state < model->numberOfStates; state++) {
        const char* stateName = EmissionDistSeries_getStateName(state);
        sprintf(table[0][2 + state],"%s", stateName);
    }
    sprintf(table[0][2 + model->numberOfStates], "End");

    // write a transition matrix for each region one at a time
    int numberOfRows = 1;
    for(int region=0; region < model->numberOfRegions; region++) {
        Transition *transition = model->transitionPerRegion[region];
        // fill the first and second columns for this region
        for(int state=0; state < model->numberOfStates; state++) {
            sprintf(table[numberOfRows + state][0], "%d", region);
            const char* stateName = EmissionDistSeries_getStateName(state);
            sprintf(table[numberOfRows + state][1], "%s", stateName);
        }
        sprintf(table[numberOfRows + model->numberOfStates][0], "%d", region);
        sprintf(table[numberOfRows + model->numberOfStates][1],"Start");
        // fill the probability matrix
        for (int preState = 0; preState < model->numberOfStates + 1; preState++) {
            for (int state = 0; state < model->numberOfStates + 1; state++) {
                sprintf(table[numberOfRows + preState][2 + state],
                        "%.2e", transition->matrix->data[preState][state]);
            }
        }
        numberOfRows += model->numberOfStates + 1; // +1 for "Start" state
    }

    for(int row=0; row < numberOfRows; row++){
        for(int col=0; col < numberOfColumns - 1; col++){
            fwrite(fout, "%s\t", table[row][col]);
        }
        // last column needs \n instead of \t
        fwrite(fout, "%s\n", table[row][numberOfColumns - 1]);
    }
}

void HMM_printEmissionParametersInTsvFormat(HMM* model, FILE* fout){
    // 100 = maximum number of rows
    // 40 = maximum number of columns
    // 200 = maximum size of each entry/word
    char table[100][40][200];
    // 4 is for the first 4 columns; State	Distribution	Components	Parameter
    int numberOfColumns = 4 + model->numberOfRegions;
    sprintf(table[0][0], "#State");
    sprintf(table[0][1], "Distribution");
    sprintf(table[0][2], "Components");
    sprintf(table[0][3], "Parameter");
    for(int region=0; region < model->numberOfRegions; region++) {
        sprintf(table[0][4 + region], "Values_Region_%d", region);
    }
    int numberOfRows=1;
    for(int region=0; region < model->numberOfRegions; region++){
        EmissionDistSeries *emissionDistSeries = model->emissionDistSeriesPerRegion[region];
        for(int state=0; state < emissionDistSeries->numberOfDists; state++) {
            int numberOfParams;
            const char* stateName = EmissionDistSeries_getStateName(state);
            const char* distributionName = EmissionDistSeries_getDistributionName(emissionDistSeries, state);
            const char** parameterNames = EmissionDistSeries_getParameterNames(emissionDistSeries, state, &numberOfParams);
            double **parameterValues = EmissionDistSeries_getParameterValues(emissionDistSeries, state, &numberOfParams);
            int numberOfComps = EmissionDistSeries_getNumberOfComps(emissionDistSeries, state);
            for(int param=0; param < numberOfParams; param++){
                char * parameterValuesStr = String_joinDoubleArray(parameterValues[param], numberOfComps, ',');
                sprintf(table[numberOfRows][0], "%s", stateName);
                sprintf(table[numberOfRows][1], "%s", distributionName);
                sprintf(table[numberOfRows][2], "%d", numberOfComps);
                sprintf(table[numberOfRows][3], "%s", parameterNames[param]);
                sprintf(table[numberOfRows][4 + region], "%s", parameterValuesStr);
                free(parameterValuesStr);
                numberOfRows += 1;
            }
            free(parameterValues);
        }
    }

    for(int row=0; row < numberOfRows; row++){
        for(int col=0; col < numberOfColumns - 1; col++){
            fwrite(fout, "%s\t", table[row][col]);
        }
        // last column needs \n instead of \t
        fwrite(fout, "%s\n", table[row][numberOfColumns - 1]);
    }
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


EM *EM_construct(CoverageInfo **coverageInfoSeq, int seqLen, HMM *model) {
    assert(seqLen > 0);
    EM *em = malloc(sizeof(EM));
    em->coverageInfoSeq = coverageInfoSeq;
    em->seqLen = seqLen;
    em->model = model;
    // Allocate and initialize forward and backward matrices
    em->f = Double_construct2DArray(em->seqLen, model->numberOfStates);
    em->b = Double_construct2DArray(em->seqLen, model->numberOfStates);
    // Initialize scale to avoid underflow
    em->scales = Double_construct1DArray(em->seqLen);
    em->px = -1.0;
    return em;
}

void EM_destruct(EM *em) {
    Double_construct2DArray(em->f, em->seqLen);
    Double_construct2DArray(em->b, em->seqLen);
    Double_construct1DArray(em->scales);
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
            EmissionDistSeries_updateEstimator(emissionDistSeries,
                                                state,
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

void EM_estimateParameters(EM *em) {
    HMM *model = em->model;
    for(int region=0; region < model->numberOfRegions; region++) {
        EmissionDistSeries_estimateParameters(model->emissionDistSeriesPerRegion[region]);
        Transition_estimateTransitionMatrix(model->transitionPerRegion[region]);
    }
}
