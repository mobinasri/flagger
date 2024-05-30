#include <math.h>
#include "float.h"
#include <assert.h>
#include "data_types.h"
#include "hmm.h"
#include <stdio.h>
#include "track_reader.h"
#include "common.h"
#include "ptBlock.h"
#include "hmm_utils.h"
#include "sonLib.h"
#include "stdlib.h"


#define TRANSITION_INITIAL_DIAG_PROB 0.99
#define TRANSITION_PSEUDO_COUNT_VALUE 0.001
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
                   char *pathToTransitionCounts,
                   ModelType modelType,
                   MatrixDouble *alpha,
                   bool excludeMisjoin) {
    HMM *model = malloc(1 * sizeof(HMM));
    model->excludeMisjoin = excludeMisjoin;
    model->emissionDistSeriesPerRegion = malloc(numberOfRegions * sizeof(EmissionDistSeries *));
    model->transitionPerRegion = malloc(numberOfRegions * sizeof(Transition *));
    int maxNumberOfComps = Int_getMaxValue1DArray(numberOfCompsPerState, numberOfStates);
    for (int region = 0; region < numberOfRegions; region++) {
        fprintf(stdout, "region=%d\n", region);
        fprintf(stdout, "construct emissions\n");
        // construct emissionDistSeries
        double **meansForRegion = Double_copy2DArray(means, numberOfStates, maxNumberOfComps);
        Double_multiply2DArray(meansForRegion,
                               numberOfStates,
                               maxNumberOfComps,
                               meanScalePerRegion[region]);
        model->emissionDistSeriesPerRegion[region] = EmissionDistSeries_constructForModel(modelType,
                                                                                          meansForRegion,
                                                                                          numberOfCompsPerState,
                                                                                          numberOfStates,
                                                                                          excludeMisjoin);
        Double_destruct2DArray(meansForRegion, numberOfStates);

        fprintf(stdout, "construct transition\n");
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

int HMM_getStartStateIndex(HMM *model) {
    return model->excludeMisjoin ? START_STATE_INDEX - 1 : START_STATE_INDEX;
}

int HMM_getEndStateIndex(HMM *model) {
    return model->excludeMisjoin ? END_STATE_INDEX - 1 : END_STATE_INDEX;
}

bool HMM_estimateParameters(HMM *model, double convergenceTol) {
    bool converged = true;
    for (int region = 0; region < model->numberOfRegions; region++) {
        converged &= EmissionDistSeries_estimateParameters(model->emissionDistSeriesPerRegion[region], convergenceTol);
        converged &= Transition_estimateTransitionMatrix(model->transitionPerRegion[region], convergenceTol);
    }
    return converged;
}

void HMM_resetEstimators(HMM *model) {
    for (int region = 0; region < model->numberOfRegions; region++) {
        EmissionDistSeries_resetParameterEstimators(model->emissionDistSeriesPerRegion[region]);
        Transition_resetCountData(model->transitionPerRegion[region]);
    }
}


void HMM_printTransitionMatrixInTsvFormat(HMM *model, FILE *fout) {
    // 100 = maximum number of rows
    // 10 = maximum number of columns
    // 20 = maximum size of each entry/word
    char table[100][10][20];
    int numberOfColumns = model->numberOfStates + 3; // +3; because "Region" + "PreState/State" and end state
    // fill first row (header)
    sprintf(table[0][0], "#Region");
    sprintf(table[0][1], "State");
    for (int state = 0; state < model->numberOfStates; state++) {
        const char *stateName = EmissionDistSeries_getStateName(state);
        sprintf(table[0][2 + state], "%s", stateName);
    }
    sprintf(table[0][2 + model->numberOfStates], "End");

    // write a transition matrix for each region one at a time
    int numberOfRows = 1;
    for (int region = 0; region < model->numberOfRegions; region++) {
        Transition *transition = model->transitionPerRegion[region];
        // fill the first and second columns for this region
        for (int state = 0; state < model->numberOfStates; state++) {
            sprintf(table[numberOfRows + state][0], "%d", region);
            const char *stateName = EmissionDistSeries_getStateName(state);
            sprintf(table[numberOfRows + state][1], "%s", stateName);
        }
        sprintf(table[numberOfRows + model->numberOfStates][0], "%d", region);
        sprintf(table[numberOfRows + model->numberOfStates][1], "Start");
        // fill the probability matrix
        for (int preState = 0; preState < model->numberOfStates + 1; preState++) {
            for (int state = 0; state < model->numberOfStates + 1; state++) {
                sprintf(table[numberOfRows + preState][2 + state],
                        "%.2e", transition->matrix->data[preState][state]);
            }
        }
        numberOfRows += model->numberOfStates + 1; // +1 for "Start" state
    }

    for (int row = 0; row < numberOfRows; row++) {
        for (int col = 0; col < numberOfColumns - 1; col++) {
            fprintf(fout, "%s\t", table[row][col]);
        }
        // last column needs \n instead of \t
        fprintf(fout, "%s\n", table[row][numberOfColumns - 1]);
    }
}

void HMM_printEmissionParametersInTsvFormat(HMM *model, FILE *fout) {
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
    for (int region = 0; region < model->numberOfRegions; region++) {
        sprintf(table[0][4 + region], "Values_Region_%d", region);
    }
    int numberOfRows = 1;
    EmissionDistSeries *emissionDistSeries = model->emissionDistSeriesPerRegion[0];
    for (int state = 0; state < emissionDistSeries->numberOfDists; state++) {
        int numberOfParams;
        const char *stateName = EmissionDistSeries_getStateName(state);
        const char *distributionName = EmissionDistSeries_getDistributionName(emissionDistSeries, state);
        const char **parameterNames = EmissionDistSeries_getParameterNames(emissionDistSeries, state, &numberOfParams);
        int numberOfComps = EmissionDistSeries_getNumberOfComps(emissionDistSeries, state);
        for (int param = 0; param < numberOfParams; param++) {
            sprintf(table[numberOfRows][0], "%s", stateName);
            sprintf(table[numberOfRows][1], "%s", distributionName);
            sprintf(table[numberOfRows][2], "%d", numberOfComps);
            sprintf(table[numberOfRows][3], "%s", parameterNames[param]);
            numberOfRows += 1;
        }
    }
    for (int region = 0; region < model->numberOfRegions; region++) {
        numberOfRows = 1;
        EmissionDistSeries *emissionDistSeries = model->emissionDistSeriesPerRegion[region];
        for (int state = 0; state < emissionDistSeries->numberOfDists; state++) {
            int numberOfParams;
            double **parameterValues = EmissionDistSeries_getParameterValues(emissionDistSeries, state,
                                                                             &numberOfParams);
            int numberOfComps = EmissionDistSeries_getNumberOfComps(emissionDistSeries, state);
            for (int param = 0; param < numberOfParams; param++) {
                char *parameterValuesStr = String_joinDoubleArray(parameterValues[param], numberOfComps, ',');
                sprintf(table[numberOfRows][4 + region], "%s", parameterValuesStr);
                free(parameterValuesStr);
                numberOfRows += 1;
            }
            free(parameterValues);
        }
    }


    for (int row = 0; row < numberOfRows; row++) {
        for (int col = 0; col < numberOfColumns - 1; col++) {
            fprintf(fout, "%s\t", table[row][col]);
        }
        // last column needs \n instead of \t
        fprintf(fout, "%s\n", table[row][numberOfColumns - 1]);
    }
}

void HMM_destruct(HMM *model) {
    for (int region = 0; region < model->numberOfRegions; region++) {
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
    em->emissionDistSeriesPerRegion = EmissionDistSeries_copy1DArray(model->emissionDistSeriesPerRegion, model->numberOfRegions);
    em->transitionPerRegion = Transition_copy1DArray(model->transitionPerRegion, model->numberOfRegions);
    // Initialize scale to avoid underflow
    em->scales = Double_construct1DArray(em->seqLen);
    em->px = -1.0;
    em->numberOfRegions = model->numberOfRegions;
    return em;
}

void EM_destruct(EM *em) {
    Double_destruct2DArray(em->f, em->seqLen);
    Double_destruct2DArray(em->b, em->seqLen);
    Double_destruct1DArray(em->scales);
    EmissionDistSeries_destruct1DArray(em->emissionDistSeriesPerRegion, em->numberOfRegions);
    Transition_destruct1DArray(em->transitionPerRegion, em->numberOfRegions);
}

void EM_renewParametersAndEstimatorsFromModel(EM *em, HMM *model){
    if(em->emissionDistSeriesPerRegion != NULL) {
        EmissionDistSeries_destruct1DArray(em->emissionDistSeriesPerRegion, em->numberOfRegions);
    }
    if(em->transitionPerRegion != NULL) {
        Transition_destruct1DArray(em->transitionPerRegion, em->numberOfRegions);
    }
    em->emissionDistSeriesPerRegion = EmissionDistSeries_copy1DArray(model->emissionDistSeriesPerRegion, model->numberOfRegions);
    em->transitionPerRegion = Transition_copy1DArray(model->transitionPerRegion, model->numberOfRegions);
}


///////////////////////////////////////
// Functions for forward algorithm   //
//////////////////////////////////////

void EM_resetAllColumnsForward(EM *em) {
    HMM *model = em->model;
    // Initialize to zero
    for (int i = 0; i < em->seqLen; i++) {
        for (int s = 0; s < model->numberOfStates; s++) {
            em->f[i][s] = 0.0;
        }
    }
}

void EM_fillFirstColumnForward(EM *em) {
    HMM *model = em->model;
    double scale = 0.0; // For scaling the forward probabilities at each location to avoid underflow
    double eProb;
    double tProb;
    uint8_t preX = 0;
    double alpha = 0.0; // for the first column alpha can not be greater than 0
    // Set the 0-th block of the forward matrix
    scale = 0.0;
    for (int state = 0; state < model->numberOfStates; state++) {
        uint8_t region = CoverageInfo_getRegionIndex(em->coverageInfoSeq[0]);
        uint8_t x = em->coverageInfoSeq[0]->coverage;
        // Emission probability
        eProb = EmissionDistSeries_getProb(model->emissionDistSeriesPerRegion[region],
                                           state,
                                           x,
                                           preX,
                                           alpha);
        // Transition probability
        tProb = Transition_getStartProb(model->transitionPerRegion[region], state);
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

void EM_fillOneColumnForward(EM *em, int columnIndex) {
    if (columnIndex == 0) {
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
    double scale = 0.0;
    region = CoverageInfo_getRegionIndex(em->coverageInfoSeq[i]);
    preRegion = CoverageInfo_getRegionIndex(em->coverageInfoSeq[i - 1]);
    x = em->coverageInfoSeq[i]->coverage;
    preX = em->coverageInfoSeq[i - 1]->coverage;
    for (int state = 0; state < model->numberOfStates; state++) {
        for (int preState = 0; preState < model->numberOfStates; preState++) { // Transition from c1 comp to c2 comp
            alpha = model->alpha->data[preState][state];
            // Emission probability
            // Not that alpha can be zero and in that case emission probability is not
            // dependent on the previous observation
            eProb = EmissionDistSeries_getProb(model->emissionDistSeriesPerRegion[region],
                                               state,
                                               x,
                                               preX,
                                               model->alpha->data[preState][state]);
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
        scale += em->f[i][state];
    }
    em->scales[i] = scale;
    if (em->scales[i] < 1e-50) {
        fprintf(stderr, "scale (= %.2e) is very low!\n", em->scales[i]);
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


void EM_resetAllColumnsBackward(EM *em) {
    HMM *model = em->model;
    // Initialize to zero
    for (int i = 0; i < em->seqLen; i++) {
        for (int s = 0; s < model->numberOfStates; s++) {
            em->b[i][s] = 0.0;
        }
    }
}


void EM_fillLastColumnBackward(EM *em) {
    HMM *model = em->model;
    double tProb;
    // Set the last block of the backward matrix
    for (int state = 0; state < model->numberOfStates; state++) {
        uint8_t region = CoverageInfo_getRegionIndex(em->coverageInfoSeq[em->seqLen - 1]);
        // Transition probability
        tProb = Transition_getTerminationProb(model->transitionPerRegion[region], state);
        // Update backward
        em->b[em->seqLen - 1][state] = tProb;
    }
    // Scale last column of the backward matrix
    for (int state = 0; state < model->numberOfStates; state++) {
        em->b[em->seqLen - 1][state] /= em->scales[em->seqLen - 1];
    }
}


void EM_fillOneColumnBackward(EM *em, int columnIndex) {
    if (columnIndex == em->seqLen - 1) {
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
    region = CoverageInfo_getRegionIndex(em->coverageInfoSeq[i + 1]);
    preRegion = CoverageInfo_getRegionIndex(em->coverageInfoSeq[i]);
    covInfo = em->coverageInfoSeq[i + 1];
    preCovInfo = em->coverageInfoSeq[i];
    x = covInfo->coverage;
    preX = preCovInfo->coverage;
    for (int state = 0; state < model->numberOfStates; state++) {
        for (int preState = 0; preState < model->numberOfStates; preState++) { // Transition from c1 comp to c2 comp
            alpha = model->alpha->data[preState][state];
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
            //if (em->b[i][preState] == 0.0){
            //fprintf(stdout, "i=%d, tProb =%.2e, eProb =%.2e,em->b[i + 1][state]=%.2e\n", i, tProb , eProb , em->b[i + 1][state]);
            //	fprintf(stdout, "x = %d, state = %d, preState=%d, eProb=%.2e\n", x, state, preState, eProb);
            //}
        }
    }
    if (em->scales[i] < 1e-50) {
        fprintf(stderr, "scale (= %.2e) is very low!\n", em->scales[i]);
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
    for (int columnIndex = em->seqLen - 1; columnIndex >= 0; columnIndex--) {
        EM_fillOneColumnBackward(em, columnIndex);
    }
    // Update P(x)
    // TODO: P(x) is not calculated here p(x) = s(1) x s(2) x ... x s(L) check Durbin's
    // It shouldn't make any problem in the training process since px is not needed in the scaled calculations
    em->px = 0;
}


void EM_updateModelEstimators(EM *em) {
    HMM *model = em->model;
    for(int region=0; region < em->numberOfRegions; region++){
        EmissionDistSeries *emissionDistSeriesForEM = em->emissionDistSeriesPerRegion[region];
        EmissionDistSeries *emissionDistSeriesForModel = model->emissionDistSeriesPerRegion[region];
        EmissionDistSeries_updateEstimatorFromOtherEstimator(emissionDistSeriesForModel,
                                                             emissionDistSeriesForEM);
        TransitionCountData *transitionDataForEM = em->transitionPerRegion[region]->transitionCountData;
        TransitionCountData *transitionDataForModel = model->transitionPerRegion[region]->transitionCountData;
        TransitionCountData_incrementFromOtherCountData(transitionDataForModel,
                                                        transitionDataForEM);
    }
}


void EM_updateEstimatorsUsingOneColumn(EM *em, int columnIndex) {
    if (columnIndex == em->seqLen - 1) {
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

    // get observations
    region = CoverageInfo_getRegionIndex(em->coverageInfoSeq[i + 1]);
    preRegion = CoverageInfo_getRegionIndex(em->coverageInfoSeq[i]);
    x = em->coverageInfoSeq[i + 1]->coverage;
    preX = em->coverageInfoSeq[i]->coverage;

    emissionDistSeries = em->emissionDistSeriesPerRegion[region];
    transition = em->transitionPerRegion[region];
    for (int state = 0; state < model->numberOfStates; state++) {
        for (int preState = 0; preState < model->numberOfStates; preState++) { // Transition from c1 comp to c2 comp
            // get model attributes
            alpha = model->alpha->data[preState][state];
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
            double adjustedCount = count / transition->terminationProb;
            EmissionDistSeries_updateEstimator(emissionDistSeries,
                                               state,
                                               x,
                                               preX,
                                               alpha,
                                               adjustedCount);

            /*if(count < 1e-3 && i % 50000  == 0){
            fprintf(stdout, "i=%d, count=%.2e, em->f[%d][%d]=%.2e, tProb=%.2e, eProb=%.2e, em->b[%d][%d]=%.2e, scale=%.2e\n", i, count, i, preState, em->f[i][preState],tProb,eProb, i+1, state, em->b[i + 1][state], em->scales[i]);
            }*/
            TransitionCountData_increment(transition->transitionCountData,
                                          adjustedCount,
                                          preState,
                                          state);
        }
    }
}

void EM_updateEstimators(EM *em) {
    // skip first column since alpha might be > 0
    for (int columnIndex = 1; columnIndex < em->seqLen; columnIndex++) {
        EM_updateEstimatorsUsingOneColumn(em, columnIndex);
    }
}

bool EM_estimateParameters(EM *em, double convergenceTol) {
    bool converged = true;
    HMM *model = em->model;
    for (int region = 0; region < model->numberOfRegions; region++) {
        converged &= EmissionDistSeries_estimateParameters(model->emissionDistSeriesPerRegion[region], convergenceTol);
        converged &= Transition_estimateTransitionMatrix(model->transitionPerRegion[region], convergenceTol);
    }
    return converged;
}

void EM_resetEstimators(EM *em) {
    HMM *model = em->model;
    for (int region = 0; region < model->numberOfRegions; region++) {
        EmissionDistSeries_resetParameterEstimators(model->emissionDistSeriesPerRegion[region]);
        Transition_resetCountData(model->transitionPerRegion[region]);
    }
}


double *EM_getPosterior(EM *em, int pos) {
    HMM *model = em->model;
    double *posterior = malloc(model->numberOfStates * sizeof(double));
    double total = 0.0;
    for (int s = 0; s < model->numberOfStates; s++) {
        posterior[s] = em->f[pos][s] * em->b[pos][s] * em->scales[pos];
        total += posterior[s];
        //fprintf(stdout, "s=%d, %.2e\t",s, posterior[s]);
    }
    for (int s = 0; s < model->numberOfStates; s++) {
        posterior[s] /= total;
    }
    //fprintf(stdout, "\n");
    return posterior;
}

int EM_getMostProbableState(EM *em, int pos) {
    double *posterior = EM_getPosterior(em, pos);
    int state = Double_getArgMaxIndex1DArray(posterior, em->model->numberOfStates);
    free(posterior);
    return state;
}

void EM_printPosteriorInTsvFormat(EM *em, FILE *fout) {
    HMM *model = em->model;
    EmissionDistSeries *emissionDistSeries = model->emissionDistSeriesPerRegion[0];
    fprintf(fout, "position\t");
    for (int state = 0; state < model->numberOfStates; state++) {
        const char *stateName = EmissionDistSeries_getStateName(state);
        fprintf(fout, "%s_%d_posterior\t", stateName, state);
    }
    fprintf(fout, "state_index_prediction\n");
    for (int pos = 0; pos < em->seqLen; pos++) {
        fprintf(fout, "%d\t", pos);
        double *posterior = EM_getPosterior(em, pos);
        int prediction = EM_getMostProbableState(em, pos);
        for (int state = 0; state < model->numberOfStates; state++) {
            fprintf(fout, "%.3f\t", posterior[state]);
        }
        fprintf(fout, "%d\n", prediction);
        free(posterior);
    }
}


void EM_runOneIterationAndUpdateEstimatorsForThreadPool(void *arg_) {
    work_arg_t *arg = arg_;
    EM *em = arg->data;

    EM_runForward(em);
    EM_runBackward(em);
    EM_updateEstimators(em);
    EM_updateModelEstimators(em);

    // update prediction labels
    for (int pos = 0; pos < em->seqLen; pos++) {
        CoverageInfo *coverageInfo = em->coverageInfoSeq[pos];
        if (coverageInfo->data != NULL) {
            Inference *inference = coverageInfo->data;
            inference->prediction = EM_getMostProbableState(em, pos);
        }
    }
}

void EM_runOneIterationForList(stList *emList, HMM *model, int threads) {
    tpool_t *tm = tpool_create(threads);
    for (int i = 0; i < stList_length(emList); i++) {

        // get EM struct for this chunk index
        EM *em = stList_get(emList, i);
        EM_renewParametersAndEstimatorsFromModel(em, model);
        // add EM to the work struct
        work_arg_t *argWork = malloc(sizeof(work_arg_t));
        argWork->data = (void *) em;
        // queue job
        tpool_add_work(tm, EM_runOneIterationAndUpdateEstimatorsForThreadPool, argWork);
    } // jobs for running EM on all chunks are queued

    // wait until all jobs are finished for this iteration
    tpool_wait(tm);
    // destroy thread pool
    tpool_destroy(tm);
}
