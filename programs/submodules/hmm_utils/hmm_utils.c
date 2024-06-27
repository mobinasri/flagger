//
// Created by mobin on 1/8/24.
//

#include "hmm_utils.h"
#include <string.h>
#include <stdint.h>
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "common.h"
#include "ptBlock.h"
#include "stdbool.h"


ParameterBinding *ParameterBinding_construct(double **coefs, int numberOfParams, int numberOfComps){
    ParameterBinding *parameterBinding = malloc(1 * sizeof(ParameterBinding));
    parameterBinding->coefs = Double_copy2DArray(coefs, numberOfParams, numberOfComps);
    parameterBinding->numberOfComps = numberOfParams;
    parameterBinding->numberOfParams = numberOfParams;
    return parameterBinding;
}

ParameterBinding *ParameterBinding_constructSequenceByStep(double* firstCoefsPerParam,
                                                           double* stepsPerParam,
                                                           int numberOfParams,
                                                           int numberOfComps){
    ParameterBinding *parameterBinding = malloc(1 * sizeof(ParameterBinding));
    parameterBinding->coefs = Double_construct2DArray(numberOfParams, numberOfComps);
    for(int p = 0; p < numberOfParams; p++) {
        for (int c = 0; c < numberOfComps; c++) {
            parameterBinding->coefs[p][c] = firstCoefsPerParam[p] + (double) c * stepsPerParam[p];
        }
    }
    return parameterBinding;
}

ParameterBinding *ParameterBinding_constructForSingleComp(double* coefsPerParam, int numberOfParams){
    double *stepsPerParam = Double_construct1DArray(numberOfParams);
    Double_fill1DArray(stepsPerParam,  numberOfParams, 0.0);
    ParameterBinding *parameterBinding =  ParameterBinding_constructSequenceByStep(coefsPerParam,
                                                                                   stepsPerParam,
                                                                                   numberOfParams,
                                                                                   1);
    Double_destruct1DArray(stepsPerParam);
    return parameterBinding;
}

double ParameterBinding_getValue(ParameterBinding* parameterBinding, int paramIndex, int compIndex){
    return parameterBinding->coefs[paramIndex][compIndex];
}

void ParameterBinding_destruct(ParameterBinding* parameterBinding){
    Double_destruct2DArray(parameterBinding->coefs, parameterBinding->numberOfParams);
    free(parameterBinding);
}


ParameterBinding **ParameterBinding_getDefault1DArrayForModel(ModelType modelType, int numberOfCollapsedComps){
    if(modelType == MODEL_GAUSSIAN){
        return ParameterBinding_getDefault1DArrayForGaussian(numberOfCollapsedComps);
    }else if (modelType == MODEL_TRUNC_EXP_GAUSSIAN){
        return ParameterBinding_getDefault1DArrayForTruncExpGaussian(numberOfCollapsedComps);
    }else if (modelType == MODEL_NEGATIVE_BINOMIAL){
        return ParameterBinding_getDefault1DArrayForNegativeBinomial(numberOfCollapsedComps);
    }
}

ParameterBinding **ParameterBinding_getDefault1DArrayForGaussian(int numberOfCollapsedComps){
    int numberOfParams = 2; // mean and var
    ParameterBinding ** parameterBinding1DArray = malloc(NUMBER_OF_STATES * sizeof(ParameterBinding*));

    // create one array and fill them for each state separately
    double *coefsPerParam = Double_construct1DArray(numberOfParams);

    coefsPerParam[GAUSSIAN_MEAN] = ERR_COMP_BINDING_COEF;
    coefsPerParam[GAUSSIAN_VAR] = ERR_COMP_BINDING_COEF;
    parameterBinding1DArray[STATE_ERR] = ParameterBinding_constructForSingleComp(coefsPerParam, numberOfParams);

    coefsPerParam[GAUSSIAN_MEAN] = 0.5;
    coefsPerParam[GAUSSIAN_VAR] = 0.5;
    parameterBinding1DArray[STATE_DUP] = ParameterBinding_constructForSingleComp(coefsPerParam, numberOfParams);

    coefsPerParam[GAUSSIAN_MEAN] = 1.0;
    coefsPerParam[GAUSSIAN_VAR] = 1.0;
    parameterBinding1DArray[STATE_HAP] = ParameterBinding_constructForSingleComp(coefsPerParam, numberOfParams);

    coefsPerParam[GAUSSIAN_MEAN] = 1.0;
    coefsPerParam[GAUSSIAN_VAR] = 1.0;
    parameterBinding1DArray[STATE_MSJ] = ParameterBinding_constructForSingleComp(coefsPerParam, numberOfParams);

    double *stepsPerParam = Double_construct1DArray(numberOfParams);
    coefsPerParam[GAUSSIAN_MEAN] = 2.0;
    coefsPerParam[GAUSSIAN_VAR] = 2.0;
    stepsPerParam[GAUSSIAN_MEAN] = 1.0;
    stepsPerParam[GAUSSIAN_VAR] = 1.0;
    parameterBinding1DArray[STATE_COL] = ParameterBinding_constructSequenceByStep(coefsPerParam,
                                                                                  stepsPerParam,
                                                                                  numberOfParams,
                                                                                  numberOfCollapsedComps);
    Double_destruct1DArray(coefsPerParam);
    Double_destruct1DArray(stepsPerParam);
    return parameterBinding1DArray;
}

ParameterBinding **ParameterBinding_getDefault1DArrayForNegativeBinomial(int numberOfCollapsedComps){
    int numberOfParams = 2; // mean and var
    ParameterBinding **parameterBinding1DArray = malloc(NUMBER_OF_STATES * sizeof(ParameterBinding*));

    // create one array and fill them for each state separately
    double *coefsPerParam = Double_construct1DArray(numberOfParams);

    coefsPerParam[NB_MEAN] = ERR_COMP_BINDING_COEF;
    coefsPerParam[NB_VAR] = ERR_COMP_BINDING_COEF;
    parameterBinding1DArray[STATE_ERR] = ParameterBinding_constructForSingleComp(coefsPerParam, numberOfParams);

    coefsPerParam[NB_MEAN] = 0.5;
    coefsPerParam[NB_VAR] = 0.5;
    parameterBinding1DArray[STATE_DUP] = ParameterBinding_constructForSingleComp(coefsPerParam, numberOfParams);

    coefsPerParam[NB_MEAN] = 1.0;
    coefsPerParam[NB_VAR] = 1.0;
    parameterBinding1DArray[STATE_HAP] = ParameterBinding_constructForSingleComp(coefsPerParam, numberOfParams);

    coefsPerParam[NB_MEAN] = 1.0;
    coefsPerParam[NB_VAR] = 1.0;
    parameterBinding1DArray[STATE_MSJ] = ParameterBinding_constructForSingleComp(coefsPerParam, numberOfParams);

    double *stepsPerParam = Double_construct1DArray(numberOfParams);
    coefsPerParam[NB_MEAN] = 2.0;
    coefsPerParam[NB_VAR] = 2.0;
    stepsPerParam[NB_MEAN] = 1.0;
    stepsPerParam[NB_VAR] = 1.0;
    parameterBinding1DArray[STATE_COL] = ParameterBinding_constructSequenceByStep(coefsPerParam,
                                                                                  stepsPerParam,
                                                                                  numberOfParams,
                                                                                  numberOfCollapsedComps);
    Double_destruct1DArray(coefsPerParam);
    Double_destruct1DArray(stepsPerParam);
    return parameterBinding1DArray;
}

ParameterBinding **ParameterBinding_getDefault1DArrayForTruncExpGaussian(int numberOfCollapsedComps){
    ParameterBinding **parameterBinding1DArray = ParameterBinding_getDefault1DArrayForGaussian(numberOfCollapsedComps);
    ParameterBinding_destruct(parameterBinding1DArray[STATE_ERR]);

    int numberOfParams = 1;
    double *coefsPerParam = Double_construct1DArray(numberOfParams);

    coefsPerParam[TRUNC_EXP_MEAN] = 0.0; // zero means no binding
    parameterBinding1DArray[STATE_ERR] = ParameterBinding_constructForSingleComp(coefsPerParam, numberOfParams);
    Double_destruct1DArray(coefsPerParam);

    return parameterBinding1DArray;
}


void ParameterBinding_destruct1DArray(ParameterBinding **parameterBinding1DArray, int arrayLength){
    for(int s=0; s < arrayLength; s++) {
        ParameterBinding_destruct(parameterBinding1DArray[s]);
    }
    free(parameterBinding1DArray);
}


CountData *CountData_construct(uint8_t length){
    CountData *countData = malloc(sizeof(CountData));
    countData->counts = malloc(length * sizeof(double));
    memset(countData->counts, 0, length * sizeof(double));
    countData->countsLength = length;
    countData->totalCount = 0;
    countData->sum = 0.0;
    return countData;
}

CountData **CountData_construct1DArray(uint8_t length, int arrayLength){
    CountData **countData1DArray = malloc(arrayLength * sizeof(CountData*));
    for(int i=0; i < arrayLength; i++){
        countData1DArray[i] = CountData_construct(length);
    }
    return countData1DArray;
}

void CountData_destruct1DArray(CountData **countData1DArray, int arrayLength){
    for(int i=0; i < arrayLength; i++){
        CountData_destruct(countData1DArray[i]);
    }
    free(countData1DArray);
}

void CountData_increment(CountData *countData, uint8_t value, double count) {
    uint8_t effValue = value;
    if (countData->countsLength <= value) {
        effValue = countData->countsLength - 1;
    }
    countData->counts[effValue] += count;
    countData->sum += effValue * count;
    countData->totalCount += count;
}

void CountData_reset(CountData* countData){
    memset(countData->counts, 0, countData->countsLength * sizeof(double));
    countData->sum = 0.0;
    countData->totalCount = 0;
}

int CountData_getMostFrequentValue(CountData *countData, int minValue, int maxValue){
    int value = -1;
    double maxCount = -1.0;
    for(int i=minValue; i < min(maxValue, countData->countsLength); i++){
        if(maxCount <= countData->counts[i]){
            value = i;
            maxCount = countData->counts[i];
        }
    }
    return value;
}

int CountData_getMaxCount(CountData *countData, int minValue, int maxValue){
    int value = -1;
    double maxCount = -1.0;
    for(int i=minValue; i < min(maxValue, countData->countsLength); i++){
        if(maxCount <= countData->counts[i]){
            value = i;
            maxCount = countData->counts[i];
        }
    }
    return maxCount;
}

void CountData_destruct(CountData* countData){
    free(countData->counts);
    free(countData);
}

/**
 * Constructs a Negative Binomial object containing the parameters for all of its mixture components
 *
 * @param mean	Initial mean values
 * @param var	Initial variance values
 * @param numberOfComps	The number of mixture components
 * @return nb
 */
NegativeBinomial *NegativeBinomial_construct(double *mean, double *var, int numberOfComps) {
    NegativeBinomial *nb = malloc(1 * sizeof(NegativeBinomial));
    memcpy(nb->mean, mean, numberOfComps * sizeof(double));
    memcpy(nb->var, var, numberOfComps * sizeof(double));
    // Allocating and initializing mixture weights
    nb->weights = malloc(numberOfComps * sizeof(double));
    nb->numberOfComps = numberOfComps;
    Double_fill1DArray(nb->weights, numberOfComps, 1.0 / numberOfComps);
    return nb;
}


/**
 * Destructs a Negative Binomial object
 *
 * @param nb
 */
void NegativeBinomial_destruct(NegativeBinomial *nb) {
    free(nb->mean);
    free(nb->var);
    free(nb->weights);
    free(nb);
}

/**
 * Constructs a NB object in which variance values
 * are initialized to the initial mean values times a constant factor greater than 1
 *
 * @param mean	                    mean values
 * @param numberOfMixtures	        The number of mixture components
 * @param factor	                A constant factor greater than 1 for computing variance values
 * @return nb
 */
NegativeBinomial *NegativeBinomial_constructByMean(double *mean, double factor, int numberOfComps){
    double *var = malloc(numberOfComps * sizeof(double));
    for (int m = 0; m < numberOfComps; m++) {
        var[m] = mean[m] * factor;
    }
    NegativeBinomial *nb = NegativeBinomial_construct(mean, var, numberOfComps);
    free(var);
    return nb;
}


double NegativeBinomial_getTheta(double mean, double var) {
    double theta = mean / var;
    return theta;
}

double NegativeBinomial_getR(double mean, double var) {
    double r = pow(mean, 2) / (var - mean);
    return r;
}


double NegativeBinomial_getMean(double theta, double r) {
    double mean = r * (1 - theta) / theta;
    return mean;
}

double NegativeBinomial_getVar(double theta, double r) {
    double var = r * (1 - theta) / pow(theta, 2);
    return var;
}

/**
 * Returns the total probability of emitting a value from the given Negative Binomial
 *
 * @param x	    The emitted value
 * @param nb	The NegativeBinomial object
 * @return prob	The emission probability
 */
double NegativeBinomial_getProb(NegativeBinomial *nb, uint8_t x) {
    double *compProbs = NegativeBinomial_getComponentProbs(nb, x);
    double totProb = Double_sum1DArray(compProbs, nb->numberOfComps);
    free(compProbs);
    return totProb;
}


/**
 * Returns the probabilities of emitting a value from the mixture components of Negative Binomial
 *
 * @param x         The emitted value
 * @param nb        The Negative Binomial object
 * @return probs     An array of probabilities for all NB components
 */
double *NegativeBinomial_getComponentProbs(NegativeBinomial *nb, uint8_t x) {
    double w;
    double theta;
    double r;
    double *probs = malloc(nb->numberOfComps * sizeof(double));
    // iterate over mixture components
    for (int m = 0; m < nb->numberOfComps; m++) {
        theta = NegativeBinomial_getTheta(nb->mean[m], nb->var[m]);
        r = NegativeBinomial_getR(nb->mean[m], nb->var[m]);
        w = nb->weights[m];
        probs[m] = w * exp(lgamma(r + x) - lgamma(r) - lgamma(x + 1) + r * log(theta) +
                           (double) x * log(1 - theta));
        if (probs[m] != probs[m]) {
            fprintf(stderr, "prob is NAN\n");
            exit(EXIT_FAILURE);
        }
        if (probs[m] < 1e-40) {
            probs[m] = 1e-40;
        }
    }
    return probs;
}




/**
 * Constructs a Gaussian object containing the parameters for all of its mixture components
 *
 * @param mean	Initial mean values
 * @param var	Initial variance values
 * @param numberOfComps	The number of mixture components
 * @return gaussian
 */
Gaussian *Gaussian_construct(double *mean, double *var, int numberOfComps) {
    Gaussian *gaussian = malloc(1 * sizeof(Gaussian));
    memcpy(gaussian->mean, mean, numberOfComps * sizeof(double));
    memcpy(gaussian->var, var, numberOfComps * sizeof(double));
    // Allocating and initializing mixture weights
    gaussian->weights = malloc(numberOfComps * sizeof(double));
    gaussian->numberOfComps = numberOfComps;
    Double_fill1DArray(gaussian->weights, numberOfComps, 1.0 / numberOfComps);
    return gaussian;
}

/**
 * Destructs a Gaussian object
 *
 * @param gaussian
 */
void Gaussian_destruct(Gaussian *gaussian) {
    free(gaussian->mean);
    free(gaussian->var);
    free(gaussian->weights);
    free(gaussian);
}

/**
 * Constructs a Gaussian object in which variance values are initialized
 * to the initial mean values times a constant factor greater than 1
 *
 * @param mean	                    Initial mean values
 * @param numberOfMixtures	        The number of mixture components
 * @param factor	                A constant factor greater than 1 for computing variance values
 * @return gaussian
 */
Gaussian *Gaussian_constructByMean(double *mean, double factor, int numberOfComps){
    double *var = malloc(numberOfComps * sizeof(double));
    for (int m = 0; m < numberOfComps; m++) {
        var[m] = mean[m] * factor;
    }
    Gaussian *gaussian = Gaussian_construct(mean, var, numberOfComps);
    free(var);
    return gaussian;
}

/**
 * Returns the total probability of emitting a value from the given Gaussian
 *
 * @param x	        The emitted value
 * @param gaussian	The gaussian object
 * @return prob	    The emission probability
 */

double Gaussian_getProb(Gaussian *gaussian, uint8_t x, uint8_t preX, double alpha) {
    double *probs = Gaussian_getComponentProbs(gaussian, x, preX, alpha);
    double totProb = Double_sum1DArray(probs, gaussian->numberOfComps);
    free(probs);
    return totProb;
}


/**
 * Returns the probabilities of emitting a value from the mixture components of the given Gaussian
 *
 * @param x             The emitted value
 * @param gaussian      The Gaussian object
 * @return probs        An array of probabilities for all Gaussian components
 */
double *Gaussian_getComponentProbs(Gaussian *gaussian, uint8_t x, uint8_t preX, double alpha) {
    double mean;
    double var;
    double w;
    double *probs = malloc(gaussian->numberOfComps * sizeof(double));
    // iterate over mixture components
    for (int m = 0; m < gaussian->numberOfComps; m++) {
        mean = (1 - alpha) * gaussian->mean[m] + alpha * preX;
        var = gaussian->var[m];
        w = gaussian->weights[m];
        // adjust the mean value based on the previous observation and alpha (dependency factor)
        probs[m] = w / (sqrt(var * 2 * PI)) * exp(-0.5 * pow((x - mean), 2) / var);
        if (probs[m] != probs[m]) {
            double u = -0.5 * pow((x - mean), 2) / var;
            fprintf(stderr, "[Error] prob is NAN exp(%.2e) mean=%.2e, var=%.2e\n", u, mean, var);
            exit(EXIT_FAILURE);
        }
        if (probs[m] < 1e-40) {
            fprintf(stderr, "[Warning] prob is lower than 1e-40. It is set to 1e-40\n");
            probs[m] = 1e-40;
        }
    }
    return probs;
}



TruncExponential *TruncExponential_construct(double lambda, double truncPoint){
    TruncExponential *truncExponential = malloc(sizeof(TruncExponential));
    truncExponential->lambda = lambda;
    truncExponential->truncPoint = truncPoint;
    return truncExponential;
}

void TruncExponential_destruct(TruncExponential* truncExponential){
    free(truncExponential);
}

double TruncExponential_getProb(TruncExponential* truncExponential, uint8_t x) {
    double lam = truncExponential->lambda;
    double b = truncExponential->truncPoint;
    return lam * exp(-lam * x) / (1 - exp(-lam * b));
}

double TruncExponential_getLogLikelihoodByParams(double lambda, double truncPoint, CountData *countData) {
    double num = countData->sum;
    double denom = countData->totalCount;
    double lam = lambda;
    double b = truncPoint;
    return denom * log(lam) - denom * log(1.0 - exp(-lam * b)) - num * lam;
}

double TruncExponential_getLogLikelihood(TruncExponential *truncExponential, CountData *countData) {
    double num = countData->sum;
    double denom = countData->totalCount;
    double lam = truncExponential->lambda;
    double b = truncExponential->truncPoint;
    return denom * log(lam) - denom * log(1.0 - exp(-lam * b)) - num * lam;
}

// taken from Jordan's script
// which is originally adapted from https://en.wikipedia.org/wiki/Golden-section_search
// tol = 1e-6
double TruncExponential_estimateLambda(TruncExponential *truncExponential, CountData *countData, double tol) {
    double a = 0.0;
    double b = truncExponential->truncPoint;

    double invphi = (sqrt(5.0) - 1.0) / 2.0;  // 1 / phi
    double invphi2 = (3.0 - sqrt(5.0)) / 2.0;  // 1 / phi^2

    double h = b - a;
    if (h <= tol) {
        return (b + a) / 2.0;
    }

    int n = ceil(log(tol / h) / log(invphi));

    double c = a + invphi2 * h;
    double d = a + invphi * h;
    double yc = TruncExponential_getLogLikelihoodByParams(c, truncExponential->truncPoint, countData); // f(c)
    double yd = TruncExponential_getLogLikelihoodByParams(d, truncExponential->truncPoint, countData); // f(d)
    for (int k = 0; k < n - 1; k++) {
        if (yc > yd) {
            b = d;
            d = c;
            yd = yc;
            h = invphi * h;
            c = a + invphi2 * h;
            yc = TruncExponential_getLogLikelihoodByParams(c, truncExponential->truncPoint, countData); // f(c)
        } else {
            a = c;
            c = d;
            yc = yd;
            h = invphi * h;
            d = a + invphi * h;
            yd = TruncExponential_getLogLikelihoodByParams(d, truncExponential->truncPoint, countData); //f(d)
        }
    }

    if (yc > yd) {
        return (a + d) / 2.0;
    } else {
        return (c + b) / 2.0;
    }
}


EmissionDist *EmissionDist_construct(void* dist, DistType distType){
    EmissionDist *emissionDist = malloc(1 * sizeof(EmissionDist));
    emissionDist->dist = dist;
    emissionDist->distType = distType;
    return emissionDist;
}


void EmissionDist_destruct(EmissionDist* emissionDist){
    if (emissionDist->distType == DIST_TRUNC_EXPONENTIAL){
        TruncExponential_destruct((TruncExponential *) emissionDist->dist);
    }
    else if (emissionDist->distType == DIST_GAUSSIAN){
        Gaussian_destruct((Gaussian *) emissionDist->dist);
    }
    else if (emissionDist->distType == DIST_NEGATIVE_BINOMIAL){
        NegativeBinomial_destruct((NegativeBinomial *) emissionDist->dist);
    }
    free(emissionDist);
}

double EmissionDist_getProb(EmissionDist *emissionDist, uint8_t x, uint8_t preX, double alpha){
    if (emissionDist->distType == DIST_TRUNC_EXPONENTIAL){
        return TruncExponential_getProb((TruncExponential *) emissionDist->dist, x);
    }
    else if (emissionDist->distType == DIST_GAUSSIAN){
        return Gaussian_getProb((Gaussian *) emissionDist->dist, x, preX, alpha);
    }
    else if (emissionDist->distType == DIST_NEGATIVE_BINOMIAL){
        return NegativeBinomial_getProb((NegativeBinomial *) emissionDist->dist, x);
    }
}


EmissionDistSeries *EmissionDistSeries_constructForModel(ModelType modelType,
                                                         double **means,  // [numberOfDists] x [maxMixtures]
                                                         int *numberOfCompsPerDist,
                                                         int numberOfDists){
    // Constructing the emission objects and setting their parameters
    // emissionDistSeries[c] is pointing to the emission distribution of the c-th dists
    EmissionDistSeries *emissionDistSeries = (EmissionDistSeries *) malloc(1 * sizeof(EmissionDistSeries));
    emissionDistSeries->emissionDists = (EmissionDist **) malloc(numberOfDists * sizeof(EmissionDist*));
    emissionDistSeries->numberOfDists = numberOfDists;

    void *dist;
    // different model types contain different emission distributions
    if (modelType == MODEL_TRUNC_EXP_GAUSSIAN){
        // For MODEL_TRUNC_EXP_GAUSSIAN  the first component is modeled by a truncated exponential distribution
        dist = TruncExponential_construct(1.0, means[STATE_HAP][0] * EXP_TRUNC_POINT_COV_FRACTION);
        emissionDistSeries->emissionDists[STATE_ERR] = EmissionDist_construct(dist, DIST_TRUNC_EXPONENTIAL);
        // Remaining components are modeled by Gaussian distributions
        // Initialize the emission parameters of each state with the given vector
        for (int s = 1; s < numberOfDists; s++) {
            dist = Gaussian_constructByMean(means[s], 1.0, numberOfCompsPerDist[s]);
            emissionDistSeries->emissionDists[s] = EmissionDist_construct(dist, DIST_GAUSSIAN);
        }
    } else if (modelType == MODEL_GAUSSIAN) {
        // Constructing the emission Gaussian and set their parameters
        // emit[r][c] is pointing to Gaussian of the c-th component of the r-th class
        for (int s = 0; s < numberOfDists; s++) {
            dist = Gaussian_constructByMean(means[s], 1.0,numberOfCompsPerDist[s]);
            emissionDistSeries->emissionDists[s] = EmissionDist_construct(dist, DIST_GAUSSIAN);
        }
    } else if (modelType == MODEL_NEGATIVE_BINOMIAL) {
        for (int s = 0; s < numberOfDists; s++) {
            dist = NegativeBinomial_constructByMean(means[s], 1.0, numberOfCompsPerDist[s]);
            emissionDistSeries->emissionDists[s] = EmissionDist_construct(dist, DIST_NEGATIVE_BINOMIAL);
        }
    }

    emissionDistSeries->parameterBindingPerDist =
            ParameterBinding_getDefault1DArrayForModel(modelType,
                                                       numberOfCompsPerDist[STATE_COL]);
    emissionDistSeries->countDataPerDist = CountData_construct1DArray(MAX_COVERAGE_VALUE, numberOfDists);
    emissionDistSeries->numberOfDists = numberOfDists;
    emissionDistSeries->modelType = modelType;
    return emissionDistSeries;
}

EmissionDist *EmissionDistSeries_getEmissionDist(EmissionDistSeries* emissionDistSeries, int distIndex){
    return emissionDistSeries->emissionDists[distIndex];
}

void EmissionDistSeries_destruct(EmissionDistSeries* emissionDistSeries){
    for(int s = 0; s < emissionDistSeries->numberOfDists; s++) {
        EmissionDist_destruct(emissionDistSeries->emissionDists[s]);
    }
    free(emissionDistSeries->emissionDists);
    CountData_destruct1DArray(emissionDistSeries->countDataPerDist,
                              emissionDistSeries->numberOfDists);
    ParameterBinding_destruct1DArray(emissionDistSeries->parameterBindingPerDist,
                                     NUMBER_OF_STATES);
    free(emissionDistSeries);
}

double *EmissionDistSeries_getAllProbs(EmissionDistSeries* emissionDistSeries,
                                       uint8_t x,
                                       uint8_t preX,
                                       double alpha){
    double *probs = malloc(emissionDistSeries->numberOfDists * sizeof(double));
    for(int s=0; s < emissionDistSeries->numberOfDists; s++){
        probs[s] = EmissionDist_getProb(emissionDistSeries->emissionDists[s], x, preX, alpha);
    }
    return probs;
}

double EmissionDistSeries_getProb(EmissionDistSeries* emissionDistSeries,
                                  int distIndex,
                                  uint8_t x,
                                  uint8_t preX,
                                  double alpha){
    return EmissionDist_getProb(emissionDistSeries->emissionDists[distIndex], x, preX, alpha);
}

TransitionRequirements *TransitionRequirements_construct(double minHighlyClippedRatio, double maxHighMapqRatio){
    TransitionRequirements *transitionRequirements = malloc(sizeof(TransitionRequirements));
    transitionRequirements->minHighlyClippedRatio = minHighlyClippedRatio;
    transitionRequirements->maxHighMapqRatio = maxHighMapqRatio;
    return transitionRequirements;
}

void TransitionRequirements_destruct(TransitionRequirements *transitionRequirements){
    free(transitionRequirements);
}


TransitionCountData *TransitionCountData_construct(int numberOfStates){
    TransitionCountData *transitionCountData = malloc(sizeof(TransitionCountData));
    transitionCountData->countMatrix = MatrixDouble_construct0(numberOfStates + 1, numberOfStates + 1);
    transitionCountData->pseudoCountMatrix = MatrixDouble_construct0(numberOfStates + 1, numberOfStates + 1);
    transitionCountData->numberOfStates = numberOfStates;
}
void TransitionCountData_parsePseudoCountFromFile(TransitionCountData *transitionCountData, char* pathToMatrix, int dim){
    MatrixDouble_destruct(transitionCountData->pseudoCountMatrix);
    transitionCountData->pseudoCountMatrix = MatrixDouble_parseFromFile(pathToMatrix, dim ,dim);
}

TransitionCountData *TransitionCountData_destruct(TransitionCountData *transitionCountData){
    free(transitionCountData->countMatrix);
    free(transitionCountData->pseudoCountMatrix);
    free(transitionCountData);
}




/**
 * Makes a transition matrix with uniform probabilities
 *
 * @param dim	The dimension of the transition matrix supposedly
 * 		should be equal to the number of states + 1
 * @return trans The transition matrix
 */
Transition *Transition_constructUniform(int numberOfStates) {
    int dimension = numberOfStates + 1;
    Transition *transition = malloc(sizeof(Transition));
    transition->matrix = MatrixDouble_construct0(dimension, dimension);
    MatrixDouble_setValue(transition->matrix, 1.0 / dimension);
    transition->transitionCountData = TransitionCountData_construct(numberOfStates);
    transition->numberOfValidityFunctions = 0;
    transition->numberOfStates = numberOfStates;
    return transition;
}

/**
 * Makes a transition matrix with probabilities biased toward
 * not changing states
 *
 * @param dim	The dimension of the transition matrix supposedly
 *              should be equal to the number of states + 1
 * @param diagonalProb The probability of staying in the same state
 * @return transition
 */
Transition *Transition_constructSymmetricBiased(int numberOfStates, double diagonalProb) {
    int dimension = numberOfStates + 1;
    Transition *transition = malloc(sizeof(Transition));
    transition->matrix = MatrixDouble_construct0(dimension, dimension);
    MatrixDouble_setValue(transition->matrix, (1.0 - diagonalProb) / (dimension - 1));
    MatrixDouble_setDiagonalValue(transition->matrix, diagonalProb);
    transition->transitionCountData = TransitionCountData_construct(numberOfStates);
    transition->numberOfValidityFunctions = 0;
    transition->numberOfStates = numberOfStates;
    transition->requirements = NULL;
    transition->validityFunctions = NULL;
    return transition;
}

void Transition_destruct(Transition *transition) {
    MatrixDouble_destruct(transition->matrix);
    TransitionCountData_destruct(transition->transitionCountData);
    TransitionRequirements_destruct(transition->requirements);
    free(transition->validityFunctions);
    free(transition);
}

void Transition_addRequirements(Transition *transition, TransitionRequirements *requirements){
    transition->requirements = requirements;
}

void Transition_addValidityFunction(Transition *transition, ValidityFunction validityFunction){
    transition->validityFunctions = realloc(transition->validityFunctions, (transition->numberOfValidityFunctions + 1) * sizeof(ValidityFunction));
    transition->validityFunctions[transition->numberOfValidityFunctions] = validityFunction;
    transition->numberOfValidityFunctions += 1;
}


bool ValidityFunction_checkDupByMapq(StateType state, CoverageInfo *coverageInfo, TransitionRequirements *requirements){
    double highMapqRatio = coverageInfo->coverage_high_mapq / coverageInfo->coverage;
    if ((state == STATE_DUP) && (highMapqRatio > requirements->maxHighMapqRatio)){
        return false;
    }
    return true;
}

bool ValidityFunction_checkMsjByClipping(StateType state, CoverageInfo *coverageInfo, TransitionRequirements *requirements){
    double highlyClippedRatio = coverageInfo->coverage_high_clip / coverageInfo->coverage;
    if ((state == STATE_MSJ) && (highlyClippedRatio < requirements->minHighlyClippedRatio)){
        return false;
    }
    return true;
}

bool Transition_isStateValid(Transition *transition, StateType state, CoverageInfo *coverageInfo){
    for(int i=0; i < transition->numberOfValidityFunctions; i++){
        ValidityFunction isValid = transition->validityFunctions[i];
        if (isValid(state, coverageInfo, transition->requirements) == false){
            return false;
        }
    }
    return true;
}

double Transition_getProb(Transition *transition, StateType preState, StateType state){
    return transition->matrix->data[preState][state];
}

double Transition_getTerminationProb(Transition *transition, StateType state){
    return transition->matrix->data[state][transition->numberOfStates];
}

double Transition_getStartProb(Transition *transition, StateType state){
    return transition->matrix->data[transition->numberOfStates][state];
}

double Transition_getProbConditional(Transition *transition, StateType preState, StateType state, CoverageInfo *coverageInfo){
    double totProbValid = 0.0;
    for(int s=0; s < transition->numberOfStates + 1; s++){
        if (Transition_isStateValid(transition, s, coverageInfo)) {
            totProbValid += Transition_getProb(transition, preState, s);
        }
    }
    double prob = Transition_getProb(transition, preState, state);
    if (Transition_isStateValid(transition, state, coverageInfo)){
        return prob / totProbValid;
    }else{
        return 0.0;
    }
}

double Exponential_getPdf(double x, double lam) {
    if (x < 0.0){
        return 0.0;
    }else{
        return lam * exp(-lam * x);
    }
}

double log_trunc_exp_pdf(double x, double lam, double b) {
    if (x < 0.0 || x > b) {
        return -DBL_MAX;
    }
    else{
        return Exponential_getPdf((double) x, lam) / (1.0 - exp(-lam * b));
    }
}
