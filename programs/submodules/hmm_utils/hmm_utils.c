//
// Created by mobin on 1/8/24.
//

#include "hmm_utils.h"
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "ptBlock.h"
#include "stdbool.h"
#include "common.h"
#include "digamma.h"


ModelType getModelTypeFromString(const char* modelString){
        if (strcmp(modelString, "nb") == 0 ||
            strcmp(modelString, "negative_binomial") == 0){
                return MODEL_NEGATIVE_BINOMIAL;
        }
        else if (strcmp(modelString, "gaussian") == 0){
                return MODEL_GAUSSIAN;
        }
        else if(strcmp(modelString, "trunc_exp_gaussian") == 0 ||
                strcmp(modelString, "truncated_exponential_gaussian") == 0 ){
                return MODEL_TRUNC_EXP_GAUSSIAN;
        }
        return MODEL_UNDEFINED;
}


//////////////////////////////////////
//// ParameterEstimator Functions ////
/////////////////////////////////////



ParameterEstimator *ParameterEstimator_construct(EmissionDist * emissionDist, int numberOfComps){
    ParameterEstimator *parameterEstimator = malloc(1 * sizeof(ParameterEstimator));
    parameterEstimator->numeratorPerComp = Double_construct1DArray(numberOfComps);
    parameterEstimator->denominatorPerComp = Double_construct1DArray(numberOfComps);
    parameterEstimator->numberOfComps = numberOfComps;
    parameterEstimator->emissionDist = emissionDist;
    parameterEstimator->mutexPtr = malloc(sizeof(pthread_mutex_t));
    pthread_mutex_init(parameterEstimator->mutexPtr, NULL);
    return parameterEstimator;
}

void ParameterEstimator_increment(ParameterEstimator *parameterEstimator, double numerator, double denominator, int compIndex){
    pthread_mutex_lock(parameterEstimator->mutexPtr);
    parameterEstimator->numeratorPerComp[compIndex] += numerator;
    parameterEstimator->denominatorPerComp[compIndex] += denominator;
    pthread_mutex_unlock(parameterEstimator->mutexPtr);
}

void ParameterEstimator_incrementDenominatorForAllComps(ParameterEstimator *parameterEstimator, double numerator, double denominator, int compIndex){
    pthread_mutex_lock(parameterEstimator->mutexPtr);
    parameterEstimator->numeratorPerComp[compIndex] += numerator;
    for(int i=0; i < parameterEstimator->numberOfComps; i++) {
        parameterEstimator->denominatorPerComp[i] += denominator;
    }
    pthread_mutex_unlock(parameterEstimator->mutexPtr);
}

double ParameterEstimator_getEstimation(ParameterEstimator *parameterEstimator, int compIndex, double *count){
    *count = parameterEstimator->denominatorPerComp[compIndex];
    if(parameterEstimator->denominatorPerComp[compIndex] == 0){
        fprintf(stderr, "Warning: the denominator is zero for comp=%d, so the estimation will be returned zero\n", compIndex);
        return 0.0;
    }
    double est;
    if (parameterEstimator->emissionDist->distType == DIST_TRUNC_EXPONENTIAL){
        est = TruncExponential_estimateLambda((TruncExponential *) parameterEstimator->emissionDist->dist,
                                              parameterEstimator,
                                              1e-6);
    }else {
        est = parameterEstimator->numeratorPerComp[compIndex] / parameterEstimator->denominatorPerComp[compIndex];
    }
    return est;
}

void ParameterEstimator_reset(ParameterEstimator *parameterEstimator){
    Double_fill1DArray(parameterEstimator->numeratorPerComp, parameterEstimator->numberOfComps, 0.0);
    Double_fill1DArray(parameterEstimator->denominatorPerComp, parameterEstimator->numberOfComps, 0.0);
}


void ParameterEstimator_destruct(ParameterEstimator *parameterEstimator){
    Double_destruct1DArray(parameterEstimator->numeratorPerComp);
    Double_destruct1DArray(parameterEstimator->denominatorPerComp);
    free(parameterEstimator->mutexPtr);
    free(parameterEstimator);
}




/////////////////////////////////////
//// ParameterBinding Functions ////
////////////////////////////////////




ParameterBinding *ParameterBinding_construct(double **coefs, int numberOfParams, int numberOfComps){
    ParameterBinding *parameterBinding = malloc(1 * sizeof(ParameterBinding));
    parameterBinding->coefs = Double_copy2DArray(coefs, numberOfParams, numberOfComps);
    parameterBinding->numberOfComps = numberOfComps;
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
    parameterBinding->numberOfComps = numberOfComps;
    parameterBinding->numberOfParams = numberOfParams;
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


ParameterBinding **ParameterBinding_getDefault1DArrayForModel(ModelType modelType, int numberOfCollapsedComps, bool excludeMisjoin){
    if(modelType == MODEL_GAUSSIAN){
        return ParameterBinding_getDefault1DArrayForGaussian(numberOfCollapsedComps, excludeMisjoin);
    }else if (modelType == MODEL_TRUNC_EXP_GAUSSIAN){
        return ParameterBinding_getDefault1DArrayForTruncExpGaussian(numberOfCollapsedComps, excludeMisjoin);
    }else if (modelType == MODEL_NEGATIVE_BINOMIAL){
        return ParameterBinding_getDefault1DArrayForNegativeBinomial(numberOfCollapsedComps, excludeMisjoin);
    }
}

ParameterBinding **ParameterBinding_getDefault1DArrayForGaussian(int numberOfCollapsedComps, bool excludeMisjoin) {
    int numberOfParams = 3; // mean and var
    int numberOfStates = excludeMisjoin ? NUMBER_OF_STATES - 1 : NUMBER_OF_STATES;
    ParameterBinding **parameterBinding1DArray = malloc(numberOfStates * sizeof(ParameterBinding *));

    // create one array and fill them for each state separately
    double *coefsPerParam = Double_construct1DArray(numberOfParams);

    coefsPerParam[GAUSSIAN_MEAN] = ERR_COMP_BINDING_COEF;
    coefsPerParam[GAUSSIAN_VAR] = ERR_COMP_BINDING_COEF;
    coefsPerParam[GAUSSIAN_WEIGHT] = 0.0;
    parameterBinding1DArray[STATE_ERR] = ParameterBinding_constructForSingleComp(coefsPerParam, numberOfParams);

    coefsPerParam[GAUSSIAN_MEAN] = 0.5;
    coefsPerParam[GAUSSIAN_VAR] = 0.5;
    coefsPerParam[GAUSSIAN_WEIGHT] = 0.0;
    parameterBinding1DArray[STATE_DUP] = ParameterBinding_constructForSingleComp(coefsPerParam, numberOfParams);

    coefsPerParam[GAUSSIAN_MEAN] = 1.0;
    coefsPerParam[GAUSSIAN_VAR] = 1.0;
    coefsPerParam[GAUSSIAN_WEIGHT] = 0.0;
    parameterBinding1DArray[STATE_HAP] = ParameterBinding_constructForSingleComp(coefsPerParam, numberOfParams);

    double *stepsPerParam = Double_construct1DArray(numberOfParams);

    coefsPerParam[GAUSSIAN_MEAN] = 2.0;
    coefsPerParam[GAUSSIAN_VAR] = 2.0;
    coefsPerParam[GAUSSIAN_WEIGHT] = 0.0;

    stepsPerParam[GAUSSIAN_MEAN] = 1.0;
    stepsPerParam[GAUSSIAN_VAR] = 1.0;
    stepsPerParam[GAUSSIAN_WEIGHT] = 0.0;
    parameterBinding1DArray[STATE_COL] = ParameterBinding_constructSequenceByStep(coefsPerParam,
                                                                                  stepsPerParam,
                                                                                  numberOfParams,
                                                                                  numberOfCollapsedComps);

    if (! excludeMisjoin) {
        coefsPerParam[GAUSSIAN_MEAN] = 1.0;
        coefsPerParam[GAUSSIAN_VAR] = 1.0;
        coefsPerParam[GAUSSIAN_WEIGHT] = 0.0;
        parameterBinding1DArray[STATE_MSJ] = ParameterBinding_constructForSingleComp(coefsPerParam, numberOfParams);
    }

    Double_destruct1DArray(coefsPerParam);
    Double_destruct1DArray(stepsPerParam);
    return parameterBinding1DArray;
}

ParameterBinding **ParameterBinding_getDefault1DArrayForNegativeBinomial(int numberOfCollapsedComps, bool excludeMisjoin){
    int numberOfParams = 3; // theta and lambda
    int numberOfStates = excludeMisjoin ? NUMBER_OF_STATES - 1 : NUMBER_OF_STATES;
    ParameterBinding **parameterBinding1DArray = malloc(numberOfStates * sizeof(ParameterBinding *));

    // create one array and fill them for each state separately
    double *coefsPerParam = Double_construct1DArray(numberOfParams);

    coefsPerParam[NB_THETA] = 1.0;
    coefsPerParam[NB_LAMBDA] = ERR_COMP_BINDING_COEF;
    coefsPerParam[NB_WEIGHT] = 0;
    parameterBinding1DArray[STATE_ERR] = ParameterBinding_constructForSingleComp(coefsPerParam, numberOfParams);

    coefsPerParam[NB_THETA] = 1.0;
    coefsPerParam[NB_LAMBDA] = 0.5;
    coefsPerParam[NB_WEIGHT] = 0;
    parameterBinding1DArray[STATE_DUP] = ParameterBinding_constructForSingleComp(coefsPerParam, numberOfParams);

    coefsPerParam[NB_THETA] = 1.0;
    coefsPerParam[NB_LAMBDA] = 1.0;
    coefsPerParam[NB_WEIGHT] = 0;
    parameterBinding1DArray[STATE_HAP] = ParameterBinding_constructForSingleComp(coefsPerParam, numberOfParams);

    double *stepsPerParam = Double_construct1DArray(numberOfParams);

    coefsPerParam[NB_THETA] = 1.0;
    coefsPerParam[NB_LAMBDA] = 2.0;
    coefsPerParam[NB_WEIGHT] = 0;

    stepsPerParam[NB_THETA] = 0.0;
    stepsPerParam[NB_LAMBDA] = 1.0;
    stepsPerParam[NB_WEIGHT] = 0;
    parameterBinding1DArray[STATE_COL] = ParameterBinding_constructSequenceByStep(coefsPerParam,
                                                                                  stepsPerParam,
                                                                                  numberOfParams,
                                                                                  numberOfCollapsedComps);

    if (! excludeMisjoin) {
        coefsPerParam[NB_THETA] = 0;
        coefsPerParam[NB_LAMBDA] = 1.0;
        coefsPerParam[NB_WEIGHT] = 0;
        parameterBinding1DArray[STATE_MSJ] = ParameterBinding_constructForSingleComp(coefsPerParam, numberOfParams);
    }

    Double_destruct1DArray(coefsPerParam);
    Double_destruct1DArray(stepsPerParam);
    return parameterBinding1DArray;
}

ParameterBinding **ParameterBinding_getDefault1DArrayForTruncExpGaussian(int numberOfCollapsedComps, bool excludeMisjoin){
    ParameterBinding **parameterBinding1DArray = ParameterBinding_getDefault1DArrayForGaussian(numberOfCollapsedComps, excludeMisjoin);
    ParameterBinding_destruct(parameterBinding1DArray[STATE_ERR]);

    int numberOfParams = 1;
    double *coefsPerParam = Double_construct1DArray(numberOfParams);

    coefsPerParam[TRUNC_EXP_LAMBDA] = 0.0; // zero means no binding
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






////////////////////////////////////
//// NegativeBinomial Functions ////
////////////////////////////////////





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
    nb->theta = Double_construct1DArray(numberOfComps);
    nb->lambda = Double_construct1DArray(numberOfComps);
    for(int comp=0;comp < numberOfComps; comp++){
        nb->theta[comp] = NegativeBinomial_getTheta(mean[comp], var[comp]);
        nb->lambda[comp] = NegativeBinomial_getLambda(mean[comp], var[comp]);
    }
    // Allocating and initializing mixture weights
    nb->weights = Double_construct1DArray(numberOfComps);
    Double_fill1DArray(nb->weights, numberOfComps, 1.0 / numberOfComps);
    nb->numberOfComps = numberOfComps;
    NegativeBinomial_fillDigammaTable(nb);
    // wrap nb in EmissionDist for initializing estimator
    nb->thetaEstimator = NULL;
    nb->lambdaEstimator = NULL;
    nb->weightsEstimator = NULL;
    return nb;
}

void NegativeBinomial_fillDigammaTable(NegativeBinomial *nb) {
    nb->digammaTable = Double_construct2DArray(nb->numberOfComps, MAX_COVERAGE_VALUE + 1);
    for (int comp = 0; comp < nb->numberOfComps; comp++) {
        double r = NegativeBinomial_getR(nb->theta[comp], nb->lambda[comp]);
        double digammal_0 = digammal(r); // digamma(r + x) for x = 0
        nb->digammaTable[comp][0] = digammal_0;
        for (int x = 1; x <= MAX_COVERAGE_VALUE; x++) {
            // digamma(1 + z) = digamma(z) + 1 / z
            // or digamma(r + x) = digamma(r + x - 1) + 1 / (r + x - 1)
            nb->digammaTable[comp][x] = nb->digammaTable[comp][x - 1] + 1.0 / (r + x - 1);
        }
    }
}


/**
 * Destructs a Negative Binomial object
 *
 * @param nb
 */
void NegativeBinomial_destruct(NegativeBinomial *nb) {
    free(nb->theta);
    free(nb->lambda);
    free(nb->weights);
    ParameterEstimator_destruct(nb->thetaEstimator);
    ParameterEstimator_destruct(nb->lambdaEstimator);
    ParameterEstimator_destruct(nb->weightsEstimator);
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

double NegativeBinomial_getLambda(double mean, double var) {
    double r = pow(mean, 2) / (var - mean);
    double lambda = -1 * r * log(NegativeBinomial_getTheta(mean, var));
    return lambda;
}

double NegativeBinomial_getR(double theta, double lambda) {
    double r = -1 * lambda / log(theta);
    return r;
}

double NegativeBinomial_getMean(double theta, double lambda) {
    double r = -1 * lambda / log(theta);
    double mean = r * (1 - theta) / theta;
    return mean;
}

double NegativeBinomial_getVar(double theta, double lambda) {
    double r = -1 * lambda / log(theta);
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
    for (int comp = 0; comp < nb->numberOfComps; comp++) {
        theta = nb->theta[comp];
        r = NegativeBinomial_getR(nb->theta[comp], nb->lambda[comp]);
        w = nb->weights[comp];
        probs[comp] = w * exp(lgamma(r + x) - lgamma(r) - lgamma(x + 1) + r * log(theta) +
                           (double) x * log(1 - theta));
        if (probs[comp] != probs[comp])
	{
            fprintf(stderr, "prob is NAN lambda=%.2e, r=%.2e, theta=%.2e, w=%.2e\n", nb->lambda[comp], r, theta, w);
            exit(EXIT_FAILURE);
        }
        if (probs[comp] < 1e-40) {
            probs[comp] = 1e-40;
        }
    }
    return probs;
}


ParameterEstimator *NegativeBinomial_getEstimator(NegativeBinomial *nb, NegativeBinomialParameterType parameterType){
    if (parameterType == NB_THETA){
        return nb->thetaEstimator;
    }else if (parameterType == NB_LAMBDA){
        return nb->lambdaEstimator;
    }else if(parameterType == NB_WEIGHT){
        return nb->weightsEstimator;
    }
}



void NegativeBinomial_updateEstimator(NegativeBinomial *nb,
                                         uint8_t x,
                                         double count){
    double *componentProbs = NegativeBinomial_getComponentProbs(nb, x);
    double totProb = Double_sum1DArray(componentProbs, nb->numberOfComps);
    for (int comp = 0; comp < nb->numberOfComps; comp++) {
        double theta = nb->theta[comp];
        double r = NegativeBinomial_getR(nb->theta[comp], nb->lambda[comp]);
        double beta = -1 * theta / (1 - theta) - 1 / log(theta);
        double w = count * componentProbs[comp] / totProb;
        //w2 = em->f[i][c] * em->b[i][c] * em->scales[i] / terminateProb;
        // Update sufficient stats for estimating mean vectors
        double delta = r * (nb->digammaTable[comp][x] - nb->digammaTable[comp][0]);
        ParameterEstimator_increment(nb->lambdaEstimator,
                                     w * delta,
                                     w,
                                     comp);
        ParameterEstimator_increment(nb->thetaEstimator,
                                     w * delta * beta,
                                     w * delta * beta + w * (x - delta),
                                     comp);
        ParameterEstimator_incrementDenominatorForAllComps(nb->weightsEstimator,
                                                           w,
                                                           w,
                                                           comp);
    }
    free(componentProbs);
}

bool NegativeBinomial_updateParameter(NegativeBinomial *nb, NegativeBinomialParameterType parameterType, int compIndex, double value, double convergenceTol){
    double oldValue = 0.0;
    if (parameterType == NB_THETA){
	oldValue = nb->theta[compIndex];
        nb->theta[compIndex] = value;
    }else if (parameterType == NB_LAMBDA){
	oldValue = nb->lambda[compIndex];
        nb->lambda[compIndex] = value;
    }else if (parameterType == NB_WEIGHT){
	oldValue = nb->weights[compIndex];
        nb->weights[compIndex] = value;
    }
    double diffRatio = 1.0e-4 < oldValue ? fabs(value / oldValue - 1.0) : 0.0;
    //fprintf(stdout, "%.4e ,%.4e, %.4e, %.4e\n", value, oldValue, diffRatio, convergenceTol);
    bool converged = diffRatio < convergenceTol;
    return converged;
}

double *NegativeBinomial_getParameterValues(NegativeBinomial *nb, NegativeBinomialParameterType parameterType) {
    double *paramValues = Double_construct1DArray(nb->numberOfComps);
    for(int comp=0; comp < nb->numberOfComps; comp++){
        if (parameterType == NB_LAMBDA)
            paramValues[comp] = nb->lambda[comp];
        else if (parameterType == NB_THETA)
            paramValues[comp] = nb->theta[comp];
        else if (parameterType == NB_WEIGHT)
            paramValues[comp] = nb->weights[comp];
        else if (parameterType == NB_MEAN)
            paramValues[comp] = NegativeBinomial_getMean(nb->theta[comp], nb->lambda[comp]);
        else if (parameterType == NB_VAR)
            paramValues[comp] = NegativeBinomial_getVar(nb->theta[comp], nb->lambda[comp]);
        else
            exit(EXIT_FAILURE);
    }
    return paramValues;
}

const char *NegativeBinomial_getParameterName(NegativeBinomialParameterType parameterType){
    return NegativeBinomialParameterToString[parameterType];
}


/////////////////////////////
//// Gaussian Functions ////
////////////////////////////





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
    gaussian->mean = Double_construct1DArray(numberOfComps);
    gaussian->var = Double_construct1DArray(numberOfComps);
    memcpy(gaussian->mean, mean, numberOfComps * sizeof(double));
    memcpy(gaussian->var, var, numberOfComps * sizeof(double));
    // Allocating and initializing mixture weights
    gaussian->weights = Double_construct1DArray(numberOfComps);
    Double_fill1DArray(gaussian->weights, numberOfComps, 1.0 / numberOfComps);
    gaussian->numberOfComps = numberOfComps;
    // wrap gaussian in EmissionDist for initializing estimators
    gaussian->meanEstimator = NULL;
    gaussian->varEstimator = NULL;
    gaussian->weightsEstimator = NULL;
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
    ParameterEstimator_destruct(gaussian->meanEstimator);
    ParameterEstimator_destruct(gaussian->varEstimator);
    ParameterEstimator_destruct(gaussian->weightsEstimator);
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
    for (int comp = 0; comp < numberOfComps; comp++) {
        var[comp] = mean[comp] * factor;
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
    for (int comp = 0;  comp < gaussian->numberOfComps; comp++) {
        mean = (1 - alpha) * gaussian->mean[comp] + alpha * preX;
        var = gaussian->var[comp];
        w = gaussian->weights[comp];
        // adjust the mean value based on the previous observation and alpha (dependency factor)
        probs[comp] = w / (sqrt(var * 2 * PI)) * exp(-0.5 * pow((x - mean), 2) / var);
        if (probs[comp] != probs[comp]) {
            double u = -0.5 * pow((x - mean), 2) / var;
            fprintf(stderr, "[Error] prob is NAN exp(%.2e) mean=%.2e, var=%.2e\n", u, mean, var);
            exit(EXIT_FAILURE);
        }
        if (probs[comp] < 1e-40) {
            //fprintf(stderr, "[Warning] prob is lower than 1e-40. It is set to 1e-40 [comp = %d] mean=%.3f, var=%.3f, alpha=%.3f, preX=%.3f, x=%.3f, w=%.3f\n", comp,mean, var, alpha, (double) preX, (double)x, w);
            probs[comp] = 1e-40;
        }
    }
    return probs;
}

ParameterEstimator *Gaussian_getEstimator(Gaussian *gaussian, GaussianParameterType parameterType){
    if (parameterType == GAUSSIAN_MEAN){
        return gaussian->meanEstimator;
    }else if (parameterType == GAUSSIAN_VAR){
        return gaussian->varEstimator;
    }else if(parameterType == GAUSSIAN_WEIGHT){
        return gaussian->weightsEstimator;
    }
}


void Gaussian_updateEstimator(Gaussian *gaussian,
                                 uint8_t x,
                                 uint8_t preX,
                                 double alpha,
                                 double count){
    double x_adjusted = (x - alpha * preX) / (1 - alpha);
    double *componentProbs = Gaussian_getComponentProbs(gaussian, x, preX, alpha);
    double totProb = Double_sum1DArray(componentProbs, gaussian->numberOfComps);
    for (int c = 0; c < gaussian->numberOfComps; c++) {
        double w = count * componentProbs[c] / totProb;
        ParameterEstimator_increment(gaussian->meanEstimator,
                                     w * x_adjusted,
                                     w,
                                     c);
        double z = x_adjusted - gaussian->mean[c];
        ParameterEstimator_increment(gaussian->varEstimator,
                                     w * z * z,
                                     w,
                                     c);
        ParameterEstimator_incrementDenominatorForAllComps(gaussian->weightsEstimator,
                                                           w,
                                                           w,
                                                           c);
    }
    free(componentProbs);
}


bool Gaussian_updateParameter(Gaussian *gaussian, GaussianParameterType parameterType, int compIndex, double value, double convergenceTol){
    double oldValue = 0.0;
    if (parameterType == GAUSSIAN_MEAN){
	oldValue = gaussian->mean[compIndex];
        gaussian->mean[compIndex] = value;
    }else if (parameterType == GAUSSIAN_VAR){
	oldValue = gaussian->var[compIndex];
        gaussian->var[compIndex] = value;
    }else if(parameterType == GAUSSIAN_WEIGHT){
	oldValue = gaussian->weights[compIndex];
        gaussian->weights[compIndex] = value;
    }
    double diffRatio = 1.0e-4 < oldValue ? fabs(value / oldValue - 1.0) : 0.0;
    //fprintf(stdout, "%.4e ,%.4e, %.4e, %.4e\n", value, oldValue, diffRatio, convergenceTol);
    bool converged = diffRatio < convergenceTol;
    return converged; 
}

double *Gaussian_getParameterValues(Gaussian *gaussian, GaussianParameterType parameterType) {
    double *paramValues = Double_construct1DArray(gaussian->numberOfComps);
    for(int comp=0; comp < gaussian->numberOfComps; comp++){
        if (parameterType == GAUSSIAN_MEAN)
            paramValues[comp] = gaussian->mean[comp];
        else if (parameterType == GAUSSIAN_VAR)
            paramValues[comp] = gaussian->var[comp];
        else if (parameterType == NB_WEIGHT)
            paramValues[comp] = gaussian->weights[comp];
        else
            exit(EXIT_FAILURE);
    }
    return paramValues;
}

const char *Gaussian_getParameterName(GaussianParameterType parameterType){
    return GaussianParameterToString[parameterType];
}


////////////////////////////////////
//// TruncExponential Functions ////
////////////////////////////////////




TruncExponential *TruncExponential_construct(double lambda, double truncPoint){
    TruncExponential *truncExponential = malloc(sizeof(TruncExponential));
    truncExponential->lambda = lambda;
    truncExponential->truncPoint = truncPoint;
    truncExponential->lambdaEstimator = NULL; // wrap TruncExponential in EmissionDist for initializing estimator
    return truncExponential;
}

void TruncExponential_destruct(TruncExponential* truncExponential){
    ParameterEstimator_destruct(truncExponential->lambdaEstimator);
    free(truncExponential);
}

double TruncExponential_getProb(TruncExponential* truncExponential, uint8_t x) {
    double lam = truncExponential->lambda;
    double b = truncExponential->truncPoint;
    if (x < 0.0 || truncExponential->truncPoint < x)
	    return 0.0;
    return lam * exp(-lam * x) / (1 - exp(-lam * b));
}

double TruncExponential_getLogLikelihoodByParams(double lambda, double truncPoint, ParameterEstimator *lambdaEstimator) {
    double num = lambdaEstimator->numeratorPerComp[0];
    double denom = lambdaEstimator->denominatorPerComp[0];
    double lam = lambda;
    double b = truncPoint;
    return denom * log(lam) - denom * log(1.0 - exp(-lam * b)) - num * lam;
}

double TruncExponential_getLogLikelihood(TruncExponential *truncExponential, ParameterEstimator *lambdaEstimator) {
    double num = lambdaEstimator->numeratorPerComp[0];
    double denom = lambdaEstimator->denominatorPerComp[0];
    double lam = truncExponential->lambda;
    double b = truncExponential->truncPoint;
    return denom * log(lam) - denom * log(1.0 - exp(-lam * b)) - num * lam;
}

// taken from Jordan's script
// which is originally adapted from https://en.wikipedia.org/wiki/Golden-section_search
// tol = 1e-6
double TruncExponential_estimateLambda(TruncExponential *truncExponential, ParameterEstimator *lambdaEstimator, double tol) {
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
    double yc = TruncExponential_getLogLikelihoodByParams(c, truncExponential->truncPoint, lambdaEstimator); // f(c)
    double yd = TruncExponential_getLogLikelihoodByParams(d, truncExponential->truncPoint, lambdaEstimator); // f(d)
    for (int k = 0; k < n - 1; k++) {
        if (yc > yd) {
            b = d;
            d = c;
            yd = yc;
            h = invphi * h;
            c = a + invphi2 * h;
            yc = TruncExponential_getLogLikelihoodByParams(c, truncExponential->truncPoint, lambdaEstimator); // f(c)
        } else {
            a = c;
            c = d;
            yc = yd;
            h = invphi * h;
            d = a + invphi * h;
            yd = TruncExponential_getLogLikelihoodByParams(d, truncExponential->truncPoint, lambdaEstimator); //f(d)
        }
    }

    if (yc > yd) {
        return (a + d) / 2.0;
    } else {
        return (c + b) / 2.0;
    }
}


ParameterEstimator *TruncExponential_getEstimator(TruncExponential *truncExponential, TruncatedExponentialParameterType parameterType){
    if (parameterType == TRUNC_EXP_LAMBDA) {
        return truncExponential->lambdaEstimator;
    }
}


void TruncExponential_updateEstimator(TruncExponential *truncExponential,
                                         uint8_t x,
                                         double count){
    ParameterEstimator_increment(truncExponential->lambdaEstimator,
                                 count * x,
                                 count,
                                 0);
}

bool TruncExponential_updateParameter(TruncExponential *truncExponential, TruncatedExponentialParameterType parameterType, double value, double convergenceTol){
    double oldLambda = 0.0;
    double newLambda = value;
    if (parameterType == TRUNC_EXP_LAMBDA){
	oldLambda = truncExponential->lambda;
        truncExponential->lambda = value;
	//truncExponential->truncPoint = value / ERR_COMP_BINDING_COEF * EXP_TRUNC_POINT_COV_FRACTION;
	//fprintf(stderr, "truncExponential->truncPoint = %.2e, truncExponential->lambda= %.2e\n", truncExponential->truncPoint, truncExponential->lambda);
    }
    if (parameterType == TRUNC_EXP_TRUNC_POINT){
	    truncExponential->truncPoint = value;
	    return true;
    }
    double diffRatio = 1.0e-4 < oldLambda ? fabs(newLambda / oldLambda - 1.0) : 0.0;
    bool converged = diffRatio < convergenceTol;
    return converged;
}

double *TruncExponential_getParameterValues(TruncExponential *truncExponential, TruncatedExponentialParameterType parameterType){
    double *paramValues = Double_construct1DArray(1);
    if (parameterType == TRUNC_EXP_LAMBDA){
        paramValues[0] = truncExponential->lambda;
    }
    else if (parameterType == TRUNC_EXP_MEAN){
        paramValues[0] = 1.0 / truncExponential->lambda;
    }
    else if (parameterType == TRUNC_EXP_TRUNC_POINT){
	paramValues[0] = truncExponential->truncPoint;
    }
    else{
        exit(EXIT_FAILURE);
    }
    return paramValues;
}

const char *TruncExponential_getParameterName(TruncatedExponentialParameterType parameterType){
    return TruncatedExponentialParameterToString[parameterType];
}


/////////////////////////////////
//// EmissionDist Functions ////
////////////////////////////////




EmissionDist *EmissionDist_construct(void* dist, DistType distType){
    EmissionDist *emissionDist = malloc(1 * sizeof(EmissionDist));
    emissionDist->dist = dist;
    emissionDist->distType = distType;
    EmissionDist_initParameterEstimators(emissionDist);
    return emissionDist;
}

void EmissionDist_initParameterEstimators(EmissionDist *emissionDist){
    if(emissionDist->distType == DIST_TRUNC_EXPONENTIAL){
        TruncExponential *truncExponential = (TruncExponential *) emissionDist->dist;
        truncExponential->lambdaEstimator = ParameterEstimator_construct(emissionDist, 1);
    }else if(emissionDist->distType == DIST_GAUSSIAN){
        Gaussian *gaussian = (Gaussian *) emissionDist->dist;
        gaussian->meanEstimator = ParameterEstimator_construct(emissionDist, gaussian->numberOfComps);
        gaussian->varEstimator = ParameterEstimator_construct(emissionDist, gaussian->numberOfComps);
        gaussian->weightsEstimator = ParameterEstimator_construct(emissionDist, gaussian->numberOfComps);
    }else if(emissionDist->distType == DIST_NEGATIVE_BINOMIAL){
        NegativeBinomial *nb = (NegativeBinomial*) emissionDist->dist;
        nb->thetaEstimator = ParameterEstimator_construct(emissionDist, nb->numberOfComps);
        nb->lambdaEstimator = ParameterEstimator_construct(emissionDist, nb->numberOfComps);
        nb->weightsEstimator = ParameterEstimator_construct(emissionDist, nb->numberOfComps);
    }
}

void EmissionDist_resetParameterEstimators(EmissionDist *emissionDist){
    if(emissionDist->distType == DIST_TRUNC_EXPONENTIAL){
        TruncExponential *truncExponential = (TruncExponential *) emissionDist->dist;
        ParameterEstimator_reset(truncExponential->lambdaEstimator);
    }else if(emissionDist->distType == DIST_GAUSSIAN){
        Gaussian *gaussian = (Gaussian *) emissionDist->dist;
        ParameterEstimator_reset(gaussian->meanEstimator);
        ParameterEstimator_reset(gaussian->varEstimator);
        ParameterEstimator_reset(gaussian->weightsEstimator);
    }else if(emissionDist->distType == DIST_NEGATIVE_BINOMIAL){
        NegativeBinomial *nb = (NegativeBinomial*) emissionDist->dist;
        ParameterEstimator_reset(nb->thetaEstimator);
        ParameterEstimator_reset(nb->lambdaEstimator);
        ParameterEstimator_reset(nb->weightsEstimator);
    }
}


int EmissionDist_getNumberOfComps(EmissionDist* emissionDist){
    if (emissionDist->distType == DIST_TRUNC_EXPONENTIAL){
        return 1;
    }
    else if (emissionDist->distType == DIST_GAUSSIAN){
        Gaussian* gaussian = (Gaussian *) emissionDist->dist;
        return gaussian->numberOfComps;
    }
    else if (emissionDist->distType == DIST_NEGATIVE_BINOMIAL){
        NegativeBinomial * nb = (NegativeBinomial *) emissionDist->dist;
        return nb->numberOfComps;
    }
}

ParameterEstimator *EmissionDist_getEstimator(EmissionDist* emissionDist, void *parameterTypePtr){
    if (emissionDist->distType == DIST_TRUNC_EXPONENTIAL){
        TruncatedExponentialParameterType *parameterTypeTruncExpPtr = parameterTypePtr;
        return TruncExponential_getEstimator((TruncExponential *) emissionDist->dist,
                                             *parameterTypeTruncExpPtr);
    }
    else if (emissionDist->distType == DIST_GAUSSIAN){
        GaussianParameterType *parameterTypeGaussianPtr = parameterTypePtr;
        return Gaussian_getEstimator((Gaussian *) emissionDist->dist,
                                     *parameterTypeGaussianPtr);
    }
    else if (emissionDist->distType == DIST_NEGATIVE_BINOMIAL){
        NegativeBinomialParameterType *parameterTypeNBPtr = parameterTypePtr;
        return NegativeBinomial_getEstimator((NegativeBinomial *) emissionDist->dist,
                                             *parameterTypeNBPtr);
    }else {
        return NULL;
    }
}

bool EmissionDist_updateParameter(EmissionDist* emissionDist, void *parameterTypePtr, int compIndex, double value, double convergenceTol){
    if (emissionDist->distType == DIST_TRUNC_EXPONENTIAL){
        TruncatedExponentialParameterType *parameterTypeTruncExpPtr = parameterTypePtr;
        return TruncExponential_updateParameter((TruncExponential *) emissionDist->dist,
                                                *parameterTypeTruncExpPtr,
                                                value,
						convergenceTol);
    }
    else if (emissionDist->distType == DIST_GAUSSIAN){
        GaussianParameterType *parameterTypeGaussianPtr = parameterTypePtr;
        return Gaussian_updateParameter((Gaussian *) emissionDist->dist,
                                        *parameterTypeGaussianPtr,
                                        compIndex,
                                        value,
					convergenceTol);
    }
    else if (emissionDist->distType == DIST_NEGATIVE_BINOMIAL){
        NegativeBinomialParameterType *parameterTypeNBPtr = parameterTypePtr;
        return NegativeBinomial_updateParameter((NegativeBinomial *) emissionDist->dist,
                                                *parameterTypeNBPtr,
                                                compIndex,
                                                value,
						convergenceTol);
    }
    return true;
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

void EmissionDist_updateEstimator(EmissionDist *emissionDist, uint8_t x, uint8_t preX, double alpha, double count) {
    if (emissionDist->distType == DIST_TRUNC_EXPONENTIAL){
        TruncExponential_updateEstimator((TruncExponential *) emissionDist->dist, x, count);
    }
    else if (emissionDist->distType == DIST_GAUSSIAN){
        Gaussian_updateEstimator((Gaussian *) emissionDist->dist, x, preX, alpha, count);
    }
    else if (emissionDist->distType == DIST_NEGATIVE_BINOMIAL){
        NegativeBinomial_updateEstimator((NegativeBinomial *) emissionDist->dist, x, count);
    }
}

double *EmissionDist_getParameterValuesForOneType(EmissionDist* emissionDist, void *parameterTypePtr){
    if (emissionDist->distType == DIST_TRUNC_EXPONENTIAL){
        TruncatedExponentialParameterType *parameterTypeTruncExpPtr = parameterTypePtr;
        return TruncExponential_getParameterValues((TruncExponential *) emissionDist->dist,
                                                   *parameterTypeTruncExpPtr);
    }
    else if (emissionDist->distType == DIST_GAUSSIAN){
        GaussianParameterType *parameterTypeGaussianPtr = parameterTypePtr;
        return Gaussian_getParameterValues((Gaussian *) emissionDist->dist,
                                           *parameterTypeGaussianPtr);
    }
    else if (emissionDist->distType == DIST_NEGATIVE_BINOMIAL){
        NegativeBinomialParameterType *parameterTypeNBPtr = parameterTypePtr;
        return NegativeBinomial_getParameterValues((NegativeBinomial *) emissionDist->dist,
                                                   *parameterTypeNBPtr);
    }else {
        return NULL;
    }
}

const char *EmissionDist_getParameterName(EmissionDist* emissionDist, void *parameterTypePtr){
    if (emissionDist->distType == DIST_TRUNC_EXPONENTIAL){
        TruncatedExponentialParameterType *parameterTypeTruncExpPtr = parameterTypePtr;
        return TruncExponential_getParameterName(*parameterTypeTruncExpPtr);
    }
    else if (emissionDist->distType == DIST_GAUSSIAN){
        GaussianParameterType *parameterTypeGaussianPtr = parameterTypePtr;
        return Gaussian_getParameterName(*parameterTypeGaussianPtr);
    }
    else if (emissionDist->distType == DIST_NEGATIVE_BINOMIAL){
        NegativeBinomialParameterType *parameterTypeNBPtr = parameterTypePtr;
        return NegativeBinomial_getParameterName(*parameterTypeNBPtr);
    }else {
        return "NA";
    }
}

int EmissionDist_getParameterIndex(EmissionDist* emissionDist, void *parameterTypePtr){
    if (emissionDist->distType == DIST_TRUNC_EXPONENTIAL){
        TruncatedExponentialParameterType *parameterTypeTruncExpPtr = parameterTypePtr;
        return (int) *parameterTypeTruncExpPtr;
    }
    else if (emissionDist->distType == DIST_GAUSSIAN){
        GaussianParameterType *parameterTypeGaussianPtr = parameterTypePtr;
        return (int) *parameterTypeGaussianPtr;
    }
    else if (emissionDist->distType == DIST_NEGATIVE_BINOMIAL){
        NegativeBinomialParameterType *parameterTypeNBPtr = parameterTypePtr;
        return (int) *parameterTypeNBPtr;
    }else {
        return -1;
    }
}


void **EmissionDist_getParameterTypePtrsForLogging(EmissionDist* emissionDist, int *length){
    if (emissionDist->distType == DIST_TRUNC_EXPONENTIAL){
        TruncatedExponentialParameterType **parameterTypeTruncExpPtrs = malloc(2*sizeof(TruncatedExponentialParameterType*));
        for(int i=0; i < 2; i++){
	    parameterTypeTruncExpPtrs[i] = malloc(sizeof(TruncatedExponentialParameterType));
	}
        parameterTypeTruncExpPtrs[0][0] = TRUNC_EXP_MEAN;
	parameterTypeTruncExpPtrs[1][0] = TRUNC_EXP_TRUNC_POINT;
        *length = 2;
        return (void**) parameterTypeTruncExpPtrs;
    }
    else if (emissionDist->distType == DIST_GAUSSIAN){
        GaussianParameterType **parameterTypeGaussianPtrs = malloc(3*sizeof(GaussianParameterType*));
        for(int i=0; i < 3; i++){
            parameterTypeGaussianPtrs[i] = malloc(sizeof(GaussianParameterType));
        }
        parameterTypeGaussianPtrs[0][0] = GAUSSIAN_MEAN;
        parameterTypeGaussianPtrs[1][0] = GAUSSIAN_VAR;
        parameterTypeGaussianPtrs[2][0] = GAUSSIAN_WEIGHT;
        *length = 3;
        return (void**) parameterTypeGaussianPtrs;
    }
    else if (emissionDist->distType == DIST_NEGATIVE_BINOMIAL){
        NegativeBinomialParameterType **parameterTypeNegativeBinomialPtrs = malloc(3*sizeof(NegativeBinomialParameterType*));
        for(int i=0; i < 3; i++){
            parameterTypeNegativeBinomialPtrs[i] = malloc(sizeof(NegativeBinomialParameterType));
        }
        parameterTypeNegativeBinomialPtrs[0][0] = NB_MEAN;
        parameterTypeNegativeBinomialPtrs[1][0] = NB_VAR;
        parameterTypeNegativeBinomialPtrs[2][0] = NB_WEIGHT;
        *length =3;
        return (void**) parameterTypeNegativeBinomialPtrs;
    }else {
        *length = 0;
        return NULL;
    }
}

const char**EmissionDist_getParameterNames(EmissionDist* emissionDist, void **parameterTypePtrs, int numberOfParams){
    const char **parameterNames = malloc(numberOfParams * sizeof(char*));
    for(int i=0; i < numberOfParams; i++){
         parameterNames[i] = EmissionDist_getParameterName(emissionDist, parameterTypePtrs[i]);
    }
    return parameterNames;
}

double **EmissionDist_getParameterValues(EmissionDist* emissionDist, void **parameterTypePtrs, int numberOfParams){
    double **parameterValues = malloc(numberOfParams * sizeof(double*));
    for(int i=0; i < numberOfParams; i++){
        parameterValues[i] = EmissionDist_getParameterValuesForOneType(emissionDist, parameterTypePtrs[i]);
    }
    return parameterValues;
}

const char* EmissionDist_getDistributionName(EmissionDist* emissionDist){
    return DistToString[emissionDist->distType];
}

//////////////////////////////////////
//// EmissionDistSeries Functions ////
//////////////////////////////////////





EmissionDistSeries *EmissionDistSeries_constructForModel(ModelType modelType,
                                                         double **means,  // [numberOfDists] x [maxMixtures]
                                                         int *numberOfCompsPerDist,
                                                         int numberOfDists,
							 bool excludeMisjoin){
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
            dist = NegativeBinomial_constructByMean(means[s], 1.5, numberOfCompsPerDist[s]);
            emissionDistSeries->emissionDists[s] = EmissionDist_construct(dist, DIST_NEGATIVE_BINOMIAL);
        }
    }

    emissionDistSeries->parameterBindingPerDist =
            ParameterBinding_getDefault1DArrayForModel(modelType,
                                                       numberOfCompsPerDist[STATE_COL],
						       excludeMisjoin);
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
                                     emissionDistSeries->numberOfDists);
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


int EmissionDistSeries_getNumberOfComps(EmissionDistSeries * emissionDistSeries, int distIndex){
    EmissionDist *emissionDist = emissionDistSeries->emissionDists[distIndex];
    return EmissionDist_getNumberOfComps(emissionDist);
}

void EmissionDistSeries_updateEstimator(EmissionDistSeries *emissionDistSeries, int distIndex, uint8_t x, uint8_t preX, double alpha, double count) {
    EmissionDist *emissionDist = emissionDistSeries->emissionDists[distIndex];
    EmissionDist_updateEstimator(emissionDist, x, preX, alpha, count);
}

void EmissionDistSeries_resetParameterEstimators(EmissionDistSeries *emissionDistSeries){
	for (int distIndex = 0; distIndex < emissionDistSeries->numberOfDists; distIndex++) {
        	EmissionDist *emissionDist = emissionDistSeries->emissionDists[distIndex];
		EmissionDist_resetParameterEstimators(emissionDist);
	}
}

ParameterEstimator *EmissionDistSeries_getBoundParameterEstimator(EmissionDistSeries *emissionDistSeries, DistType distType, void *parameterTypePtr) {
    EmissionDist *mockEmissionDist = EmissionDist_construct(NULL, DIST_UNDEFINED);
    ParameterEstimator *boundParameterEstimator = ParameterEstimator_construct(mockEmissionDist, 1);
    for (int distIndex = 0; distIndex < emissionDistSeries->numberOfDists; distIndex++) {
        EmissionDist *emissionDist = emissionDistSeries->emissionDists[distIndex];
        if (emissionDist->distType != distType){
            continue;
        }
	int paramIndex = EmissionDist_getParameterIndex(emissionDist, parameterTypePtr);
        ParameterBinding *parameterBinding = emissionDistSeries->parameterBindingPerDist[distIndex];
        ParameterEstimator *parameterEstimator = EmissionDist_getEstimator(emissionDist, parameterTypePtr);
        for (int comp = 0; comp < parameterBinding->numberOfComps; comp++) {
            double factor = parameterBinding->coefs[paramIndex][comp];
            if (0.0 < factor) {
                ParameterEstimator_increment(boundParameterEstimator,
                                             parameterEstimator->numeratorPerComp[comp] / factor,
                                             parameterEstimator->denominatorPerComp[comp],
                                             0);
            }
        }
    }
    return boundParameterEstimator;
}

bool EmissionDistSeries_estimateOneParameterType(EmissionDistSeries *emissionDistSeries, DistType distType, void *parameterTypePtr, double convergenceTol) {
    bool converged = true;
    ParameterEstimator *parameterEstimatorBound = EmissionDistSeries_getBoundParameterEstimator(emissionDistSeries, distType, parameterTypePtr);
    double boundCount;
    double boundEstimation = ParameterEstimator_getEstimation(parameterEstimatorBound,
                                                              0,
							      &boundCount);
    for(int distIndex=0; distIndex < emissionDistSeries->numberOfDists; distIndex++){
        EmissionDist *emissionDist = emissionDistSeries->emissionDists[distIndex];
        if (emissionDist->distType != distType){
            continue;
        }
	int paramIndex = EmissionDist_getParameterIndex(emissionDist, parameterTypePtr);
        ParameterBinding *parameterBinding = emissionDistSeries->parameterBindingPerDist[distIndex];
        for(int comp=0; comp < parameterBinding->numberOfComps; comp++) {
            double factor = parameterBinding->coefs[paramIndex][comp];
            double estimation;
	    double count;
            if (0.0 <  factor) {
                estimation = boundEstimation * factor;
		count = boundCount;
            }else{
                ParameterEstimator *parameterEstimator = EmissionDist_getEstimator(emissionDist, parameterTypePtr);
                estimation = ParameterEstimator_getEstimation(parameterEstimator, comp, &count);
            }
	    //fprintf(stderr, "distIndex = %d, param=%s -> count = %.2e\n", distIndex, EmissionDist_getParameterName(emissionDist, parameterTypePtr), count );
	    if (MIN_COUNT_FOR_PARAMETER_UPDATE < count){
		// if at least one parameter is not converged then "converged" variable will become false
                converged &= EmissionDist_updateParameter(emissionDist,
                                             parameterTypePtr,
                                             comp,
                                             estimation,
					     convergenceTol);
	    }
        }
    }
    return converged;
}

bool EmissionDistSeries_estimateParameters(EmissionDistSeries *emissionDistSeries, double convergenceTol){
    bool converged = true;
    if(emissionDistSeries->modelType == MODEL_GAUSSIAN || emissionDistSeries->modelType == MODEL_TRUNC_EXP_GAUSSIAN){
        GaussianParameterType *parameterTypeArrayGaussian = malloc(3 * sizeof(GaussianParameterType));
        parameterTypeArrayGaussian[0] = GAUSSIAN_MEAN;
        parameterTypeArrayGaussian[1] = GAUSSIAN_VAR;
        parameterTypeArrayGaussian[2] = GAUSSIAN_WEIGHT;
        for (int i=0; i < 3; i++) {
            converged &= EmissionDistSeries_estimateOneParameterType(emissionDistSeries, DIST_GAUSSIAN, &parameterTypeArrayGaussian[i], convergenceTol);
        }
        free(parameterTypeArrayGaussian);
        if(emissionDistSeries->modelType == MODEL_TRUNC_EXP_GAUSSIAN){
            TruncatedExponentialParameterType *parameterTypeArrayTruncExp = malloc(1 * sizeof(TruncatedExponentialParameterType));
            parameterTypeArrayTruncExp[0] = TRUNC_EXP_LAMBDA;
            converged &= EmissionDistSeries_estimateOneParameterType(emissionDistSeries, DIST_TRUNC_EXPONENTIAL, &parameterTypeArrayTruncExp[0], convergenceTol);
	    TruncExponential *errDist = (TruncExponential *) emissionDistSeries->emissionDists[STATE_ERR]->dist;
	    Gaussian *hapDist = (Gaussian *) emissionDistSeries->emissionDists[STATE_HAP]->dist;
	    double newTruncPoint = hapDist->mean[0] * EXP_TRUNC_POINT_COV_FRACTION;
	    //fprintf(stderr, "NEW TRUNC POINT = %.2e, %.2e, %.2e\n", newTruncPoint, hapDist->mean[0], EXP_TRUNC_POINT_COV_FRACTION); 
	    TruncExponential_updateParameter(errDist, TRUNC_EXP_TRUNC_POINT, newTruncPoint, convergenceTol);
            free(parameterTypeArrayTruncExp);
        }
    }else if(emissionDistSeries->modelType == MODEL_NEGATIVE_BINOMIAL){
        NegativeBinomialParameterType *parameterTypeArrayNB = malloc(3 * sizeof(NegativeBinomialParameterType));
        parameterTypeArrayNB[0] = NB_THETA;
        parameterTypeArrayNB[1] = NB_LAMBDA;
        parameterTypeArrayNB[2] = NB_WEIGHT;
        for (int i=0; i < 3; i++) {
            converged &= EmissionDistSeries_estimateOneParameterType(emissionDistSeries, DIST_NEGATIVE_BINOMIAL, &parameterTypeArrayNB[i], convergenceTol);
        }
	// update digamma table
	for(int distIndex=0; distIndex < emissionDistSeries->numberOfDists; distIndex++){
		EmissionDist *emissionDist = emissionDistSeries->emissionDists[distIndex];
		NegativeBinomial *nb = emissionDist->dist;
		NegativeBinomial_fillDigammaTable(nb);
	}
        free(parameterTypeArrayNB);
    }
    return converged;
}

const char**EmissionDistSeries_getParameterNames(EmissionDistSeries *emissionDistSeries, int distIndex, int* numberOfParams) {
    EmissionDist *emissionDist = emissionDistSeries->emissionDists[distIndex];
    void **parameterTypePtrs = EmissionDist_getParameterTypePtrsForLogging(emissionDist,numberOfParams);
    return EmissionDist_getParameterNames(emissionDist, parameterTypePtrs, *numberOfParams);
}

double **EmissionDistSeries_getParameterValues(EmissionDistSeries *emissionDistSeries, int distIndex, int* numberOfParams) {
    EmissionDist *emissionDist = emissionDistSeries->emissionDists[distIndex];
    void **parameterTypePtrs = EmissionDist_getParameterTypePtrsForLogging(emissionDist,numberOfParams);
    return EmissionDist_getParameterValues(emissionDist, parameterTypePtrs, *numberOfParams);
}

const char* EmissionDistSeries_getDistributionName(EmissionDistSeries *emissionDistSeries, int distIndex) {
    return EmissionDist_getDistributionName(emissionDistSeries->emissionDists[distIndex]);
}

const char* EmissionDistSeries_getStateName(int distIndex) {
    return StateToString[distIndex];
}

//////////////////////////////////////////
//// TransitionRequirements Functions ////
/////////////////////////////////////////



TransitionRequirements *TransitionRequirements_construct(double minHighlyClippedRatio, double maxHighMapqRatio){
    TransitionRequirements *transitionRequirements = malloc(sizeof(TransitionRequirements));
    transitionRequirements->minHighlyClippedRatio = minHighlyClippedRatio;
    transitionRequirements->maxHighMapqRatio = maxHighMapqRatio;
    return transitionRequirements;
}

void TransitionRequirements_destruct(TransitionRequirements *transitionRequirements){
    free(transitionRequirements);
}




////////////////////////////////////////
//// TransitionCountData Functions ////
///////////////////////////////////////




TransitionCountData *TransitionCountData_construct(int numberOfStates){
    TransitionCountData *transitionCountData = malloc(sizeof(TransitionCountData));
    transitionCountData->countMatrix = MatrixDouble_construct0(numberOfStates + 1, numberOfStates + 1);
    transitionCountData->pseudoCountMatrix = MatrixDouble_construct0(numberOfStates + 1, numberOfStates + 1);
    transitionCountData->numberOfStates = numberOfStates;
    transitionCountData->mutexPtr = malloc(sizeof(pthread_mutex_t));
    pthread_mutex_init(transitionCountData->mutexPtr, NULL);
    return transitionCountData;
}

void TransitionCountData_increment(TransitionCountData *transitionCountData, double count, StateType preState, StateType state){
    pthread_mutex_lock(transitionCountData->mutexPtr);
    transitionCountData->countMatrix->data[preState][state] += count;
    pthread_mutex_unlock(transitionCountData->mutexPtr);
}


void TransitionCountData_resetCountMatrix(TransitionCountData *transitionCountData){
    MatrixDouble_setValue(transitionCountData->countMatrix, 0.0);
}

void TransitionCountData_parsePseudoCountFromFile(TransitionCountData *transitionCountData, char* pathToMatrix, int dim){
    assert(dim == (transitionCountData->numberOfStates + 1));
    MatrixDouble_destruct(transitionCountData->pseudoCountMatrix);
    transitionCountData->pseudoCountMatrix = MatrixDouble_parseFromFile(pathToMatrix, dim ,dim, false);
}

void TransitionCountData_setPseudoCountMatrix(TransitionCountData *transitionCountData, double value){
    MatrixDouble_setValue(transitionCountData->pseudoCountMatrix, value);
}

TransitionCountData *TransitionCountData_destruct(TransitionCountData *transitionCountData){
    free(transitionCountData->countMatrix);
    free(transitionCountData->pseudoCountMatrix);
    free(transitionCountData->mutexPtr);
    free(transitionCountData);
}




//////////////////////////////
//// Transition Functions ////
//////////////////////////////




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
    transition->terminationProb = 1e-3;
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
    transition->terminationProb = 1e-3;
    return transition;
}


void Transition_resetCountData(Transition *transition) {
    TransitionCountData_resetCountMatrix(transition->transitionCountData);
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

bool Transition_estimateTransitionMatrix(Transition *transition, double convergenceTol){
    bool converged = true;
    double terminationProb = transition->terminationProb;
    MatrixDouble *countMatrix = transition->transitionCountData->countMatrix;
    MatrixDouble *pseudoCountMatrix = transition->transitionCountData->pseudoCountMatrix;
    for(int i1=0; i1 < transition->numberOfStates; i1++){
        double rowSum = 0.0;
        // take sum per row
        for(int i2=0; i2 < transition->numberOfStates; i2++){
            rowSum += countMatrix->data[i1][i2] + pseudoCountMatrix->data[i1][i2];
        }
        // estimate transition
        for(int i2=0; i2 < transition->numberOfStates; i2++){
	    double oldValue = transition->matrix->data[i1][i2];
	    double newValue = (countMatrix->data[i1][i2] + pseudoCountMatrix->data[i1][i2]) / rowSum * (1.0 - terminationProb);
            transition->matrix->data[i1][i2] = newValue;
	    double diffRatio = 1.0e-6 < oldValue ? fabs(newValue / oldValue - 1.0) : 0.0;
	    converged &= diffRatio < convergenceTol;

        }
    }

    // termination probs
    for(int i1=0; i1 < transition->numberOfStates; i1++){
            transition->matrix->data[i1][transition->numberOfStates] = terminationProb;
    }
    // start probs
    for(int i2=0; i2 < transition->numberOfStates; i2++){
	    transition->matrix->data[transition->numberOfStates][i2] = 1.0 / transition->numberOfStates;
    }
    // the probability of empty sequence
    transition->matrix->data[transition->numberOfStates][transition->numberOfStates] = 0.0;
    return converged;
}


/////////////////////////////
//// Validity Functions ////
////////////////////////////




bool ValidityFunction_checkDupByMapq(StateType state, CoverageInfo *coverageInfo, TransitionRequirements *requirements){
    double highMapqRatio = (double) coverageInfo->coverage_high_mapq / (0.1 + coverageInfo->coverage);
    if ((state == STATE_DUP) && (highMapqRatio > requirements->maxHighMapqRatio)){
        return false;
    }
    return true;
}

bool ValidityFunction_checkMsjByClipping(StateType state, CoverageInfo *coverageInfo, TransitionRequirements *requirements){
    double highlyClippedRatio = (double) coverageInfo->coverage_high_clip / (0.1 + coverageInfo->coverage);
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
