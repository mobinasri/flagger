#ifndef HMM_UTILS_H
#define HMM_UTILS_H

#include <stdint.h>
#include "data_types.h"
#include "ptBlock.h"
#include "common.h"
#include "count_data.h"



#define MIN_COUNT_FOR_PARAMETER_UPDATE 10
#define EXP_TRUNC_POINT_COV_FRACTION 0.25
#define NUMBER_OF_STATES 5
#define ERR_COMP_BINDING_COEF 0.1
#define MAX_COVERAGE_VALUE 2048
#define START_STATE_INDEX 5
#define END_STATE_INDEX 5

typedef enum StateType {
    STATE_ERR = 0,
    STATE_DUP = 1,
    STATE_HAP = 2,
    STATE_COL = 3,
    STATE_MSJ = 4
} StateType;

typedef enum DependencyType {
    ALPHA_ERR = 0,
    ALPHA_DUP = 1,
    ALPHA_HAP = 2,
    ALPHA_COL = 3,
    ALPHA_TRN = 4,
    ALPHA_MSJ = 5
} DependencyType;

typedef enum DistType {
    DIST_TRUNC_EXPONENTIAL = 0,
    DIST_GAUSSIAN = 1,
    DIST_NEGATIVE_BINOMIAL = 2,
    DIST_UNDEFINED = 3
}DistType;

typedef enum ModelType {
    MODEL_TRUNC_EXP_GAUSSIAN = 0,
    MODEL_GAUSSIAN = 1,
    MODEL_NEGATIVE_BINOMIAL = 2,
    MODEL_UNDEFINED = 3
} ModelType;

typedef enum GaussianParameterType{
    GAUSSIAN_MEAN = 0,
    GAUSSIAN_VAR = 1,
    GAUSSIAN_WEIGHT = 2
} GaussianParameterType;

typedef enum NegativeBinomialParameterType{
    NB_THETA = 0,
    NB_LAMBDA = 1,
    NB_WEIGHT = 2,
    NB_MEAN = 3,
    NB_VAR = 4
} NegativeBinomialParameterType;

typedef enum TruncatedExponentialParameterType{
    TRUNC_EXP_LAMBDA = 0,
    TRUNC_EXP_MEAN = 1,
    TRUNC_EXP_TRUNC_POINT = 2
} TruncatedExponentialParameterType;

static const char* NegativeBinomialParameterToString[5] = {"Theta", "Lambda", "Weight", "Mean", "Var"};
static const char* GaussianParameterToString[3] = {"Mean", "Var", "Weight"};
static const char* TruncatedExponentialParameterToString[3] = {"Rate", "Mean", "Trunc_Point"};
static const char* DistToString[4] = {"Truncated Exponential", "Gaussian", "Negative Binomial", "Undefined"};
static const char* StateToString[5] = {"Err", "Dup", "Hap", "Col" ,"Msj"};

ModelType getModelTypeFromString(const char* modelString);

/*! @typedef
 * @abstract Structure for storing a distribution structure. It can be either Gaussian,
 * Negative Binominal, or TruncExponential.
 * @field dist          The distribution structure
 * @field distType      The type of distribution
 */
typedef struct EmissionDist {
    void *dist;
    DistType distType;
}EmissionDist;


typedef struct ParameterEstimator{
    double *numeratorPerComp;
    double *denominatorPerComp;
    int numberOfComps;
    EmissionDist *emissionDist;
    // ParameterEstimator can be updated with any order
    // of observations
    // since ParameterEstimator might be shared between
    // multiple threads for faster runtime it is necessary
    // to keep a mutex per estimator
    pthread_mutex_t *mutexPtr;
}ParameterEstimator;

ParameterEstimator *ParameterEstimator_construct(EmissionDist * emissionDist, int numberOfComps);
void ParameterEstimator_increment(ParameterEstimator *parameterEstimator, double numerator, double denominator, int compIndex);
double ParameterEstimator_getEstimation(ParameterEstimator *parameterEstimator, int compIndex, double *count);
void ParameterEstimator_reset(ParameterEstimator *parameterEstimator);
void ParameterEstimator_destruct(ParameterEstimator *parameterEstimator);


/*! @typedef
 * @abstract Structure for representing the relationship between parameters of different components
 *           within a distribution. It also facilitates binding parameters of one distribution to another.
 *           This requires a 1D array of ParameterBinding structs, each storing coefficients for one
 *           distribution.
 * @field coefs             A 2D array of coefficients. The first dimension represents the parameters,
 *                          and the second dimension represents the components in the distribution.
 * @field numberOfComps     Number of components in the distribution.
 * @field numberOfParams    Number of parameters in the distribution.
 */
typedef struct ParameterBinding{
    double **coefs; // [numberOfParams] x [numberOfComps]
    int numberOfComps;
    int numberOfParams;
} ParameterBinding;

/*
 * Constructs a ParameterBinding struct
 */
ParameterBinding *ParameterBinding_construct(double **coefs, int numberOfParams, int numberOfComps);

/*
 * Constructs a ParameterBinding struct
 */
ParameterBinding *ParameterBinding_constructSequenceByStep(double* firstCoefsPerParam,double* stepsPerParam,
                                                           int numberOfParams,int numberOfComps);

/*
 * Constructs a ParameterBinding struct for a distribution with a single component
 */
ParameterBinding *ParameterBinding_constructForSingleComp(double* coefsPerParam, int numberOfParams);

/*
 * Get the binding coefficient for a parameter of one component
 */
double ParameterBinding_getValue(ParameterBinding* parameterBinding, int paramIndex, int compIndex);

/*
 * Make a 1D array of ParameterBinding struct for the specified model type. The only state with more than
 * component is the "collapsed" component so that it needs the number of components for that.
 */
ParameterBinding **ParameterBinding_getDefault1DArrayForModel(ModelType modelType, int numberOfCollapsedComps, bool excludeMisjoin);

/*
 * Make a 1D array of ParameterBinding struct for the Gaussian model. The only state with more than
 * component is the "collapsed" component so that it needs the number of components for that.
 */
ParameterBinding **ParameterBinding_getDefault1DArrayForGaussian(int numberOfCollapsedComps, bool excludeMisjoin);

/*
 * Make a 1D array of ParameterBinding struct for the Negative Binomial model. The only state with more than
 * component is the "collapsed" component so that it needs the number of components for that.
 */
ParameterBinding **ParameterBinding_getDefault1DArrayForNegativeBinomial(int numberOfCollapsedComps, bool excludeMisjoin);

/*
 * Make a 1D array of ParameterBinding struct for the Truncated Exponential + Guassial model.
 * The only state with more than component is the "collapsed" component so that it needs the
 * number of components for that.
 */
ParameterBinding **ParameterBinding_getDefault1DArrayForTruncExpGaussian(int numberOfCollapsedComps, bool excludeMisjoin);

/*
 * Destruct a ParameterBinding struct
 */
void ParameterBinding_destruct(ParameterBinding* parameterBinding);

/*
 * Destruct a 1D array of ParameterBinding struct
 */
void ParameterBinding_destruct1DArray(ParameterBinding **parameterBinding1DArray, int arrayLength);


/*! @typedef
 * @abstract Structure for storing the parameters of a mixture of Negative Binomial distributions
 * @field mean              mean values of components
 * @field var               variance values of components
 * @field theta             theta values of components
 * @field r                 r values of components
 * @field weights           weights of components
 * @field numberOfComps     number of mixture components
 */
typedef struct NegativeBinomial {
    double *theta; // theta
    double *lambda; // lambda
    double *weights; // weights of components
    ParameterEstimator *lambdaEstimator;
    ParameterEstimator *thetaEstimator; // r = -1 * lambda / log(theta);
    ParameterEstimator *weightsEstimator;
    int numberOfComps; // number of mixture components
    double **digammaTable; // [number of comps] x [max coverage + 1]
} NegativeBinomial;

/*
 * Construct a Negative Binomial structure
 */
NegativeBinomial *NegativeBinomial_construct(double *mean, double *var, int numberOfComps);

void NegativeBinomial_fillDigammaTable(NegativeBinomial *nb);

/*
 * Construct a Negative Binomial structure using an array of means and a factor to scale mean for getting variances
 */
NegativeBinomial *NegativeBinomial_constructByMean(double *mean, double factor, int numberOfComps);

/*
 * Calculate theta given mean and variance
 */
double NegativeBinomial_getTheta(double mean, double var);

/*
 * Calculate lambda given mean and variance
 */
double NegativeBinomial_getLambda(double mean, double var);

/*
 * Calculate R given theta and lambda
 */
double NegativeBinomial_getR(double theta, double lambda);

/*
 * Calculate mean given theta and lambda
 */
double NegativeBinomial_getMean(double theta, double lambda);

/*
 * Calculate variance given theta and lambda
 */
double NegativeBinomial_getVar(double theta, double lambda);

/*
 * Get the probability of observing x
 */
double NegativeBinomial_getProb(NegativeBinomial *nb, uint8_t x);

/*
 * Get the probabilities of observing x from different components
 */
double *NegativeBinomial_getComponentProbs(NegativeBinomial *nb, uint8_t x);

ParameterEstimator *NegativeBinomial_getEstimator(NegativeBinomial *nb, NegativeBinomialParameterType parameterType);

void NegativeBinomial_updateEstimator(NegativeBinomial *nb, uint8_t x, double count);

bool NegativeBinomial_updateParameter(NegativeBinomial *nb, NegativeBinomialParameterType parameterType, int compIndex, double value, double convergenceTol);

double *NegativeBinomial_getParameterValues(NegativeBinomial *nb,NegativeBinomialParameterType parameterType);

const char *NegativeBinomial_getParameterName(NegativeBinomialParameterType parameterType);

/*
 * Destruct a Negative Binomial structure
 */
void NegativeBinomial_destruct(NegativeBinomial *nb);


/*! @typedef
 * @abstract Structure for storing the parameters of a mixture of Gaussian distributions
 * @field mean              mean values of components
 * @field var               variance values of components
 * @field weights           weights of components
 * @field numberOfComps     number of mixture components
 */
typedef struct Gaussian {
    double *mean;
    double *var;
    double *weights; // weights of components
    ParameterEstimator *meanEstimator;
    ParameterEstimator *varEstimator;
    ParameterEstimator *weightsEstimator;
    int numberOfComps; // number of mixture components
} Gaussian;

/*
 * Construct a Gaussian structure
 */
Gaussian *Gaussian_construct(double *mean, double *var, int numberOfComps);

/*
 * Construct a Gaussian structure using an array of means and a factor to scale mean for getting variances
 */
Gaussian *Gaussian_constructByMean(double *mean, double factor, int numberOfComps);

/*
 * Get the probability of observing x given the previous observation preX and alpha the factor for adjusting
 * mean of Gaussian with the previous observation
 */
double Gaussian_getProb(Gaussian *gaussian, uint8_t x, uint8_t preX, double alpha);

/*
 * Get the probability of observing x from different components given the previous observation preX
 * and alpha the factor for adjusting mean of Gaussian with the previous observation
 */
double *Gaussian_getComponentProbs(Gaussian *gaussian, uint8_t x, uint8_t preX, double alpha);

ParameterEstimator *Gaussian_getEstimator(Gaussian *gaussian, GaussianParameterType parameterType);

void Gaussian_updateEstimator(Gaussian *gaussian,
                              uint8_t x,
                              uint8_t preX,
                              double alpha,
                              double count);

bool Gaussian_updateParameter(Gaussian *gaussian, GaussianParameterType parameterType, int compIndex, double value, double convergenceTol);

double *Gaussian_getParameterValues(Gaussian *gaussian, GaussianParameterType parameterType);

const char *Gaussian_getParameterName(GaussianParameterType parameterType);


/*
 * Destruct a Gaussian structure
 */
void Gaussian_destruct(Gaussian *gaussian);


/*! @typedef
 * @abstract Structure for storing the parameters of a truncated exponential distribution
 * @field lambda
 * @field truncPoint    The truncation point of the distribution
 */
typedef struct TruncExponential {
    double lambda; // lambda or rate = 1/mean
    double truncPoint; // truncation point
    ParameterEstimator *lambdaEstimator;
} TruncExponential;

/*
 * Destruct a TruncExponential structure
 */
TruncExponential * TruncExponential_construct(double lambda, double truncPoint);

/*
 * Get the probability of observing x
 */
double TruncExponential_getProb(TruncExponential* truncExponential, uint8_t x);

/*
 * Get the log likelihood of observing the count data (equivalent to a list of values)
 * given a TruncExponential structure
 */
double TruncExponential_getLogLikelihood(TruncExponential *truncExponential, ParameterEstimator *parameterEstimator);

/*
 * Get the log likelihood of observing the count data (equivalent to a list of values) given
 *  the parameters of a truncated exponential distribution
 */
double TruncExponential_getLogLikelihoodByParams(double lambda, double truncPoint, ParameterEstimator *parameterEstimator);


/*
 * Estimate lambda given a TruncExponential structure, a countData and a tolerance for specifying the accuracy
 * of the estimated lambda. It will use golden section search for estimating lambda.
 */
double TruncExponential_estimateLambda(TruncExponential *truncExponential, ParameterEstimator *parameterEstimator, double tol);

ParameterEstimator *TruncExponential_getEstimator(TruncExponential *truncExponential, TruncatedExponentialParameterType parameterType);

void TruncExponential_updateEstimator(TruncExponential *truncExponential,
                                      uint8_t x,
                                      double count);

bool TruncExponential_updateParameter(TruncExponential *truncExponential, TruncatedExponentialParameterType parameterType, double value, double convergenceTol);

double *TruncExponential_getParameterValues(TruncExponential *truncExponential, TruncatedExponentialParameterType parameterType);
const char *TruncExponential_getParameterName(TruncatedExponentialParameterType parameterType);


/*
 * Destruct a TruncExponential structure
 */
void TruncExponential_destruct(TruncExponential* truncExponential);


/*
 * Construct a EmissionDist structure
 */
EmissionDist *EmissionDist_construct(void* dist, DistType distType);

void EmissionDist_initParameterEstimators(EmissionDist *emissionDist);

int EmissionDist_getNumberOfComps(EmissionDist* emissionDist);

/*
 * Get the probability of observing x given preX, the previous observation and alpha, the factor for adjusting
 * mean of distribution with the previous observation (if applicable). Alpha and preX are only considered for
 * Gaussian distribution and ignored for NegativeBinomial and TruncExponential
 */
double EmissionDist_getProb(EmissionDist *emissionDist, uint8_t x, uint8_t preX, double alpha);

bool EmissionDist_updateParameter(EmissionDist* emissionDist, void *parameterTypePtr, int compIndex, double value, double convergenceTol);

ParameterEstimator *EmissionDist_getEstimator(EmissionDist* emissionDist, void *parameterTypePtr);

void EmissionDist_updateEstimator(EmissionDist *emissionDist, uint8_t x, uint8_t preX, double alpha, double count);

void **EmissionDist_getParameterTypePtrsForLogging(EmissionDist* emissionDist, int *length);

double *EmissionDist_getParameterValuesForOneType(EmissionDist* emissionDist, void *parameterTypePtr);

const char *EmissionDist_getParameterName(EmissionDist* emissionDist, void *parameterTypePtr);

const char**EmissionDist_getParameterNames(EmissionDist* emissionDist, void **parameterTypePtrs, int numberOfParams);

double **EmissionDist_getParameterValues(EmissionDist* emissionDist, void **parameterTypePtrs, int numberOfParams);

const char* EmissionDist_getDistributionName(EmissionDist* emissionDist);

void EmissionDist_resetParameterEstimators(EmissionDist *emissionDist);

/*
 * Destruct a EmissionDist structure
 */
void EmissionDist_destruct(EmissionDist* emissionDist);


/*! @typedef
 * @abstract Structure for storing a series of EmissionDist structures and its related attributes
 * for estimating parameters
 * @field emissionDists                 An array of EmissionDist structures
 * @field countDataPerDist              An array of CountData structures, each of which is holding observed counts for
 *                                      a distribution saved in emissionDists.
 * @field parameterBindingPerDist       An array of ParameterBinding structures, each of which is holding the binding
 *                                      coefficients for parameters and components of a distribution saved in
 *                                      emissionDists.
 * @field numberOfDists                 Number of distributions
 * @field modelType                     The type of the model. It can be either
 *                                      MODEL_TRUNC_EXP_GAUSSIAN, MODEL_GAUSSIAN, or MODEL_NEGATIVE_BINOMIAL
 */
typedef struct EmissionDistSeries {
    EmissionDist **emissionDists;
    CountData **countDataPerDist;
    ParameterBinding **parameterBindingPerDist;
    int numberOfDists;
    ModelType modelType;
}EmissionDistSeries;


/*
 * Construct a EmissionDistSeries structure given the type of the model, the mean values of different
 * distributions and the number of components per distribution
 */
EmissionDistSeries *EmissionDistSeries_constructForModel(ModelType modelType,
                                                         double **means,  // [numberOfDists] x [maxMixtures]
                                                         int *numberOfCompsPerDist,
                                                         int numberOfDists,
							 bool excludeMisjoin);

/*
 * Get a emissionDist structure from the emissionDistSeries structure
 */
EmissionDist *EmissionDistSeries_getEmissionDist(EmissionDistSeries* emissionDistSeries, int distIndex);

void EmissionDistSeries_resetParameterEstimators(EmissionDistSeries *emissionDistSeries);

/*
 * Destruct an EmissionDistSeries structure
 */
void EmissionDistSeries_destruct(EmissionDistSeries* emissionDistSeries);

/*
 * Get the array of probabilities of observing x for all distributions given preX, the previous observation
 * and alpha, the factor for adjusting mean of distribution with the previous observation (if applicable).
 * Alpha and preX are only considered for Gaussian distribution and ignored for NegativeBinomial and TruncExponential
 */
double *EmissionDistSeries_getAllProbs(EmissionDistSeries* emissionDistSeries,
                                       uint8_t x,
                                       uint8_t preX,
                                       double alpha);

/*
 * Get the probability of observing x for a specific distribution given preX, the previous observation
 * and alpha, the factor for adjusting mean of distribution with the previous observation (if applicable).
 * Alpha and preX are only considered for Gaussian distribution and ignored for NegativeBinomial and TruncExponential
 */
double EmissionDistSeries_getProb(EmissionDistSeries* emissionDistSeries,
                                  int distIndex,
                                  uint8_t x,
                                  uint8_t preX,
                                  double alpha);



int EmissionDistSeries_getNumberOfComps(EmissionDistSeries * emissionDistSeries, int distIndex);

void EmissionDistSeries_updateEstimator(EmissionDistSeries *emissionDistSeries,
                                         int distIndex,
                                         uint8_t x,
                                         uint8_t preX,
                                         double alpha,
                                         double count);

ParameterEstimator *EmissionDistSeries_getBoundParameterEstimator(EmissionDistSeries *emissionDistSeries,
                                                                  DistType distType,
                                                                  void *parameterTypePtr);

bool EmissionDistSeries_estimateOneParameterType(EmissionDistSeries *emissionDistSeries,
                                                 DistType distType,
                                                 void *parameterTypePtr,
						 double convergenceTol);

bool EmissionDistSeries_estimateParameters(EmissionDistSeries *emissionDistSeries, double convergenceTol);

const char ** EmissionDistSeries_getParameterNames(EmissionDistSeries *emissionDistSeries, int distIndex, int* numberOfParams);

double ** EmissionDistSeries_getParameterValues(EmissionDistSeries *emissionDistSeries, int distIndex, int* numberOfParams);

const char * EmissionDistSeries_getDistributionName(EmissionDistSeries *emissionDistSeries, int distIndex);

const char* EmissionDistSeries_getStateName(int distIndex);

/*! @typedef
 * @abstract Structure for holding the attributes that are used for restricting the transitions in an HMM model
 * @field minHighlyClippedRatio     The minimum ratio of the coverage of highly clipped alignments over total coverage
 *                                  for allowing a transition to the misjoin state
 * @field maxHighMapqRatio          The maximum ratio of the coverage of alignments with high mapq over total coverage
 *                                  for allowing a transition to the duplicated state
 */
typedef struct TransitionRequirements{
    double minHighlyClippedRatio;
    double maxHighMapqRatio;
}TransitionRequirements;

/*
 * Construct a TransitionRequirements structure
 */
TransitionRequirements *TransitionRequirements_construct(double minHighlyClippedRatio, double maxHighMapqRatio);

/*
 * Destruct a TransitionRequirements structure
 */
void TransitionRequirements_destruct(TransitionRequirements *transitionRequirements);

/*
 * Definition for a validity function that returns true if the state is valid given a CoverageInfo
 * and a TransitionRequirements structure
 */
typedef bool (*ValidityFunction)(StateType , CoverageInfo*, TransitionRequirements *requirements);

/*
 * A validity function that returns true if it is valid to have a transition to the duplicated state otherwise false
 * It returns true if state is anything other than duplicated
 */
bool ValidityFunction_checkDupByMapq(StateType state, CoverageInfo *coverageInfo, TransitionRequirements *requirements);

/*
 * A validity function that returns true if it is valid to have a transition to the misjoin state otherwise false
 * It returns true if the given state is anything other than misjoin
 */
bool ValidityFunction_checkMsjByClipping(StateType state, CoverageInfo *coverageInfo, TransitionRequirements *requirements);

/*! @typedef
 * @abstract Structure for holding the counts that are used for estimating transition probabilities
 * @field countMatrix           A double matrix containing the observed counts for each transition
 *                              (also include start and end state)
 * @field pseudoCountMatrix     A double matrix containing the pseudo counts for each transition
 *                              (also include start and end state)
 * @field numberOfStates        Number of states (note that the matrices have the dimension of numberOfStats + 1)
 *
 */
typedef struct TransitionCountData{
    MatrixDouble *countMatrix; //[numberOfStates + 1] x [numberOfStates + 1]
    MatrixDouble *pseudoCountMatrix; //[numberOfStates + 1] x [numberOfStates + 1]
    int numberOfStates;
    // TransitionCountData can be updated with observations
    // in any arbitrary order
    // since TransitionCountData might be shared between
    // multiple threads for faster runtime it is necessary
    // to keep a mutex
    pthread_mutex_t *mutexPtr;
} TransitionCountData;

/*
 * Construct a TransitionCountData structure
 */
TransitionCountData *TransitionCountData_construct(int numberOfStates);

/*
 * Parse a pseudo count matrix from a file
 * dim should be equal to numberOfStates + 1
 */
void TransitionCountData_parsePseudoCountFromFile(TransitionCountData *transitionCountData, char* pathToMatrix, int dim);

void TransitionCountData_setPseudoCountMatrix(TransitionCountData *transitionCountData, double value);

/*
 * Increment counts data for transition from preState to state
 */
void TransitionCountData_increment(TransitionCountData *transitionCountData, double count, StateType preState, StateType state);

/*
 * Set all elements in countMatrix to zero
 */
void TransitionCountData_resetCountMatrix(TransitionCountData *transitionCountData);

/*
 * Destruct a TransitionCountData structure
 */
TransitionCountData *TransitionCountData_destruct(TransitionCountData *transitionCountData);

/*! @typedef
 * @abstract Structure for holding the current transition probabilities, additional attributes for
 *           conditional transitions and estimating transition probabilities
 * @field matrix                        A double matrix containing the current transition probabilities
 * @field transitionCountData           A TransitionCountData structure keeping counts for estimating
 *                                      transition probabilities
 * @field numberOfStates                Number of states
 * @field validityFunctions             An array of validity functions that has to be checked for conditional transition
 * @field numberOfValidityFunctions     Number of validity functions
 * @field requirements                  A TransitionRequirements structure which is necessary for calling
 *                                      validity functions
 */
typedef struct Transition{
    MatrixDouble *matrix;
    TransitionCountData *transitionCountData;
    int numberOfStates;
    ValidityFunction *validityFunctions;
    int numberOfValidityFunctions;
    TransitionRequirements *requirements;
    double terminationProb;
}Transition;

/*
 * Construct a Transition structure with uniform probabilities
 */
Transition *Transition_constructUniform(int numberOfStates);

/*
 * Construct a Transition structure with a symmetric matrix. The diagonal probability is the probability
 * of staying in a state. The non-diagonal values are set uniformly to make the sum of each row equal to 1
 */
Transition *Transition_constructSymmetricBiased(int numberOfStates, double diagonalProb);

/*
 * Add a TransitionRequirements structure, which will be used by validity functions
 */
void Transition_addRequirements(Transition *transition, TransitionRequirements *requirements);

/*
 * Add a TransitionRequirements structure, which will be used by validity functions
 */
void Transition_addValidityFunction(Transition *transition, ValidityFunction validityFunction);

bool Transition_estimateTransitionMatrix(Transition *transition, double convergenceTol);

void Transition_resetCountData(Transition *transition);

/*
 * Destruct a Transition structure with uniform probabilities
 */
void Transition_destruct(Transition *transition);

/*
 * Get the probability of transitioning from preState to state
 */
double Transition_getProb(Transition *transition, StateType preState, StateType state);

/*
 * Get the conditional probability of transitioning from preState to state by applying the
 * previously added validity functions
 */
double Transition_getProbConditional(Transition *transition, StateType preState, StateType state, CoverageInfo *coverageInfo);

/*
 * Get the probability of ending with state
 */
double Transition_getTerminationProb(Transition *transition, StateType state);

/*
 * Get the probability of starting with state
 */
double Transition_getStartProb(Transition *transition, StateType state);

/*
 * Check if it is valid to be in state using the previously added validity functions
 */
bool Transition_isStateValid(Transition *transition, StateType state, CoverageInfo *coverageInfo);

#endif //HMM_UTILS_H
