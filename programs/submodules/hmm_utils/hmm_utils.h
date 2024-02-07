#ifndef HMM_UTILS_H
#define HMM_UTILS_H

#include <stdint.h>
#include "data_types.h"
#include "ptBlock.h"

#define EXP_TRUNC_POINT_COV_FRACTION 0.2
#define NUMBER_OF_STATES 5
#define ERR_COMP_BINDING_COEF 0.1
#define MAX_COVERAGE_VALUE 250
#define START_STATE_INDEX 5
#define END_STATE_INDEX 5

typedef enum StateType {
    STATE_ERR, STATE_DUP, STATE_HAP, STATE_MSJ, STATE_COL
} StateType;

typedef enum DependencyType {
    ALPHA_ERR, ALPHA_DUP, ALPHA_HAP, ALPHA_MSJ, ALPHA_COL, ALPHA_TRN
} DependencyType;

typedef enum DistType {
    DIST_TRUNC_EXPONENTIAL, DIST_GAUSSIAN, DIST_NEGATIVE_BINOMIAL
}DistType;

typedef enum ModelType {
    MODEL_TRUNC_EXP_GAUSSIAN, MODEL_GAUSSIAN, MODEL_NEGATIVE_BINOMIAL
} ModelType;

typedef enum GaussianParameterType{
    GAUSSIAN_MEAN, GAUSSIAN_VAR
} GaussianParameterType;

typedef enum NegativeBinomialParameterType{
    NB_MEAN, NB_VAR
} NegativeBinomialParameterType;

typedef enum TruncatedExponentialParameterType{
    TRUNC_EXP_MEAN
} TruncatedExponentialParameterType;

typedef struct ParameterBinding{
    double **coefs; // [numberOfParams] x [numberOfComps]
    int numberOfComps;
    int numberOfParams;
} ParameterBinding;

ParameterBinding *ParameterBinding_construct(double **coefs, int numberOfParams, int numberOfComps);
ParameterBinding *ParameterBinding_constructSequenceByStep(double* firstCoefsPerParam,double* stepsPerParam,
                                                           int numberOfParams,int numberOfComps);
ParameterBinding *ParameterBinding_constructForSingleComp(double* coefsPerParam, int numberOfParams);

double ParameterBinding_getValue(ParameterBinding* parameterBinding, int paramIndex, int compIndex);
ParameterBinding **ParameterBinding_getDefault1DArrayForModel(ModelType modelType, int numberOfCollapsedComps);
ParameterBinding **ParameterBinding_getDefault1DArrayForGaussian(int numberOfCollapsedComps);
ParameterBinding **ParameterBinding_getDefault1DArrayForNegativeBinomial(int numberOfCollapsedComps);
ParameterBinding **ParameterBinding_getDefault1DArrayForTruncExpGaussian(int numberOfCollapsedComps);
void ParameterBinding_destruct1DArray(ParameterBinding **parameterBinding1DArray, int arrayLength);



typedef struct CountData{
    double *counts;
    int countsLength;
    int totalCount;
    double sum;
}CountData;

CountData *CountData_construct(uint8_t length);
CountData **CountData_construct1DArray(uint8_t length, int arrayLength);
void CountData_destruct1DArray(CountData **countData1DArray, int arrayLength);
void CountData_increment(CountData * countData, uint8_t value, double count);
void CountData_reset(CountData* countData);
int CountData_getMostFrequentValue(CountData * countData, int minValue, int maxValue);
int CountData_getMaxCount(CountData *countData, int minValue, int maxValue);
void CountData_destruct(CountData* countData);



typedef struct NegativeBinomial {
    double *mean; // mean
    double *var; // variance
    double *theta; // theta
    double *r; // r
    double *weights; // weights of components
    int numberOfComps; // number of mixture components
} NegativeBinomial;

NegativeBinomial *NegativeBinomial_construct(double *mean, double *var, int numberOfComps);
NegativeBinomial *NegativeBinomial_constructByMean(double *mean, double factor, int numberOfComps);
double NegativeBinomial_getTheta(double mean, double var);
double NegativeBinomial_getR(double mean, double var);
double NegativeBinomial_getMean(double theta, double r);
double NegativeBinomial_getVar(double theta, double r);
double NegativeBinomial_getProb(NegativeBinomial *nb, uint8_t x);
double *NegativeBinomial_getComponentProbs(NegativeBinomial *nb, uint8_t x);
void NegativeBinomial_destruct(NegativeBinomial *nb);



typedef struct Gaussian {
    double *mean;
    double *var;
    double *weights; // weights of components
    int numberOfComps; // number of mixture components
} Gaussian;

Gaussian *Gaussian_construct(double *mean, double *var, int numberOfComps);
Gaussian *Gaussian_constructByMean(double *mean, double factor, int numberOfComps);
double Gaussian_getProb(Gaussian *gaussian, uint8_t x, uint8_t preX, double alpha);
double *Gaussian_getComponentProbs(Gaussian *gaussian, uint8_t x, uint8_t preX, double alpha);
void Gaussian_destruct(Gaussian *gaussian);



typedef struct TruncExponential {
    double lambda; // lambda or rate = 1/mean
    double truncPoint; // truncation point
} TruncExponential;

TruncExponential * TruncExponential_construct(double lambda, double truncPoint);
double TruncExponential_getProb(TruncExponential* truncExponential, uint8_t x);
double TruncExponential_getLogLikelihood(TruncExponential *truncExponential, CountData *countData);
double TruncExponential_getLogLikelihoodByParams(double lambda, double truncPoint, CountData *countData);
double TruncExponential_estimateLambda(TruncExponential *truncExponential, CountData *countData, double tol);
void TruncExponential_destruct(TruncExponential* truncExponential);


typedef struct EmissionDist {
    void *dist;
    DistType distType;
}EmissionDist;

EmissionDist *EmissionDist_construct(void* dist, DistType distType);
double EmissionDist_getProb(EmissionDist *emissionDist, uint8_t x, uint8_t preX, double alpha);
void EmissionDist_destruct(EmissionDist* emissionDist);

typedef struct EmissionDistSeries {
    EmissionDist **emissionDists;
    CountData **countDataPerDist;
    ParameterBinding **parameterBindingPerDist;
    int numberOfDists;
    ModelType modelType;
}EmissionDistSeries;

EmissionDistSeries *EmissionDistSeries_constructForModel(ModelType modelType,
                                                         double **means,  // [numberOfDists] x [maxMixtures]
                                                         int *numberOfCompsPerDist,
                                                         int numberOfDists);
EmissionDist *EmissionDistSeries_getEmissionDist(EmissionDistSeries* emissionDistSeries, int distIndex);
void EmissionDistSeries_destruct(EmissionDistSeries* emissionDistSeries);
double *EmissionDistSeries_getAllProbs(EmissionDistSeries* emissionDistSeries,
                                       uint8_t x,
                                       uint8_t preX,
                                       double alpha);
double EmissionDistSeries_getProb(EmissionDistSeries* emissionDistSeries,
                                  int distIndex,
                                  uint8_t x,
                                  uint8_t preX,
                                  double alpha);


typedef struct TransitionRequirements{
    double minHighlyClippedRatio;
    double maxHighMapqRatio;
}TransitionRequirements;

TransitionRequirements *TransitionRequirements_construct(double minHighlyClippedRatio, double maxHighMapqRatio);
void TransitionRequirements_destruct(TransitionRequirements *transitionRequirements);

typedef bool (*ValidityFunction)(StateType , CoverageInfo*, TransitionRequirements *requirements);

bool ValidityFunction_checkDupByMapq(StateType state, CoverageInfo *coverageInfo, TransitionRequirements *requirements);
bool ValidityFunction_checkMsjByClipping(StateType state, CoverageInfo *coverageInfo, TransitionRequirements *requirements);

typedef struct TransitionCountData{
    MatrixDouble *countMatrix;
    MatrixDouble *pseudoCountMatrix;
    int numberOfStates;
} TransitionCountData;

TransitionCountData *TransitionCountData_construct(int numberOfStates);
void TransitionCountData_parsePseudoCountFromFile(TransitionCountData *transitionCountData, char* pathToMatrix, int dim);
TransitionCountData *TransitionCountData_destruct(TransitionCountData *transitionCountData);


typedef struct Transition{
    MatrixDouble *matrix;
    TransitionCountData *transitionCountData;
    int numberOfStates;
    ValidityFunction *validityFunctions;
    int numberOfValidityFunctions;
    TransitionRequirements *requirements;
}Transition;

Transition *Transition_constructUniform(int numberOfStates);
Transition *Transition_constructSymmetricBiased(int numberOfStates, double diagonalProb);
void Transition_addRequirements(Transition *transition, TransitionRequirements *requirements);
void Transition_addValidityFunction(Transition *transition, ValidityFunction validityFunction);
void Transition_destruct(Transition *transition);
double Transition_getProb(Transition *transition, StateType preState, StateType state);
double Transition_getProbConditional(Transition *transition, StateType preState, StateType state, CoverageInfo *coverageInfo);
double Transition_getTerminationProb(Transition *transition, StateType state);
double Transition_getStartProb(Transition *transition, StateType state);
bool Transition_isStateValid(Transition *transition, StateType state, CoverageInfo *coverageInfo);

#endif //HMM_UTILS_H
