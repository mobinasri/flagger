#include <math.h>
#include <assert.h>
#include "data_types.h"
#include "hmm.h"
#include <stdio.h>
#include "block_it.h"
#include "common.h"
#include "sonLib.h"

#define PI 3.14159

static pthread_mutex_t chunkMutex;


/**
 * Constructs a Negative Binomial object containing the mean vectors, covariance
 * matrices and the weights for all of its mixture components
 *
 * @param mu	An array of double vectors containing the initial mean vectors
 * @param cov	An array of double matrices containing the initial covariance matrices
 * @param n	The number of mixture components
 * @return gaussian
 */
NegativeBinomial *NegativeBinomial_construct(VectorDouble **mu, MatrixDouble **cov, int n) {
    NegativeBinomial *nb = malloc(1 * sizeof(NegativeBinomial));
    // Mixture Gaussian Parameters:
    //
    // Allocating arrays for NB parameters
    nb->mu = VectorDouble_constructArray1D(n, mu[0]->dim);
    nb->cov = MatrixDouble_constructArray1D(n, mu[0]->dim, mu[0]->dim);
    // Allocating and initializing mixture weights
    nb->weights = malloc(n * sizeof(double));
    nb->n = n;
    for (int m = 0; m < n; m++) {
        // init mean vector
        VectorDouble_copyInPlace(nb->mu[m], mu[m]);
        // init covariance matrix
        MatrixDouble_copyInPlace(nb->cov[m], cov[m]);
        // init weights uniformly
        nb->weights[m] = 1.0 / n;
    }
    // Estimation Statistics:
    //
    // Allocating and initializing arrays for estimating parameters
    nb->thetaNum = VectorDouble_constructArray1D(n, mu[0]->dim);
    nb->thetaDenom = VectorDouble_constructArray1D(n, mu[0]->dim);
    nb->lambdaNum = VectorDouble_constructArray1D(n, mu[0]->dim);
    nb->lambdaDenom = VectorDouble_constructArray1D(n, mu[0]->dim);
    nb->weightNum = (double *) malloc(n * sizeof(double));
    memset(nb->weightNum, 0, n * sizeof(double));
    return nb;
}

NegativeBinomial * NegativeBinomial_constructSufficientStats(int nEmit, int nMixtures) {
    NegativeBinomial *nb = malloc(1 * sizeof(NegativeBinomial));
    // Allocating arrays for gaussian parameters
    nb->mu = VectorDouble_constructArray1D(nMixtures, nEmit);
    nb->cov = MatrixDouble_constructArray1D(nMixtures, nEmit, nEmit);
    // Allocating and initializing mixture weights
    nb->weights = malloc(nMixtures * sizeof(double));
    nb->n = nMixtures;
    // Estimation Statistics:
    //
    // Allocating and initializing arrays for estimating parameters
    nb->thetaNum = VectorDouble_constructArray1D(nMixtures, nEmit);
    nb->thetaDenom = VectorDouble_constructArray1D(nMixtures, nEmit);
    nb->lambdaNum = VectorDouble_constructArray1D(nMixtures, nEmit);
    nb->lambdaDenom = VectorDouble_constructArray1D(nMixtures, nEmit);
    nb->weightNum = (double *) malloc(nMixtures * sizeof(double));
    memset(nb->weightNum, 0, nMixtures * sizeof(double));
    return nb;
}

/**
 * Destructs a Gaussian object
 *
 * @param gaussian
 */
void NegativeBinomial_destruct(NegativeBinomial *nb) {
    VectorDouble_destructArray1D(nb->mu, nb->n);
    VectorDouble_destructArray1D(nb->thetaNum, nb->n);
    VectorDouble_destructArray1D(nb->thetaDenom, nb->n);
    MatrixDouble_destructArray1D(nb->cov, nb->n);
    VectorDouble_destructArray1D(nb->lambdaNum, nb->n);
    VectorDouble_destructArray1D(nb->lambdaDenom, nb->n);
    free(nb->weights);
    free(nb->weightNum);
    free(nb);
}

/**
 * Constructs a NB object for which the diagonal covariance entries
 * are initialized to the initial mean values
 *
 * @param mu	An array of double vectors containing the initial mean vectors
 * @param n	The number of mixture components
 * @return gaussian
 */
NegativeBinomial *NegativeBinomial_constructSpecial(VectorDouble **mu, int n) {
    // Allocate and initialize a special covariance matrix
    // It would be a diagonal matrix with autovariances same as means
    MatrixDouble **cov = MatrixDouble_constructArray1D(n, mu[0]->dim, mu[0]->dim);
    for (int m = 0; m < n; m++) {
        for (int j = 0; j < mu[m]->dim; j++) {
            cov[m]->data[j][j] = mu[m]->data[j] * 1.1;
        }
    }
    return NegativeBinomial_construct(mu, cov, n);
}


double NegativeBinomial_getTheta(double mean, double var){
    double theta = mean / var;
    return theta;
}

double NegativeBinomial_getR(double mean, double var){
    double r = pow(mean, 2) / (var - mean);
    return r;
}


double NegativeBinomial_getMean(double theta, double r){
    double mean = r * (1 - theta) / theta;
    return mean;
}

double NegativeBinomial_getVar(double theta, double r){
    double var = r * (1 - theta) / pow(theta, 2);
    return var;
}

/**
 * Returns the total probability of emitting a vector from the given Negative Binomial
 *
 * @param vec	The emitted vector
 * @param nb	The NegativeBinomial object
 * @return prob	The emission probability
 */

double NegativeBinomial_getProb(VectorChar *vec, NegativeBinomial *nb, int comp) {
    double *probs = NegativeBinomial_getMixtureProbs(vec, nb, comp);
    double totProb = 0.0;
    for (int m = 0; m < nb->n; m++) {
        totProb += probs[m];
    }
    free(probs);
    return totProb;
}


/**
 * Returns the probabilities of emitting a vector from the mixture components of Negative Binomial
 *
 * @param vec   The emitted vector
 * @param nb      The Negative Binomial object
 * @return prob An array of probabilities for all Gaussian components
 */

// TODO: Support higher dimension if needed. This function only supports dimensions of 1
double *NegativeBinomial_getMixtureProbs(VectorChar *vec, NegativeBinomial *nb, int comp) {
    uint8_t *x = vec->data;
    int dim = 1;//vec->dim;
    double w;
    double *mu;
    double **c;
    double *probs = malloc(nb->n * sizeof(double));
    memset(probs, 0.0, nb->n * sizeof(double));
    // iterate over mixture components
    for (int m = 0; m < nb->n; m++) {
        mu = nb->mu[m]->data;
        c = nb->cov[m]->data;
        w = nb->weights[m];
        if (dim == 1) {
            double theta = NegativeBinomial_getTheta(mu[0], c[0][0]);
            double r = NegativeBinomial_getR(mu[0],c[0][0]);
            probs[m] = w * exp(lgamma(r + x[0]) - lgamma(r) - lgamma(x[0] + 1) + r * log(theta) + (double) x[0] * log(1-theta));
        }
        if (probs[m] != probs[m]) {
            fprintf(stderr, "prob is NAN\n");
            fprintf(stderr, "%s\n", MatrixDouble_toString(nb->cov[m]));
            fprintf(stderr, "%s\n", VectorDouble_toString(nb->mu[m]));
            exit(EXIT_FAILURE);
        }
        if (probs[m] < 1e-40) {
            //fprintf(stderr, "COMP=%d, prob is lower than 1e-40\n", comp);
            probs[m] = 1e-40;
        }
    }
    return probs;
}

/**
 * Constructs a Gaussian object containing the mean vectors, covariance
 * matrices and the weights for all of its mixture components
 *
 * @param mu	An array of double vectors containing the initial mean vectors
 * @param cov	An array of double matrices containing the initial covariance matrices
 * @param n	The number of mixture components 
 * @return gaussian
 */
Gaussian *Gaussian_construct(VectorDouble **mu, MatrixDouble **cov, int n) {
    Gaussian *gaussian = malloc(1 * sizeof(Gaussian));
    // Mixture Gaussian Parameters:
    //
    // Allocating arrays for gaussian parameters
    gaussian->mu = VectorDouble_constructArray1D(n, mu[0]->dim);
    gaussian->cov = MatrixDouble_constructArray1D(n, mu[0]->dim, mu[0]->dim);
    // Allocating and initializing mixture weights
    gaussian->weights = malloc(n * sizeof(double));
    gaussian->n = n;
    for (int m = 0; m < n; m++) {
        // init mean vector
        VectorDouble_copyInPlace(gaussian->mu[m], mu[m]);
        // init covariance matrix
        MatrixDouble_copyInPlace(gaussian->cov[m], cov[m]);
        // init weights uniformly
        gaussian->weights[m] = 1.0 / n;
    }
    // Estimation Statistics:
    //
    // Allocating and initializing arrays for estimating parameters
    gaussian->muNum = VectorDouble_constructArray1D(n, mu[0]->dim);
    gaussian->muDenom = VectorDouble_constructArray1D(n, mu[0]->dim);
    gaussian->covNum = MatrixDouble_constructArray1D(n, mu[0]->dim, mu[0]->dim);
    gaussian->covDenom = MatrixDouble_constructArray1D(n, mu[0]->dim, mu[0]->dim);
    gaussian->weightNum = (double *) malloc(n * sizeof(double));
    memset(gaussian->weightNum, 0, n * sizeof(double));
    return gaussian;
}


Gaussian *Gaussian_constructSufficientStats(int nEmit, int nMixtures) {
    Gaussian *gaussian = malloc(1 * sizeof(Gaussian));
    // Allocating arrays for gaussian parameters
    gaussian->mu = VectorDouble_constructArray1D(nMixtures, nEmit);
    gaussian->cov = MatrixDouble_constructArray1D(nMixtures, nEmit, nEmit);
    // Allocating and initializing mixture weights
    gaussian->weights = malloc(nMixtures * sizeof(double));
    gaussian->n = nMixtures;
    // Estimation Statistics:
    //
    // Allocating and initializing arrays for estimating parameters
    gaussian->muNum = VectorDouble_constructArray1D(nMixtures, nEmit);
    gaussian->muDenom = VectorDouble_constructArray1D(nMixtures, nEmit);
    gaussian->covNum = MatrixDouble_constructArray1D(nMixtures, nEmit, nEmit);
    gaussian->covDenom = MatrixDouble_constructArray1D(nMixtures, nEmit, nEmit);
    gaussian->weightNum = (double *) malloc(nMixtures * sizeof(double));
    memset(gaussian->weightNum, 0, nMixtures * sizeof(double));
    return gaussian;
}

/**
 * Destructs a Gaussian object
 *
 * @param gaussian
 */
void Gaussian_destruct(Gaussian *gaussian) {
    VectorDouble_destructArray1D(gaussian->mu, gaussian->n);
    VectorDouble_destructArray1D(gaussian->muNum, gaussian->n);
    VectorDouble_destructArray1D(gaussian->muDenom, gaussian->n);
    MatrixDouble_destructArray1D(gaussian->cov, gaussian->n);
    MatrixDouble_destructArray1D(gaussian->covNum, gaussian->n);
    MatrixDouble_destructArray1D(gaussian->covDenom, gaussian->n);
    free(gaussian->weights);
    free(gaussian->weightNum);
    free(gaussian);
}

/** 
 * Constructs a gaussian object for which the diagonal covariance entries 
 * are initialized to the initial mean values
 *
 * @param mu	An array of double vectors containing the initial mean vectors
 * @param n	The number of mixture components
 * @return gaussian
 */
Gaussian *Gaussian_constructSpecial(VectorDouble **mu, int n) {
    // Allocate and initialize a special covariance matrix
    // It would be a diagonal matrix with autovariances same as means
    MatrixDouble **cov = MatrixDouble_constructArray1D(n, mu[0]->dim, mu[0]->dim);
    for (int m = 0; m < n; m++) {
        for (int j = 0; j < mu[m]->dim; j++) {
            cov[m]->data[j][j] = mu[m]->data[j];
        }
    }
    return Gaussian_construct(mu, cov, n);
}


/**
 * Makes a transition matrix with uniform probabilities
 *
 * @param dim	The dimension of the transition matrix supposedly 
 * 		should be equal to the number of states + 1
 * @return trans The transition matrix
 */
MatrixDouble *makeUniformTransition(int dim) {
    MatrixDouble *trans = MatrixDouble_construct0(dim, dim);
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            trans->data[i][j] = 1.0 / dim;
        }
    }
    return trans;
}

/**
 * Makes a transition matrix with probabilites biased toward
 * not changing states
 *
 * @param dim	The dimension of the transition matrix supposedly 
 *              should be equal to the number of states + 1
 * @return trans The transition matrix 
 */
MatrixDouble *makeBiasedTransition(int dim) {
    MatrixDouble *trans = MatrixDouble_construct0(dim, dim);
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            trans->data[i][j] = 1.0 / dim * 1.0e-7;
        }
        // the diagonal entries have higher probs
        trans->data[i][i] = 1.0 - 1.0e-7 * (dim - 1.0) / dim;
    }
    trans->data[2][3] = 1.0e-3;
    trans->data[3][2] = 1.0e-3;
    trans->data[2][2] -= (1.0e-3 - 1.0 / dim * 1.0e-7);
    trans->data[3][3] -= (1.0e-3 - 1.0 / dim * 1.0e-7);
    return trans;
}


/**
 * Returns the total probability of emitting a vector from the given Gaussian
 *
 * @param vec	The emitted vector
 * @param gaussian	The gaussian object
 * @return prob	The emission probability
 */

double getGaussianProb(VectorChar *vec, Gaussian *gaussian, int comp) {
    double *probs = getGaussianMixtureProbs(vec, gaussian, comp);
    double totProb = 0.0;
    for (int m = 0; m < gaussian->n; m++) {
        totProb += probs[m];
    }
    free(probs);
    return totProb;
}


/**
 * Returns the probabilties of emitting a vector from the mixture components of the given Gaussian
 *
 * @param vec   The emitted vector
 * @param gaussian      The gaussian object
 * @return prob An array of probabilities for all Gaussian components
 */

// TODO: Support higher dimension if needed. This function only supports dimensions less than or equal to 2.
double *getGaussianMixtureProbs(VectorChar *vec, Gaussian *gaussian, int comp) {
    uint8_t *x = vec->data;
    int dim = 1;//vec->dim;
    double w;
    double *mu;
    MatrixDouble *covCopy;
    double **c;
    double *probs = malloc(gaussian->n * sizeof(double));
    memset(probs, 0, gaussian->n * sizeof(double));
    double det;
    double a0;
    double a1;
    // iterate over mixture components
    for (int m = 0; m < gaussian->n; m++) {
        mu = gaussian->mu[m]->data;
        covCopy = MatrixDouble_copy(gaussian->cov[m]);
        c = covCopy->data;
        w = gaussian->weights[m];
        if (dim == 1) {
            probs[m] = w / (sqrt(c[0][0] * 2 * PI)) * exp(-0.5 * pow((x[0] - mu[0]), 2) / c[0][0]);
        } else if (dim == 2) {
            det = c[1][1] * c[0][0] - c[1][0] * c[0][1];
            // TODO: Find a better solution for zero determinant
            if (det < 1e-4) {
                //fprintf(stderr, "COMP=%d, The determinant is very low! %.3e\n", comp, det);
                c[1][1] += c[1][1] < 0.01 ? 0.01 : 0;
                c[0][0] += c[0][0] < 0.01 ? 0.01 : 0;
                det = c[1][1] * c[0][0] - c[1][0] * c[0][1];
                //fprintf(stderr, "The determinant changed! %.3e\n", det);
            }
            if (det < 0) {
                fprintf(stderr, "The determinant has a negative value! %.3e\n", det);
                exit(EXIT_FAILURE);
            }
            a0 = (x[0] - mu[0]) * c[1][1] - (x[1] - mu[1]) * c[1][0];
            a1 = (x[1] - mu[1]) * c[0][0] - (x[0] - mu[0]) * c[0][1];
            probs[m] = w / (sqrt(2.0 * PI) * sqrt(det)) * exp(-0.5 / det * ((x[0] - mu[0]) * a0 + (x[1] - mu[1]) * a1));
        }
        MatrixDouble_destruct(covCopy);
        if (probs[m] != probs[m]) {
            double u = -0.5 / det * ((x[0] - mu[0]) * a0 + (x[1] - mu[1]) * a1);
            fprintf(stderr, "prob is NAN exp(%.2e) det=%.2e, w=%.2e, exp()=%.2e!\n", u, det, w, exp(u));
            fprintf(stderr, "%s\n", MatrixDouble_toString(gaussian->cov[m]));
            fprintf(stderr, "%s\n", VectorDouble_toString(gaussian->mu[m]));
            exit(EXIT_FAILURE);
        }
        if (probs[m] < 1e-40) {
            //fprintf(stderr, "COMP=%d, prob is lower than 1e-40\n", comp);
            probs[m] = 1e-40;
        }
    }
    /*if(prob < 1e-100){
        double u = -0.5 / det * ((x[0] - mu[0]) * a0 + (x[1] - mu[1]) * a1);
        fprintf(stderr, "#%s\n", VectorDouble_toString(gaussian->mu));
        VectorDouble* muVect = gaussian->mu;
        fprintf(stderr, "x=[%d %d], mu = [%.2f %.2f], det=%.2e, u=%.2e\n", x[0], x[1], muVect->data[0], muVect->data[1], det, u);
    }*/
    return probs;
}

double getExpProb(double x, double lambda) {
    double prob = x < 0 ? 0 : lambda * exp(-1 * lambda * x);
    return prob;
}

/**
 * Constructs an HMM object 
 *
 * @param nClasses	The number of region classes like HSats-1,2,and 3
 * @param nComps	The number of components (states) for each region class
 * @param nEmit		The dimension of each emitted vector
 * @param nMixtures	An array that contains the number of mixture components for each Gaussian
 * 			nMixrures[HMM components]
 * @param mu		mu[classes][HMM components][mixture components][emission dim]
 * 			mu[r][c] ---> [mixture components][emission dim]
 * @return hmm		The HMM object
 */

HMM *HMM_construct(int nClasses, int nComps, int nEmit, int *nMixtures, VectorDouble ****mu, VectorDouble ***muFactors,
                   MatrixDouble ***covFactors, double maxHighMapqRatio, MatrixDouble **transNum,
                   MatrixDouble **transDenom, ModelType modelType) {
    HMM *model = malloc(sizeof(HMM));
    model->modelType = modelType;
    model->nClasses = nClasses;
    model->nComps = nComps;
    model->nEmit = nEmit;
    model->nMixtures = nMixtures;
    model->muFactors = muFactors;
    model->covFactors = covFactors;
    model->maxHighMapqRatio = maxHighMapqRatio;
    model->terminateProb = 1e-2;
    if (modelType == GAUSSIAN) {
        // Constructing the emission Gaussians and set their parameters
        // emit[r][c] is pointing to the Gaussian of the c-th component of the r-th class
        model->emit = (Gaussian ***) malloc(nClasses * sizeof(Gaussian **));
        for (int r = 0; r < nClasses; r++) {
            model->emit[r] = (Gaussian **) malloc(nComps * sizeof(Gaussian *));
            // Initialize the emission parameters of each state with the given vector
            for (int c = 0; c < nComps; c++) {
                model->emit[r][c] = Gaussian_constructSpecial(mu[r][c], nMixtures[c]);
            }
        }
    }
    else if (modelType == NEGATIVE_BINOMIAL){
        model->emit = (NegativeBinomial ***) malloc(nClasses * sizeof(NegativeBinomial **));
        for (int r = 0; r < nClasses; r++) {
            model->emit[r] = (NegativeBinomial **) malloc(nComps * sizeof(NegativeBinomial *));
            // Initialize the emission parameters of each state with the given vector
            for (int c = 0; c < nComps; c++) {
                model->emit[r][c] = NegativeBinomial_constructSpecial(mu[r][c], nMixtures[c]);
            }
        }
    }

    model->transNum = transNum;
    model->transDenom = transDenom;

    // The last row of the transition matrix is pointing to the probabilities of
    // starting the sequence from each possible state
    // The last column of the transition matrix is pointing to the probabilities of
    // ending the sequence with each possible state
    // The dimension of the transition matrix is (#nComps + 1) x (#nComps + 1)
    // Each class of region will have its own transition matrix
    model->trans = (MatrixDouble **) malloc(nClasses * sizeof(MatrixDouble * ));
    for (int r = 0; r < nClasses; r++) {
        model->trans[r] = MatrixDouble_construct0(nComps + 1, nComps + 1);
    }
    for (int r = 0; r < nClasses; r++) {
        for (int c1 = 0; c1 < nComps; c1++) {
            for (int c2 = 0; c2 < nComps; c2++) {
                model->trans[r]->data[c1][c2] = (1 - model->terminateProb) * model->transNum[r]->data[c1][c2] /
                                                model->transDenom[r]->data[c1][c2];
            }
            model->trans[r]->data[c1][nComps] = model->terminateProb;
        }
        for (int c2 = 0; c2 < nComps; c2++) {
            model->trans[r]->data[nComps][c2] = (1 - model->terminateProb) * 0.25;
        }
        model->trans[r]->data[nComps][nComps] = model->terminateProb;
    }
    fprintf(stderr, "%s\n", MatrixDouble_toString(model->trans[0]));

    // Initialize sufficient stats for estimating transition probs
    model->transCounts = (MatrixDouble **) malloc(nClasses * sizeof(MatrixDouble * ));
    for (int r = 0; r < nClasses; r++) {
        model->transCounts[r] = MatrixDouble_construct0(nComps + 1, nComps + 1);
    }

    model->mutexPtr = malloc(sizeof(pthread_mutex_t));
    pthread_mutex_init(model->mutexPtr, NULL);
    return model;
}

/**
 * Destructs an HMM object
 *
 * @param model	an HMM object
 */
void HMM_destruct(HMM *model) {

    int nClasses = model->nClasses;
    int nComps = model->nComps;
    int nEmit = model->nEmit;

    //free emit
    for (int r = 0; r < nClasses; r++) {
        for (int c = 0; c < nComps; c++) {
            if(model->modelType == GAUSSIAN) {
                Gaussian_destruct(model->emit[r][c]);
            }
            else if(model->modelType == NEGATIVE_BINOMIAL){
                NegativeBinomial_destruct(model->emit[r][c]);
            }
        }
        free(model->emit[r]);
    }
    free(model->emit);

    //free transition matrices
    for (int r = 0; r < nClasses; r++) {
        MatrixDouble_destruct(model->trans[r]);
    }
    free(model->trans);

    // free sufficient stats for transition matrix
    for (int r = 0; r < nClasses; r++) {
        MatrixDouble_destruct(model->transCounts[r]);
    }
    free(model->transCounts);
    free(model->mutexPtr);

    VectorDouble_destructArray2D(model->muFactors, model->nComps, 4);
    MatrixDouble_destructArray2D(model->covFactors, model->nComps, 4);
    free(model);
}

EM *EM_construct(VectorChar **seqEmit, uint8_t *seqClass, int seqLength, HMM *model) {
    //fprintf(stderr, "Constructing EM\n");
    assert(seqLength > 0);
    EM *em = malloc(sizeof(EM));
    em->nClasses = model->nClasses;
    em->nComps = model->nComps;
    em->nEmit = model->nEmit;
    // Copy the emission sequence
    em->seqEmit = (VectorChar **) malloc(seqLength * sizeof(VectorChar * ));
    for (int i = 0; i < seqLength; i++) {
        em->seqEmit[i] = VectorChar_copy(seqEmit[i]);
    }
    // Copy the classes sequence
    em->seqClass = (uint8_t *) malloc(seqLength * sizeof(uint8_t));
    memcpy(em->seqClass, seqClass, seqLength);
    if(model->modelType == GAUSSIAN) {
        // Construct Gaussian objects for saving sufficient stats
        em->emit = (Gaussian ***) malloc(model->nClasses * sizeof(Gaussian **));
        for (int r = 0; r < model->nClasses; r++) {
            em->emit[r] = (Gaussian **) malloc(model->nComps * sizeof(Gaussian *));
            // Initialize the emission parameters of each state with the given vector
            for (int c = 0; c < model->nComps; c++) {
                Gaussian *gaussian = model->emit[r][c];
                em->emit[r][c] = Gaussian_constructSufficientStats(model->nEmit, gaussian->n);
            }
        }
    }
    else if(model->modelType == NEGATIVE_BINOMIAL){
        em->emit = (NegativeBinomial ***) malloc(model->nClasses * sizeof(NegativeBinomial **));
        for (int r = 0; r < model->nClasses; r++) {
            em->emit[r] = (NegativeBinomial **) malloc(model->nComps * sizeof(NegativeBinomial *));
            // Initialize the emission parameters of each state with the given vector
            for (int c = 0; c < model->nComps; c++) {
                NegativeBinomial *nb = model->emit[r][c];
                em->emit[r][c] = NegativeBinomial_constructSufficientStats(model->nEmit, nb->n);
            }
        }
    }
    // Allocate and initialize forward and backward matrices
    em->seqLength = seqLength;
    em->f = (double **) malloc(seqLength * sizeof(double *));
    em->b = (double **) malloc(seqLength * sizeof(double *));
    for (int i = 0; i < seqLength; i++) {
        em->f[i] = (double *) malloc(model->nComps * sizeof(double));
        em->b[i] = (double *) malloc(model->nComps * sizeof(double));
        for (int c = 0; c < model->nComps; c++) {
            em->f[i][c] = 0.0;
            em->b[i][c] = 0.0;
        }
    }
    // Initialize scale to avoid underflow
    em->scales = (double *) malloc(seqLength * sizeof(double));
    for (int i = 0; i < seqLength; i++) {
        em->scales[i] = 1.0;
    }
    em->px = -1.0;
    em->modelType = model->modelType;
    //fprintf(stderr, "Constructed\n");
    return em;
}

void EM_destruct(EM *em) {
    int nClasses = em->nClasses;
    int nComps = em->nComps;
    int nEmit = em->nEmit;

    int seqLength = em->seqLength;
    for (int i = 0; i < seqLength; i++) {
        VectorChar_destruct(em->seqEmit[i]);
    }
    free(em->seqEmit);
    free(em->seqClass);

    for (int r = 0; r < nClasses; r++) {
        for (int c = 0; c < nComps; c++) {
            if(em->modelType == GAUSSIAN) {
                Gaussian_destruct(em->emit[r][c]);
            }
            else if(em->modelType == NEGATIVE_BINOMIAL){
                NegativeBinomial_destruct(em->emit[r][c]);
            }
        }
        free(em->emit[r]);
    }
    free(em->emit);

    for (int i = 0; i < seqLength; i++) {
        free(em->f[i]);
        free(em->b[i]);
    }
    free(em->f);
    free(em->b);

    free(em->scales);
    free(em);
}

void runForward(HMM *model, EM *em) {
    // Get em sequence attributes
    VectorChar **seqEmit = em->seqEmit; // Sequence of emission vectors (can be vector of length 1)
    uint8_t *seqClass = em->seqClass; // Sequence of region class numbers
    int seqLength = em->seqLength; // Length of sequence
    // Get model attributes
    int nComps = model->nComps; // Number of components (usually 4)
    int nClasses = model->nClasses; // Number of region classes
    int nEmit = 1;//model->nEmit; // The emission dimension
    MatrixDouble **trans = model->trans; // Transition matrices. For each class we have a separate matrix
    //Gaussian ***emit = model->emit; // A 2D array [nClasses x nComps] of Gaussian objects
    double scale = 0.0; // For scaling the forward probabilities at each location to avoid underflow
    // Initialize to zero
    for (int i = 0; i < seqLength; i++) {
        for (int c = 0; c < nComps; c++) {
            em->f[i][c] = 0.0;
        }
    }
    double eProb;
    double tProb;
    // Set the 0-th block of the forward matrix
    scale = 0.0;
    for (int c = 0; c < nComps; c++) {
        // Emission probability
        if(model->modelType == GAUSSIAN) {
            eProb = getGaussianProb(seqEmit[0], model->emit[seqClass[0]][c], c);
        }
        else if(model->modelType == NEGATIVE_BINOMIAL){
            eProb = NegativeBinomial_getProb(seqEmit[0], model->emit[seqClass[0]][c], c);
        }
        // Transition probability
        //
        // trans[]->data[nComps][c] holds the prob of starting with c
        // This ratio is the number of alignments with high mapq against the total alignment in this locus
        double highMapqRatio = (double) seqEmit[0]->data[1] / seqEmit[0]->data[0];
        if (model->maxHighMapqRatio < highMapqRatio) { // Due to high number of high mapq it cannot be false duplication
            // probs are getting renormalized exlcuding transision to false duplication (state index = 1)
            tProb = c == 1 ? 0.0 : trans[seqClass[0]]->data[nComps][c] / (1 - trans[seqClass[0]]->data[nComps][c]);
        } else {// ratio is low so it can be false duplication so all probs are being used
            tProb = trans[seqClass[0]]->data[nComps][c];
        }
        // Update forward
        em->f[0][c] = eProb * tProb;
        scale += em->f[0][c];
    }
    em->scales[0] = scale;
    // Scale f 0-th column
    for (int c = 0; c < nComps; c++) {
        em->f[0][c] /= scale;
    }
    // Fill the remaining rows of the forward matrix
    for (int i = 1; i < seqLength; i++) {
        scale = 0.0;
        for (int c2 = 0; c2 < nComps; c2++) {
            for (int c1 = 0; c1 < nComps; c1++) { // Transition from c1 comp to c2 comp
                if (seqClass[i - 1] != seqClass[i]) { // if the class has changed
                    // Make the transition prob uniform
                    tProb = 1.0 / (nComps + 1);
                } else {
                    // This ratio is the number of alignments with high mapq against the total alignment in this locus
                    double highMapqRatio = (double) seqEmit[i]->data[1] / seqEmit[i]->data[0];
                    if (model->maxHighMapqRatio <
                        highMapqRatio) { // Due to high number of high mapq it cannot be false duplication
                        // probs are getting renormalized exlcuding transision to false duplication (state index = 1)
                        tProb = c2 == 1 ? 0.0 : trans[seqClass[i]]->data[c1][c2] /
                                                (1 - trans[seqClass[i]]->data[c1][1]);
                    } else {// ratio is low so it can be false duplication so all probs are being used
                        tProb = trans[seqClass[i]]->data[c1][c2];
                    }
                }
                em->f[i][c2] += (em->f[i - 1][c1] * tProb);
            }
            if(model->modelType == GAUSSIAN) {
                eProb = getGaussianProb(seqEmit[i], model->emit[seqClass[i]][c2], c2);
            }
            if(model->modelType == NEGATIVE_BINOMIAL) {
                eProb = NegativeBinomial_getProb(seqEmit[i], model->emit[seqClass[i]][c2], c2);
            }
            em->f[i][c2] *= eProb;
            scale += em->f[i][c2];
        }
        if (scale < 1e-100) {
            fprintf(stderr, "scale is very low! %.2e\n", scale);
            /*for (int c2 = 0; c2 < nComps; c2++) {
                eProb = getGaussianProb(seqEmit[i], emit[seqClass[i]][c2], c2);
                //fprintf(stderr, "vec=%s,\n mu=%s\nc2=%d, eprob=%.2e\n", VectorChar_toString(seqEmit[i]), VectorDouble_toString(emit[seqClass[i]][c2]), c2, eProb);
            }*/
            exit(EXIT_FAILURE);
        }
        // Scale f
        for (int c = 0; c < nComps; c++) {
            em->f[i][c] /= scale;
        }
        em->scales[i] = scale;
    }
    // Update P(x)
    // TODO: P(x) is calculated wrongly here p(x) = s(1) x s(2) x ... x s(L) check Durbin's
    // It shouldn't make any problem in the training process since px is not needed in the scaled calculations
    em->px = 0;
    for (int c = 0; c < nComps; c++) {
        tProb = trans[seqClass[seqLength - 1]]->data[c][nComps];
        // data[c][nComps] is the probability of ending with state c
        em->px += em->f[seqLength - 1][c] * tProb;
    }
}


// runForward should be run before runBackward because of scales
// after running runForward scales will be save in em->scales
// so they can be used in backward
void runBackward(HMM *model, EM *em) {
    // Get em attributes
    VectorChar **seqEmit = em->seqEmit;
    uint8_t *seqClass = em->seqClass;
    int seqLength = em->seqLength;
    // Get model attributes
    int nComps = model->nComps;
    int nClasses = model->nClasses;
    int nEmit = 1;//model->nEmit;
    MatrixDouble **trans = model->trans;
    //Gaussian ***emit = model->emit;
    // Initialize to zero
    for (int i = 0; i < seqLength; i++) {
        for (int c = 0; c < nComps; c++) {
            em->b[i][c] = 0.0;
        }
    }
    // Set the last column of the backward matrix
    for (int c = 0; c < nComps; c++) {
        // trans[]->data[c][nComps] holds the probability of ending with c
        em->b[seqLength - 1][c] = trans[seqClass[seqLength - 1]]->data[c][nComps];
        // Scale b
        em->b[seqLength - 1][c] /= em->scales[seqLength - 1];
        if (em->b[seqLength - 1][c] != em->b[seqLength - 1][c]) {
            fprintf(stderr, "b[seqLength -1] is NAN!\n");
            fprintf(stderr, "scale = %.2e\n", em->scales[seqLength - 1]);
            exit(EXIT_FAILURE);
        }
    }
    double tProb;
    double eProb;
    // Fill the remaining columns of the backward matrix
    for (int i = seqLength - 2; 0 <= i; i--) {
        for (int c1 = 0; c1 < nComps; c1++) {
            for (int c2 = 0; c2 < nComps; c2++) { // Transition from c1 t c2
                if (seqClass[i] != seqClass[i + 1]) {// if the class has changed
                    // Make the transition prob uniform
                    tProb = 1.0 / (nComps + 1);
                } else {
                    double highMapqRatio = (double) seqEmit[i + 1]->data[1] / seqEmit[i + 1]->data[0];
                    if (model->maxHighMapqRatio <
                        highMapqRatio) { // Due to high number of high mapq it cannot be false duplication
                        // probs are getting renormalized exlcuding transition to false duplication (state index = 1)
                        tProb = c2 == 1 ? 0.0 : trans[seqClass[i]]->data[c1][c2] /
                                                (1 - trans[seqClass[i]]->data[c1][1]);
                    } else {
                        tProb = trans[seqClass[i]]->data[c1][c2];
                    }
                }
                if(model->modelType == GAUSSIAN) {
                    eProb = getGaussianProb(seqEmit[i + 1], model->emit[seqClass[i + 1]][c2], c2);
                }
                else if(model->modelType == NEGATIVE_BINOMIAL) {
                    eProb = NegativeBinomial_getProb(seqEmit[i + 1], model->emit[seqClass[i + 1]][c2], c2);
                }
                em->b[i][c1] += tProb * eProb * em->b[i + 1][c2];
            }
            // Scale b
            em->b[i][c1] /= em->scales[i];
            if (em->b[i][c1] != em->b[i][c1]) {
                fprintf(stderr, "b is NAN!\n");
                for (int c2 = 0; c2 < nComps; c2++) {
                    fprintf(stderr, "c2=%d -> b[i+1][c2] = %.2e, scale = %.2e\n", c2, em->b[i + 1][c2], em->scales[i]);
                }
                exit(EXIT_FAILURE);
            }
        }
    }
    // Update P(x)
    // TODO: P(x) is calculated wrongly here p(x) = s(1) x s(2) x ... x s(L) check Durbin's
    em->px = 0;
    for (int c = 0; c < nComps; c++) {
        // starting with c
        tProb = trans[seqClass[0]]->data[nComps][c];
        if(model->modelType == GAUSSIAN) {
            eProb = getGaussianProb(seqEmit[0], model->emit[seqClass[0]][c], c);
        }
        else if(model->modelType == NEGATIVE_BINOMIAL) {
            eProb = NegativeBinomial_getProb(seqEmit[0], model->emit[seqClass[0]][c], c);
        }
        em->px += em->b[0][c] * eProb * tProb;
        if (em->px != em->px) {
            fprintf(stderr, "px is NAN!\n");
            exit(EXIT_FAILURE);
        }
    }
}

double *getForward(EM *em, int pos) {
    double *f = malloc(em->nComps * sizeof(double));
    for (int c = 0; c < em->nComps; c++) {
        f[c] = em->f[pos][c];
    }
    return f;
}

double *getBackward(EM *em, int pos) {
    double *b = malloc(em->nComps * sizeof(double));
    for (int c = 0; c < em->nComps; c++) {
        b[c] = em->b[pos][c];
    }
    return b;
}

double *getPosterior(EM *em, int pos) {
    double *posterior = malloc(em->nComps * sizeof(double));
    for (int c = 0; c < em->nComps; c++) {
        posterior[c] = em->f[pos][c] * em->b[pos][c] / em->px * em->scales[pos];
    }
    return posterior;
}

void NegativeBinomial_resetSufficientStats(HMM *model) {

    for (int r = 0; r < model->nClasses; r++) {
        MatrixDouble_setValue(model->transCounts[r], 0);
    }

    for (int r = 0; r < model->nClasses; r++) {
        for (int c = 0; c < model->nComps; c++) {
            NegativeBinomial *nb = model->emit[r][c];
            for (int m = 0; m < nb->n; m++) {
                VectorDouble_setValue(nb->thetaNum[m], 0.0);
                VectorDouble_setValue(nb->thetaDenom[m], 0.0);
                VectorDouble_setValue(nb->lambdaNum[m], 0.0);
                VectorDouble_setValue(nb->lambdaDenom[m], 0.0);
                nb->weightNum[m] = 0.0;
            }
        }
    }
}

void NegativeBinomial_estimateParameters(HMM *model) {
    // Get model attributes
    int nComps = model->nComps;
    int nClasses = model->nClasses;
    int nEmit = 1;//model->nEmit;
    MatrixDouble **trans = model->trans;
    // A 1D array of count matrices [nClasses]
    MatrixDouble **transCounts = model->transCounts;
    // A 1D array of pseudo-count matrices for denominator
    MatrixDouble **transDenom = model->transDenom;
    // A 1D array of pseudo-count matrices for numarator
    MatrixDouble **transNum = model->transNum;

    // Update transition probs
    double terminateProb = model->terminateProb;
    for (int r = 0; r < nClasses; r++) {
        fprintf(stderr, "transition counts\nr=%d\n", r);
        fprintf(stderr, "%s\n", MatrixDouble_toString(transCounts[r]));
        for (int c1 = 0; c1 < nComps; c1++) {
            double norm = 0;
            for (int c2 = 0; c2 < nComps; c2++) {
                // transition counts for region class of r
                // transition from state c1 to state c2
                norm += transCounts[r]->data[c1][c2];
            }
            if (norm < 10) continue;
            for (int c2 = 0; c2 < nComps; c2++) {
                if (c1 == 2) {
                    fprintf(stdout, "r=%d, 2->%d: %.2e/%.2e\n", r, c2, transCounts[r]->data[c1][c2], norm);
                }
                trans[r]->data[c1][c2] = (transCounts[r]->data[c1][c2] + transNum[r]->data[c1][c2]) /
                                         (norm + transDenom[r]->data[c1][c2]) * (1.0 - terminateProb);
            }
            trans[r]->data[c1][nComps] = terminateProb;
        }
        for (int c2 = 0; c2 < nComps; c2++) {
            trans[r]->data[nComps][c2] = 1.0 / nComps * (1.0 - terminateProb);
        }
        trans[r]->data[nComps][nComps] = terminateProb;
    }


    // Update mu
    for (int r = 0; r < nClasses; r++) {
        VectorDouble *thetaPosFactor = VectorDouble_construct0(nEmit);
        VectorDouble *thetaPosFactorWeight = VectorDouble_construct0(nEmit);
        VectorDouble *lambdaPosFactor = VectorDouble_construct0(nEmit);
        VectorDouble *lambdaPosFactorWeight = VectorDouble_construct0(nEmit);
        // Sum all numerator and denomerator for estimating
        // the theta_0 with positive factors
        for (int c = 0; c < nComps; c++) {
            NegativeBinomial *nb = model->emit[r][c];
            for (int m = 0; m < nb->n; m++) {
                VectorDouble *muFactor = model->muFactors[c][m];
                for (int j = 0; j < nEmit; j++) {
                    if (0 < muFactor->data[j]) {
                        // for theta
                        thetaPosFactor->data[j] += nb->thetaNum[m]->data[j];
                        thetaPosFactorWeight->data[j] += nb->thetaDenom[m]->data[j];
                        // for lambda
                        lambdaPosFactor->data[j] += nb->lambdaNum[m]->data[j];
                        lambdaPosFactorWeight->data[j] += nb->lambdaDenom[m]->data[j];
                    }
                }
            }
        }
        // Update the theta vectors
        for (int c = 0; c < nComps; c++) {
            NegativeBinomial *nb = model->emit[r][c];
            for (int m = 0; m < nb->n; m++) {
                VectorDouble *muFactor = model->muFactors[c][m];
                for (int j = 0; j < nEmit; j++) {
                    //if (muFactor->data[j] < 0) continue;
                    // for theta
                    double theta_0=0.0;
                    if (0 < muFactor->data[j] && 1e2 < thetaPosFactorWeight->data[j]) {
                        theta_0 = thetaPosFactor->data[j] / thetaPosFactorWeight->data[j];
                    } else if (1e2 < nb->thetaDenom[m]->data[j]) {
                        theta_0 = nb->thetaNum[m]->data[j] / nb->thetaDenom[m]->data[j];
                    }
                    // for lambda
                    double lambda_0 =0.0;
                    if (0 < muFactor->data[j] && 1e2 < lambdaPosFactorWeight->data[j]) {
                        lambda_0 = lambdaPosFactor->data[j] / lambdaPosFactorWeight->data[j];
                    } else if (1e2 < nb->lambdaDenom[m]->data[j]) {
                        lambda_0 = nb->lambdaNum[m]->data[j] / nb->lambdaDenom[m]->data[j];
                    }
                    if((theta_0 > 0.0) && (lambda_0>0.0)) {
                        // get update theta and r
                        double theta_m = theta_0;
                        double lambda_m = lambda_0 * muFactor->data[j];
                        double r_m = -1 * lambda_m / log(theta_m);
                        assert(j == 0);
                        // update mu and cov
                        nb->mu[m]->data[j] = NegativeBinomial_getMean(theta_m, r_m);
                        nb->cov[m]->data[j][j] = NegativeBinomial_getVar(theta_m, r_m);
                    }
                }
            }
        }
        VectorDouble_destruct(thetaPosFactor);
        VectorDouble_destruct(thetaPosFactorWeight);
        VectorDouble_destruct(lambdaPosFactor);
        VectorDouble_destruct(lambdaPosFactorWeight);
    }


    // Update mixture weights
    double weightDenom;
    for (int r = 0; r < nClasses; r++) {
        for (int c = 0; c < nComps; c++) {
            NegativeBinomial *nb = model->emit[r][c];
            weightDenom = 0.0;
            for (int m = 0; m < nb->n; m++) {
                weightDenom += nb->weightNum[m];
            }
            if (1.0 < weightDenom) {
                for (int m = 0; m < nb->n; m++) {
                    nb->weights[m] = nb->weightNum[m] / weightDenom;
                }
            }
        }
    }
}


void estimateParameters(HMM *model) {
    // Get model attributes
    int nComps = model->nComps;
    int nClasses = model->nClasses;
    int nEmit = 1;//model->nEmit;
    MatrixDouble **trans = model->trans;
    // A 1D array of count matrices [nClasses]
    MatrixDouble **transCounts = model->transCounts;
    // A 1D array of pseudo-count matrices for denominator
    MatrixDouble **transDenom = model->transDenom;
    // A 1D array of pseudo-count matrices for numarator
    MatrixDouble **transNum = model->transNum;

    // Update transition probs
    double terminateProb = model->terminateProb;
    for (int r = 0; r < nClasses; r++) {
        fprintf(stderr, "transition counts\nr=%d\n", r);
        fprintf(stderr, "%s\n", MatrixDouble_toString(transCounts[r]));
        for (int c1 = 0; c1 < nComps; c1++) {
            double norm = 0;
            for (int c2 = 0; c2 < nComps; c2++) {
                // transition counts for region class of r
                // transition from state c1 to state c2
                norm += transCounts[r]->data[c1][c2];
            }
            if (norm < 10) continue;
            for (int c2 = 0; c2 < nComps; c2++) {
                if (c1 == 2) {
                    fprintf(stdout, "r=%d, 2->%d: %.2e/%.2e\n", r, c2, transCounts[r]->data[c1][c2], norm);
                }
                trans[r]->data[c1][c2] = (transCounts[r]->data[c1][c2] + transNum[r]->data[c1][c2]) /
                                         (norm + transDenom[r]->data[c1][c2]) * (1.0 - terminateProb);
            }
            trans[r]->data[c1][nComps] = terminateProb;
        }
        for (int c2 = 0; c2 < nComps; c2++) {
            trans[r]->data[nComps][c2] = 1.0 / nComps * (1.0 - terminateProb);
        }
        trans[r]->data[nComps][nComps] = terminateProb;
    }


    // Update mu
    for (int r = 0; r < nClasses; r++) {
        VectorDouble *muPosFactor = VectorDouble_construct0(nEmit);
        VectorDouble *muPosFactorWeight = VectorDouble_construct0(nEmit);
        // Sum all numerator and denomerator for estimating
        // the mean vectors with positive factors
        for (int c = 0; c < nComps; c++) {
            Gaussian *gaussian = model->emit[r][c];
            for (int m = 0; m < gaussian->n; m++) {
                VectorDouble *muFactor = model->muFactors[c][m];
                for (int j = 0; j < nEmit; j++) {
                    if (0 < muFactor->data[j]) {
                        muPosFactor->data[j] += gaussian->muNum[m]->data[j];
                        muPosFactorWeight->data[j] += gaussian->muDenom[m]->data[j];
                    }
                }
            }
        }
        // Update the mean vectors
        for (int c = 0; c < nComps; c++) { // skipping erroneous component
            Gaussian *gaussian = model->emit[r][c];
            for (int m = 0; m < gaussian->n; m++) {
                VectorDouble *muFactor = model->muFactors[c][m];
                for (int j = 0; j < nEmit; j++) {
                    //if (muFactor->data[j] < 0) continue;
                    if (0 < muFactor->data[j] && 1e2 < muPosFactorWeight->data[j]) {
                        gaussian->mu[m]->data[j] =
                                muPosFactor->data[j] / muPosFactorWeight->data[j] * muFactor->data[j];
                    } else if (1e2 < gaussian->muDenom[m]->data[j]) {
                        gaussian->mu[m]->data[j] = gaussian->muNum[m]->data[j] / gaussian->muDenom[m]->data[j];
                    }
                }
            }
        }
        VectorDouble_destruct(muPosFactor);
        VectorDouble_destruct(muPosFactorWeight);
    }

    // Update cov
    for (int r = 0; r < nClasses; r++) {
        MatrixDouble *covPosFactor = MatrixDouble_construct0(nEmit, nEmit);
        MatrixDouble *covPosFactorWeight = MatrixDouble_construct0(nEmit, nEmit);
        for (int c = 0; c < nComps; c++) {
            Gaussian *gaussian = model->emit[r][c];
            for (int m = 0; m < gaussian->n; m++) {
                MatrixDouble *covFactor = model->covFactors[c][m];
                for (int j1 = 0; j1 < nEmit; j1++) {
                    for (int j2 = 0; j2 < nEmit; j2++) {
                        double factor = covFactor->data[j1][j2];
                        if (0 < factor) {
                            covPosFactor->data[j1][j2] += gaussian->covNum[m]->data[j1][j2];
                            covPosFactorWeight->data[j1][j2] += gaussian->covDenom[m]->data[j1][j2];
                        }
                    }
                }
            }
        }
        for (int c = 0; c < nComps; c++) {
            Gaussian *gaussian = model->emit[r][c];
            for (int m = 0; m < gaussian->n; m++) {
                MatrixDouble *covFactor = model->covFactors[c][m];
                for (int j1 = 0; j1 < nEmit; j1++) {
                    for (int j2 = 0; j2 < nEmit; j2++) {
                        double factor = covFactor->data[j1][j2];
                        //if (factor < 0) continue;
                        if (0 < factor && 1e2 < covPosFactorWeight->data[j1][j2]) {
                            gaussian->cov[m]->data[j1][j2] =
                                    covPosFactor->data[j1][j2] / covPosFactorWeight->data[j1][j2] * factor;
                        } else if (1e2 < gaussian->covDenom[m]->data[j1][j2]) {
                            gaussian->cov[m]->data[j1][j2] =
                                    gaussian->covNum[m]->data[j1][j2] / gaussian->covDenom[m]->data[j1][j2];
                        }
                    }
                }
            }
        }
        MatrixDouble_destruct(covPosFactor);
        MatrixDouble_destruct(covPosFactorWeight);
    }


    // Update mixture weights
    double weightDenom;
    for (int r = 0; r < nClasses; r++) {
        for (int c = 0; c < nComps; c++) {
            Gaussian *gaussian = model->emit[r][c];
            weightDenom = 0.0;
            for (int m = 0; m < gaussian->n; m++) {
                weightDenom += gaussian->weightNum[m];
            }
            if (1.0 < weightDenom) {
                for (int m = 0; m < gaussian->n; m++) {
                    gaussian->weights[m] = gaussian->weightNum[m] / weightDenom;
                }
            }
        }
    }
}


void resetSufficientStats(HMM *model) {

    for (int r = 0; r < model->nClasses; r++) {
        MatrixDouble_setValue(model->transCounts[r], 0);
    }

    for (int r = 0; r < model->nClasses; r++) {
        for (int c = 0; c < model->nComps; c++) {
            Gaussian *gaussian = model->emit[r][c];
            for (int m = 0; m < gaussian->n; m++) {
                VectorDouble_setValue(gaussian->muNum[m], 0.0);
                VectorDouble_setValue(gaussian->muDenom[m], 0.0);
                MatrixDouble_setValue(gaussian->covNum[m], 0.0);
                MatrixDouble_setValue(gaussian->covDenom[m], 0.0);
                gaussian->weightNum[m] = 0.0;
            }
        }
    }
}

void NegativeBinomial_updateSufficientStats(HMM *model, EM *em) {
    // Get em attributes
    VectorChar **seqEmit = em->seqEmit;
    uint8_t *seqClass = em->seqClass;
    int seqLength = em->seqLength;
    // Get model attributes
    int nComps = model->nComps;
    int nClasses = model->nClasses;
    int nEmit = 1;//model->nEmit;
    MatrixDouble **trans = model->trans;
    // A 1D array of count matrices [nClasses]
    MatrixDouble **transCounts = model->transCounts;
    MatrixDouble **transCountsTemp = MatrixDouble_constructArray1D(nClasses, nComps + 1, nComps + 1);

    double tProb;
    double eProb;
    uint8_t r1;
    uint8_t r2;
    uint8_t r;
    double terminateProb = model->terminateProb;



    // Update transCounts
    for (int i = 0; i < seqLength - 1; i++) {
        r1 = seqClass[i];
        r2 = seqClass[i + 1];
        for (int c1 = 0; c1 < nComps; c1++) {
            for (int c2 = 0; c2 < nComps; c2++) {
                if (r1 == r2) {
                    double highMapqRatio = (double) seqEmit[i + 1]->data[1] / seqEmit[i + 1]->data[0];
                    if (model->maxHighMapqRatio < highMapqRatio) {
                        tProb = c2 == 1 ? 0 : trans[r1]->data[c1][c2] / (1 - trans[r1]->data[c1][1]);
                    } else {
                        tProb = trans[r1]->data[c1][c2];
                    }
                    eProb = NegativeBinomial_getProb(seqEmit[i + 1], model->emit[r2][c2], c2);
                    // division by terminateProb is just for making the total count
                    // equal to the number of windows, which is ~ genome_size / window_size
                    double u = em->f[i][c1] * tProb * eProb * em->b[i + 1][c2] / terminateProb;
                    if (u > 1) {
                        fprintf(stderr, "Error!\n");
                        fprintf(stderr, "i=%d,c1=%d,c2=%d,r1=%d,r2=%d\n", i, c1, c2, r1, r2);
                        fprintf(stderr, "u= %.3e\n", u);
                        fprintf(stderr, "f= %.3e\n", em->f[i][c1]);
                        fprintf(stderr, "b= %.3e\n", em->b[i + 1][c2]);
                        fprintf(stderr, "t= %.3e\n", tProb);
                        fprintf(stderr, "e= %.3e\n", eProb);
                    }
                    transCountsTemp[r1]->data[c1][c2] += u;
                    if (u != u) {
                        fprintf(stderr, "Updating transcounts NAN observed\n");
                        fprintf(stderr, "f= %2.e\n", em->f[i][c1]);
                        fprintf(stderr, "b= %.2e\n", em->b[i + 1][c2]);
                        fprintf(stderr, "px= %.2e\n", em->px);
                        fprintf(stderr, "tprob=%.2e\n", tProb);
                        fprintf(stderr, "eprob=%.2e\n", eProb);
                        exit(EXIT_FAILURE);
                    }
                    if (transCountsTemp[r1]->data[c1][c2] != transCountsTemp[r1]->data[c1][c2]) {
                        fprintf(stderr, "Updating transcounts NAN observed\n");
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }
    }

    pthread_mutex_lock(model->mutexPtr);
    for (int r = 0; r < nClasses; r++) {
        for (int c1 = 0; c1 < nComps; c1++) {
            for (int c2 = 0; c2 < nComps; c2++) {
                transCounts[r]->data[c1][c2] += transCountsTemp[r]->data[c1][c2];
            }
        }
    }
    pthread_mutex_unlock(model->mutexPtr);
    MatrixDouble_destructArray1D(transCountsTemp, nClasses);


    //fprintf(stderr, "Copy gaussian\n");
    // Copy Gaussian params from model to em
    for (int r = 0; r < nClasses; r++) {
        for (int c = 0; c < nComps; c++) {
            NegativeBinomial *nbEm = em->emit[r][c];
            NegativeBinomial *nbModel = model->emit[r][c];
            for (int m = 0; m < nbEm->n; m++) {
                VectorDouble_copyInPlace(nbEm->mu[m], nbModel->mu[m]);
                MatrixDouble_copyInPlace(nbEm->cov[m], nbModel->cov[m]);
                nbEm->weights[m] = nbModel->weights[m];
            }
        }
    }

    //fprintf(stderr, "Update means\n");
    // TODO: Add weightComps[] to make means dependent on each other
    double w = 0.0;
    double w1 = 0.0;
    double w2 = 0.0;
    // Update the numerators and denomerators for estimating mu
    // Update the numerator for estimating mixture weights
    for (int i = 0; i < seqLength - 1; i++) {
        r = seqClass[i];
        for (int c = 0; c < nComps; c++) { // skip erroneous component
            //fprintf(stderr, "c=%d\n",c);
            NegativeBinomial *nb = em->emit[r][c];
            double *mixtureProbs = NegativeBinomial_getMixtureProbs(seqEmit[i], nb, c);
            double totProb = 0;
            for (int m = 0; m < nb->n; m++) {
                totProb += mixtureProbs[m];
            }
            for (int m = 0; m < nb->n; m++) {
                double theta = NegativeBinomial_getTheta(nb->mu[m]->data[0], nb->cov[m]->data[0][0]);
                double r = NegativeBinomial_getR(nb->mu[m]->data[0], nb->cov[m]->data[0][0]);
                double beta = -1 * theta / (1 - theta) - 1 /log(theta);
                //fprintf(stderr, "get factors\n");
                VectorDouble *muFactor = model->muFactors[c][m];
                MatrixDouble *covFactor = model->covFactors[c][m];
                w1 = mixtureProbs[m] / totProb;
                w2 = em->f[i][c] * em->b[i][c] * em->scales[i] / terminateProb;
                w = w1 * w2;
                //fprintf(stderr, "means\n");
                // Update sufficient stats for estimating mean vectors
                for (int j = 0; j < nEmit; j++) {
                    double factor = muFactor->data[j] <= 0 ? 1 : muFactor->data[j];
                    double delta = r*(digammal(r + seqEmit[i]->data[j]) - digammal(r));
                    nb->lambdaNum[m]->data[j] += w * delta;
                    nb->lambdaDenom[m]->data[j] += w * factor;
                    nb->thetaNum[m]->data[j] += w * delta * beta;
                    nb->thetaDenom[m]->data[j] += w * delta * beta + w * (seqEmit[i]->data[j] - delta);
                }
                // Update sufficient stats for estimating mixture weights
                nb->weightNum[m] += w;
            }
            free(mixtureProbs);
        }
    }


    //fprintf(stderr, "Update means 2\n");
    pthread_mutex_lock(model->mutexPtr);
    for (int r = 0; r < nClasses; r++) {
        for (int c = 0; c < nComps; c++) {
            NegativeBinomial *nbEm = em->emit[r][c];
            NegativeBinomial *nbModel = model->emit[r][c];
            for (int m = 0; m < nbEm->n; m++) {
                // Update model gaussian sufficient stats
                //
                // For means
                for (int j = 0; j < nEmit; j++) {
                    nbModel->thetaNum[m]->data[j] += nbEm->thetaNum[m]->data[j];
                    nbModel->thetaDenom[m]->data[j] += nbEm->thetaDenom[m]->data[j];
                    nbModel->lambdaNum[m]->data[j] += nbEm->lambdaNum[m]->data[j];
                    nbModel->lambdaDenom[m]->data[j] += nbEm->lambdaDenom[m]->data[j];
                }
                // For weights
                nbModel->weightNum[m] += nbEm->weightNum[m];
            }
        }
    }
    pthread_mutex_unlock(model->mutexPtr);
}


void updateSufficientStats(HMM *model, EM *em) {
    // Get em attributes
    VectorChar **seqEmit = em->seqEmit;
    uint8_t *seqClass = em->seqClass;
    int seqLength = em->seqLength;
    // Get model attributes
    int nComps = model->nComps;
    int nClasses = model->nClasses;
    int nEmit = 1;//model->nEmit;
    MatrixDouble **trans = model->trans;
    // A 1D array of count matrices [nClasses]
    MatrixDouble **transCounts = model->transCounts;
    MatrixDouble **transCountsTemp = MatrixDouble_constructArray1D(nClasses, nComps + 1, nComps + 1);

    double tProb;
    double eProb;
    uint8_t r1;
    uint8_t r2;
    uint8_t r;
    double terminateProb = model->terminateProb;



    // Update transCounts
    for (int i = 0; i < seqLength - 1; i++) {
        r1 = seqClass[i];
        r2 = seqClass[i + 1];
        for (int c1 = 0; c1 < nComps; c1++) {
            for (int c2 = 0; c2 < nComps; c2++) {
                if (r1 == r2) {
                    double highMapqRatio = (double) seqEmit[i + 1]->data[1] / seqEmit[i + 1]->data[0];
                    if (model->maxHighMapqRatio < highMapqRatio) {
                        tProb = c2 == 1 ? 0 : trans[r1]->data[c1][c2] / (1 - trans[r1]->data[c1][1]);
                    } else {
                        tProb = trans[r1]->data[c1][c2];
                    }
                    eProb = getGaussianProb(seqEmit[i + 1], model->emit[r2][c2], c2);
                    // division by terminateProb is just for making the total count
                    // equal to the number of windows, which is ~ genome_size / window_size
                    double u = em->f[i][c1] * tProb * eProb * em->b[i + 1][c2] / terminateProb;
                    if (u > 1) {
                        fprintf(stderr, "Error!\n");
                        fprintf(stderr, "i=%d,c1=%d,c2=%d,r1=%d,r2=%d\n", i, c1, c2, r1, r2);
                        fprintf(stderr, "u= %.3e\n", u);
                        fprintf(stderr, "f= %.3e\n", em->f[i][c1]);
                        fprintf(stderr, "b= %.3e\n", em->b[i + 1][c2]);
                        fprintf(stderr, "t= %.3e\n", tProb);
                        fprintf(stderr, "e= %.3e\n", eProb);
                    }
                    transCountsTemp[r1]->data[c1][c2] += u;
                    if (u != u) {
                        fprintf(stderr, "Updating transcounts NAN observed\n");
                        fprintf(stderr, "f= %2.e\n", em->f[i][c1]);
                        fprintf(stderr, "b= %.2e\n", em->b[i + 1][c2]);
                        fprintf(stderr, "px= %.2e\n", em->px);
                        fprintf(stderr, "tprob=%.2e\n", tProb);
                        fprintf(stderr, "eprob=%.2e\n", eProb);
                        exit(EXIT_FAILURE);
                    }
                    if (transCountsTemp[r1]->data[c1][c2] != transCountsTemp[r1]->data[c1][c2]) {
                        fprintf(stderr, "Updating transcounts NAN observed\n");
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }
    }

    pthread_mutex_lock(model->mutexPtr);
    for (int r = 0; r < nClasses; r++) {
        for (int c1 = 0; c1 < nComps; c1++) {
            for (int c2 = 0; c2 < nComps; c2++) {
                transCounts[r]->data[c1][c2] += transCountsTemp[r]->data[c1][c2];
            }
        }
    }
    pthread_mutex_unlock(model->mutexPtr);
    MatrixDouble_destructArray1D(transCountsTemp, nClasses);


    //fprintf(stderr, "Copy gaussian\n");
    // Copy Gaussian params from model to em
    for (int r = 0; r < nClasses; r++) {
        for (int c = 0; c < nComps; c++) {
            Gaussian *gaussianEm = em->emit[r][c];
            Gaussian *gaussianModel = model->emit[r][c];
            for (int m = 0; m < gaussianEm->n; m++) {
                VectorDouble_copyInPlace(gaussianEm->mu[m], gaussianModel->mu[m]);
                MatrixDouble_copyInPlace(gaussianEm->cov[m], gaussianModel->cov[m]);
                gaussianEm->weights[m] = gaussianModel->weights[m];
            }
        }
    }

    //fprintf(stderr, "Update means\n");
    // TODO: Add weightComps[] to make means dependent on each other
    double w = 0.0;
    double w1 = 0.0;
    double w2 = 0.0;
    // Update the numerators and denomerators for estimating mu
    // Update the numerator for estimating mixture weights
    for (int i = 0; i < seqLength - 1; i++) {
        r = seqClass[i];
        for (int c = 0; c < nComps; c++) { // skip erroneous component
            //fprintf(stderr, "c=%d\n",c);
            Gaussian *gaussian = em->emit[r][c];
            double *mixtureProbs = getGaussianMixtureProbs(seqEmit[i], gaussian, c);
            double totGaussianProb = getGaussianProb(seqEmit[i], gaussian, c);
            for (int m = 0; m < gaussian->n; m++) {
                //fprintf(stderr, "get factors\n");
                VectorDouble *muFactor = model->muFactors[c][m];
                MatrixDouble *covFactor = model->covFactors[c][m];
                w1 = mixtureProbs[m] / totGaussianProb;
                w2 = em->f[i][c] * em->b[i][c] * em->scales[i] / terminateProb;
                w = w1 * w2;
                //fprintf(stderr, "means\n");
                // Update sufficient stats for estimating mean vectors
                for (int j = 0; j < nEmit; j++) {
                    double factor = muFactor->data[j] <= 0 ? 1 : muFactor->data[j];
                    gaussian->muNum[m]->data[j] += w * seqEmit[i]->data[j] / factor;
                    gaussian->muDenom[m]->data[j] += w;
                }
                // Update sufficient stats for estimating mixture weights
                gaussian->weightNum[m] += w;
                //fprintf(stderr, "covs\n");
                // Update sufficient stats for estimating covariance matrices
                for (int j1 = 0; j1 < nEmit; j1++) {
                    for (int j2 = 0; j2 < nEmit; j2++) {
                        double factor = covFactor->data[j1][j2] <= 0 ? 1 : covFactor->data[j1][j2];
                        double z1 = seqEmit[i]->data[j1] - gaussian->mu[m]->data[j1];
                        double z2 = seqEmit[i]->data[j2] - gaussian->mu[m]->data[j2];
                        gaussian->covNum[m]->data[j1][j2] += w * z1 * z2 / factor;
                        //fprintf(stderr, "cov denom\n");
                        gaussian->covDenom[m]->data[j1][j2] += w;
                    }
                }
            }
            free(mixtureProbs);
        }
    }


    //fprintf(stderr, "Update means 2\n");
    pthread_mutex_lock(model->mutexPtr);
    for (int r = 0; r < nClasses; r++) {
        for (int c = 0; c < nComps; c++) {
            Gaussian *gaussianEm = em->emit[r][c];
            Gaussian *gaussianModel = model->emit[r][c];
            for (int m = 0; m < gaussianEm->n; m++) {
                // Update model gaussian sufficient stats
                //
                // For means
                for (int j = 0; j < nEmit; j++) {
                    gaussianModel->muNum[m]->data[j] += gaussianEm->muNum[m]->data[j];
                    gaussianModel->muDenom[m]->data[j] += gaussianEm->muDenom[m]->data[j];
                }
                // For covs
                for (int j1 = 0; j1 < nEmit; j1++) {
                    for (int j2 = 0; j2 < nEmit; j2++) {
                        gaussianModel->covNum[m]->data[j1][j2] += gaussianEm->covNum[m]->data[j1][j2];
                        gaussianModel->covDenom[m]->data[j1][j2] += gaussianEm->covDenom[m]->data[j1][j2];
                    }
                }
                // For weights
                gaussianModel->weightNum[m] += gaussianEm->weightNum[m];
            }
        }
    }
    pthread_mutex_unlock(model->mutexPtr);
}


Chunk *Chunk_construct1(int chunkLen) {
    Chunk *chunk = malloc(sizeof(Chunk));
    chunk->chunkLen = chunkLen;
    chunk->seqEmit = NULL;
    chunk->seqClass = NULL;
    chunk->seqLen = 0;
    chunk->ctg[0] = '\0';
    chunk->ctgLen = 0;
    chunk->s = -1;
    chunk->e = -1;
    chunk->windowLen = -1;
    chunk->windowItr = -1;
    chunk->windowSumEmit = NULL;
    chunk->windowClass = NULL;
    chunk->fileOffset = 0;
    return chunk;
}

Chunk *Chunk_construct3(int chunkLen, int emissionDim, int windowLen) {
    Chunk *chunk = malloc(sizeof(Chunk));
    chunk->chunkLen = chunkLen;
    chunk->maxSeqSize = chunkLen * 2 / windowLen + 1;
    chunk->seqEmit = VectorChar_constructArray1D(chunk->maxSeqSize, emissionDim);
    chunk->seqClass = malloc(chunk->maxSeqSize * sizeof(uint8_t));
    chunk->seqLen = 0;
    chunk->ctg[0] = '\0';
    chunk->ctgLen = 0;
    chunk->s = -1;
    chunk->e = -1;
    chunk->windowLen = windowLen;
    chunk->windowItr = -1;
    chunk->windowSumEmit = VectorDouble_construct0(emissionDim);
    chunk->windowClass = malloc(windowLen * sizeof(uint8_t));
    chunk->fileOffset = 0;
    return chunk;
}

void Chunk_destruct(Chunk *chunk) {
    if (chunk->seqEmit) {
        VectorChar_destructArray1D(chunk->seqEmit, chunk->maxSeqSize);
    }
    free(chunk->seqClass);
    if (chunk->windowSumEmit) {
        VectorDouble_destruct(chunk->windowSumEmit);
    }
    free(chunk->windowClass);
    free(chunk);
}

Batch *Batch_construct(char *covPath, int chunkLen, int nThreads, int emissionDim, int windowLen) {
    Batch *batch = malloc(sizeof(Batch));
    batch->nThreads = nThreads;
    batch->chunkLen = chunkLen;
    batch->nEmit = emissionDim;
    batch->threadChunks = (Chunk **) malloc(nThreads * sizeof(Chunk *));
    for (int i = 0; i < nThreads; i++) {
        batch->threadChunks[i] = Chunk_construct3(chunkLen, emissionDim, windowLen);
    }
    char covIndexPath[200];
    sprintf(covIndexPath, "%s.index", covPath);
    batch->templateChunks = parseCovIndex(covIndexPath);
    batch->nThreadChunks = 0;
    batch->templateChunkIdx = 0;
    strcpy(batch->covPath, covPath);
    batch->windowLen = windowLen;
    batch->mutex = malloc(sizeof(pthread_mutex_t));
    pthread_mutex_init(batch->mutex, NULL);
    return batch;
}

void Batch_destruct(Batch *batch) {
    for (int i = 0; i < batch->nThreads; i++) {
        Chunk_destruct(batch->threadChunks[i]);
    }
    free(batch->threadChunks);
    stList_destruct(batch->templateChunks);
    free(batch->mutex);
    free(batch);
}

uint8_t Chunk_getWindowClass(Chunk *chunk) {
    uint8_t nClass = maxCharArray(chunk->windowClass, chunk->windowItr + 1);
    int *counts = malloc((nClass + 1) * sizeof(int));
    memset(counts, 0, (nClass + 1) * sizeof(int));
    for (int i = 0; i < chunk->windowItr + 1; i++) {
        counts[chunk->windowClass[i]] += 1;
    }
    uint8_t classMax = 0;
    int maxCount = counts[0];
    for (int i = 1; i < nClass + 1; i++) {
        if (maxCount < counts[i]) {
            classMax = i;
            maxCount = counts[i];
        }
    }
    free(counts);
    return classMax;
}

int Chunk_addWindow(Chunk *chunk) {
    if (chunk->windowItr == -1) return 0;
    VectorDouble *windowSum = chunk->windowSumEmit;
    VectorDouble_divideByValue(windowSum, chunk->windowItr + 1); // take average
    for (int j = 0; j < windowSum->dim; j++) {
        chunk->seqEmit[chunk->seqLen]->data[j] =
                255 < round(windowSum->data[j]) ? 255 : round(windowSum->data[j]); // windowSum is actually average here
    }
    VectorDouble_setValue(chunk->windowSumEmit, 0.0);
    chunk->seqClass[chunk->seqLen] = Chunk_getWindowClass(chunk);
    chunk->seqLen += 1;
    chunk->windowItr = -1;
    return 1;
}


int Chunk_addBlock(Chunk *chunk, Block_t *block) {
    // check if the contig name matches
    assert(strcmp(chunk->ctg, block->ctg) == 0);
    int nucLen = chunk->windowLen * chunk->seqLen + chunk->windowItr + 1;
    // check if the overlap starts exactly one base after where the already added sequence ends
    //fprintf(stderr, "%d + %d * %d + %d + 1 == %d\n", chunk->s , chunk->windowLen, chunk->seqLen, chunk->windowItr, max(block->s, chunk->s));
    assert(chunk->s + nucLen == max(block->s, chunk->s));
    int nucLenToAdd = min(block->e, chunk->e) - max(block->s, chunk->s) + 1;
    if (nucLenToAdd <= 0) return 0;
    VectorDouble *windowSum;
    for (int i = 0; i < nucLenToAdd; i++) {
        chunk->windowItr += 1;
        chunk->windowItr %= chunk->windowLen;
        // emitted vector at each location
        windowSum = chunk->windowSumEmit;
        for (int j = 0; j < block->attrbsLen - 1; j++) {
            // the attrbs in the block include emission entries and class of region
            windowSum->data[j] += atoi(block->attrbs[j]);
        }
        // the first attribute right after emission is the class of region
        chunk->windowClass[chunk->windowItr] = atoi(block->attrbs[block->attrbsLen - 1]);
        if (chunk->windowItr == chunk->windowLen - 1) { // the window is fully iterated
            Chunk_addWindow(chunk);
        }
    }
    return nucLenToAdd;
}


void writeCovIndex(stList *chunks, char *indexPath) {
    FILE *filePtr = fopen(indexPath, "w+");
    if (filePtr == NULL) {
        printf("Couldn't open %s\n", indexPath);
        exit(EXIT_FAILURE);
    }
    assert(0 < stList_length(chunks));
    Chunk *chunk = stList_get(chunks, 0);
    fprintf(filePtr, "%d\n", chunk->chunkLen);
    for (int i = 0; i < stList_length(chunks); i++) {
        chunk = stList_get(chunks, i);
        fprintf(filePtr, "%s\t%d\t%d\t%d\t%ld\n", chunk->ctg,
                chunk->ctgLen,
                chunk->s,
                chunk->e,
                chunk->fileOffset);
    }
    fflush(filePtr);
    fclose(filePtr);
    fprintf(stderr, "Index file (%s) is written \n", indexPath);
}

stList *createCovIndex(char *covPath, int chunkLen) {
    FILE *filePtr = fopen(covPath, "r+");
    if (filePtr == NULL) {
        printf("Couldn't open %s\n", covPath);
        exit(EXIT_FAILURE);
    }
    Block_t *block = Block_construct(COV);
    uint64_t preFileOffset = ftell(filePtr);
    Chunk *chunk = NULL;
    Chunk *preChunk = NULL;
    stList *chunks = stList_construct3(0, (void (*)(void *)) Chunk_destruct);
    while (Block_next(filePtr, block) == 1) {
        block->s -= 1;
        block->e -= 1; // make them 0-based
        if (chunk == NULL || strcmp(chunk->ctg, block->ctg) != 0) { // contig has changed
            chunk = Chunk_construct1(chunkLen);
            chunk->s = 0;
            chunk->e = block->ctgLen < 2 * chunkLen ? block->ctgLen - 1 : chunkLen - 1;
            strcpy(chunk->ctg, block->ctg);
            chunk->ctgLen = block->ctgLen;
            chunk->fileOffset = preFileOffset;
            stList_append(chunks, chunk);
            //fprintf(stderr, "%d\t%s\t%d\t%d\t%d\t%ld\n", stList_length(chunks), chunk->ctg, chunk->ctgLen, chunk->s, chunk->e, chunk->fileOffset);
        } else if (chunk->e <= block->e && chunk->e < block->ctgLen - 1) { // has reached the end of the chunk
            preChunk = chunk;
            chunk = Chunk_construct1(chunkLen);
            chunk->s = preChunk->e + 1;
            chunk->e = block->ctgLen < preChunk->e + 2 * chunkLen ? block->ctgLen - 1 : preChunk->e + chunkLen;
            strcpy(chunk->ctg, block->ctg);
            chunk->ctgLen = block->ctgLen;
            if (chunk->s <= block->e)
                chunk->fileOffset = preFileOffset;
            else
                chunk->fileOffset = ftell(filePtr);
            stList_append(chunks, chunk);
            //fprintf(stderr, "%d\t%s\t%d\t%d\t%d\t%ld\n", stList_length(chunks), chunk->ctg, chunk->ctgLen, chunk->s, chunk->e, chunk->fileOffset);
        }
        preFileOffset = ftell(filePtr);
    }
    Block_destruct(block);
    return chunks; // The code should reach here when there is no more block left in the cov file
}

stList *parseCovIndex(char *covIndexPath) {
    FILE *filePtr = fopen(covIndexPath, "r+");
    if (filePtr == NULL) {
        printf("Couldn't open %s\n", covIndexPath);
        exit(EXIT_FAILURE);
    }
    size_t len = 0;
    char *line = NULL;
    char *token;
    ssize_t read;
    stList *chunks = stList_construct3(0, (void (*)(void *)) Chunk_destruct);
    read = getline(&line, &len, filePtr);
    int chunkLen = atoi(line);
    while ((read = getline(&line, &len, filePtr)) != -1) {
        Chunk *chunk = Chunk_construct1(chunkLen);
        // contig name
        token = strtok(line, "\t");
        strcpy(chunk->ctg, token);
        // contig length
        token = strtok(NULL, "\t");
        chunk->ctgLen = atoi(token);
        // start
        token = strtok(NULL, "\t");
        chunk->s = atoi(token);
        // end
        token = strtok(NULL, "\t");
        chunk->e = atoi(token);
        // file offset
        token = strtok(NULL, "\t");
        chunk->fileOffset = atol(token);
        // add chunk to the list
        stList_append(chunks, chunk);
    }
    return chunks;
}


int Batch_readThreadChunks(Batch *batch) {
    if (batch->templateChunkIdx == stList_length(batch->templateChunks)) {
        return 0;
    }
    batch->nThreadChunks = 0;
    pthread_t *tids = malloc(batch->nThreads * sizeof(pthread_t));
    int nRunningThreads = 0;
    int nRemainingTemplateChunks = stList_length(batch->templateChunks) - batch->templateChunkIdx;
    while (nRunningThreads < min(batch->nThreads, nRemainingTemplateChunks)) {
        pthread_create(&tids[nRunningThreads], NULL, Batch_readNextChunk, (void *) batch);
        //Batch_readNextChunk((void*) batch);
        nRunningThreads += 1;
    }
    for (int t = 0; t < nRunningThreads; t++) {
        assert(pthread_join(tids[t], NULL) == 0);
    }
    free(tids);
    return 1;
}

// It may return three numbers:
// 1 means that one new chunk is added successfully
// 0 means that the coverage file has no more chunk to add
int Batch_readNextChunk(void *batch_) {
    Batch *batch = (Batch *) batch_;
    // Lock the mutex to stop other threads from reading a new chunk
    pthread_mutex_lock(batch->mutex);
    Chunk *templateChunk = stList_get(batch->templateChunks, batch->templateChunkIdx);
    // Construct a block for iteration
    Block_t *block = Block_construct(COV);
    strcpy(block->ctg, templateChunk->ctg);
    block->ctgLen = templateChunk->ctgLen;
    // Open cov file and jump to the first block of the chunk
    FILE *filePtr = fopen(batch->covPath, "r");
    if (filePtr == NULL) {
        printf("Couldn't open %s\n", batch->covPath);
        exit(EXIT_FAILURE);
    }
    fseek(filePtr, templateChunk->fileOffset, SEEK_SET);
    // Get the chunk that has to be filled with the emitted sequence
    // and set start and end coordinates and also the contig name from the given template chunk
    Chunk *chunk = batch->threadChunks[batch->nThreadChunks];
    batch->templateChunkIdx += 1;
    batch->nThreadChunks += 1;
    // Unlock the mutex to let other threads fetch a chunk to read
    pthread_mutex_unlock(batch->mutex);
    chunk->seqLen = 0;
    chunk->windowItr = -1;
    chunk->s = templateChunk->s;
    chunk->e = templateChunk->e;
    strcpy(chunk->ctg, templateChunk->ctg);
    chunk->ctgLen = templateChunk->ctgLen;
    if( chunk->ctgLen < (5 * chunk->windowLen)){ // make smaller windows for short contigs
        chunk->windowLen = max(1, chunk->windowLen / 5);
    }
    // iterate over the blocks in the cov file
    while (Block_next(filePtr, block) == 1) {
        block->s -= 1;
        block->e -= 1; // make them 0-based
        // if the block overlaps with the chunk
        if (chunk->s <= block->e && block->s <= chunk->e) {
            assert(0 < Chunk_addBlock(chunk, block));
        }
        if (chunk->e <= block->e) { // chunk is read completely
            if (chunk->windowItr != -1) { // handle a partially iterated window
                Chunk_addWindow(chunk);
            }
            Block_destruct(block);
            fclose(filePtr);
            //pthread_mutex_unlock(batch->mutex);
            return 1;
        }
    }
    Block_destruct(block);
    fclose(filePtr);
    //pthread_mutex_unlock(batch->mutex);
    return 0; // The code should reach here when there is no more block left in the cov file
}

stList* Chunk_readAllChunksFromBin(char* covPath, int chunkLen, int windowLen, int nEmit){
    char binPath[1000];
    sprintf(binPath, "%s.chunks.l_%d.w_%d.bin", covPath, chunkLen, windowLen);
    if(! file_exists(binPath)){
        fprintf(stderr, "Error: The bin file %s does not exist. Please use create_bin_chunks for creating the bin file", binPath);
        exit(EXIT_FAILURE);
    }
    FILE *fp = fopen(binPath, "rb");
    int32_t chunkLenInBin;
    int32_t windowLenInBin;
    fread(&(chunkLenInBin), sizeof(int32_t), 1, fp); // first 4 bytes
    fread(&(windowLenInBin), sizeof(int32_t), 1, fp); // second 4 bytes
    if(chunkLenInBin != chunkLen){
        fprintf(stderr, "Error: chunkLen = %d in bin file but it is set to %d by the -l parameter, please use -l %d or fix the bin file", chunkLenInBin, chunkLen, chunkLenInBin);
        exit(EXIT_FAILURE);
    }
    if(windowLenInBin != windowLen){
        fprintf(stderr, "Error: windowLen = %d in bin file but it in set to %d by the -w parameter, please use -w %d or fix the bin file", windowLenInBin, windowLen, windowLenInBin);
        exit(EXIT_FAILURE);
    }
    stList* allChunks = stList_construct3(0, Chunk_destruct);
    while(!feof(fp)) {
        Chunk* chunk = Chunk_construct3(chunkLen, nEmit, windowLen);
        // read the length of contig name + null character
        int32_t ctgNameLen;
        fread(&ctgNameLen, sizeof(int32_t), 1, fp);
        // read the contig name
        fread(chunk->ctg, sizeof(char), ctgNameLen, fp);
        // read start and end
        fread(&(chunk->s), sizeof(int32_t), 1, fp);
        fread(&(chunk->e), sizeof(int32_t), 1, fp);
        // read actual size of chunk
        fread(&(chunk->seqLen), sizeof(int32_t), 1, fp);
        int32_t* cov1Array = malloc(sizeof(int32_t) * chunk->seqLen);
        int32_t* cov2Array = malloc(sizeof(int32_t) * chunk->seqLen);
        // read cov1 and then cov2 array
        fread(cov1Array, sizeof(int32_t), chunk->seqLen, fp);
        fread(cov2Array, sizeof(int32_t), chunk->seqLen, fp);
        for(int i=0; i < chunk->seqLen; i++){
            chunk->seqEmit[i]->data[0] = cov1Array[i];
            chunk->seqEmit[i]->data[1] = cov2Array[i];
        }
        fread(chunk->seqClass, sizeof(int8_t), chunk->seqLen, fp);
        stList_append(allChunks, chunk);
    }
    fclose(fp);
    return allChunks;
}