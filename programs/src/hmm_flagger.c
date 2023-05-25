#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include "sonLib.h"
#include "common.h"
#include "hmm.h"
#include "tpool.h"
#include <pthread.h>

typedef enum StatType {
    FORWARD, BACKWARD, POSTERIOR, EMISSION, TRANSITION
} StatType;
const char *const COMP_COLORS[] = {"162,0,37", "250,104,0", "0,138,0", "170,0,255"};
const char *const COMP_NAMES[] = {"Err", "Dup", "Hap", "Col"};

struct stat st = {0};

pthread_mutex_t mutex;

typedef struct Arguments {
    HMM *model;
    EM *em;
    char name[200];
    char path[200];
    pthread_mutex_t *mutexPtr;
} Arguments;

Arguments *Arguments_construct(HMM *model, EM *em, char *name, char *path, pthread_mutex_t *mutexPtr) {
    Arguments *args = malloc(sizeof(Arguments));
    args->model = model;
    args->em = em;
    args->name[0] = '\0';
    args->path[0] = '\0';
    if (name != NULL)
        strcpy(args->name, name);
    if (path != NULL)
        strcpy(args->path, path);
    args->mutexPtr = mutexPtr;
}


void saveHMMStats(StatType statType, char *path, HMM *model) {
    FILE *fp;

    fp = fopen(path, "w+");
    if (fp == NULL) {
        printf("Couldn't open %s\n", path);
        exit(EXIT_FAILURE);
    }

    // Print transition or emission probs
    for (int r = 0; r < model->nClasses; r++) {
        fprintf(fp, "#r=%d\n------------\n", r);
        if (statType == TRANSITION) {
            char *matrixStr = MatrixDouble_toString(model->trans[r]);
            fprintf(fp, "%s\n\n", matrixStr);
        } else if (statType == EMISSION) {
            for (int c = 0; c < model->nComps; c++) {
                fprintf(fp, "##c=%d\n", c);

                if (model->modelType == GAUSSIAN) {
                    Gaussian *gaussian = model->emit[r][c];
                    for (int m = 0; m < gaussian->n; m++) {
                        fprintf(fp, "##m=%d\n", m);
                        char *vecStr = VectorDouble_toString(gaussian->mu[m]);
                        fprintf(fp, "mu = \n%s\n\n", vecStr);
                        char *matrixStr = MatrixDouble_toString(gaussian->cov[m]);
                        fprintf(fp, "cov = \n%s\n\n", matrixStr);
                        fprintf(fp, "w = %.2e\n\n", gaussian->weights[m]);
                    }
                }
                if (model->modelType == NEGATIVE_BINOMIAL) {
                    NegativeBinomial *nb = model->emit[r][c];
                    for (int m = 0; m < nb->n; m++) {
                        fprintf(fp, "##m=%d\n", m);
                        char *vecStr = VectorDouble_toString(nb->mu[m]);
                        fprintf(fp, "mu = \n%s\n\n", vecStr);
                        char *matrixStr = MatrixDouble_toString(nb->cov[m]);
                        fprintf(fp, "cov = \n%s\n\n", matrixStr);
                        fprintf(fp, "w = %.2e\n\n", nb->weights[m]);
                    }
                }
            }
        }
    }
    fclose(fp);
}


void saveEMStats(StatType statType, char *path, EM *em) {
    FILE *fp;

    fp = fopen(path, "w+");
    if (fp == NULL) {
        printf("Couldn't open %s\n", path);
        exit(EXIT_FAILURE);
    }

    // Print header
    fprintf(fp, "#i\t");
    for (int c = 0; c < em->nComps; c++) {
        fprintf(fp, "comp_%d\t", c);
    }
    fprintf(fp, "Scale\n");

    double *p;
    // Print numbers one position per each line
    for (int i = 0; i < em->seqLength; i++) {
        switch (statType) {
            case FORWARD:
                p = getForward(em, i);
                break;
            case BACKWARD:
                p = getBackward(em, i);
                break;
            case POSTERIOR:
                p = getPosterior(em, i);
                break;
        }
        fprintf(fp, "%d\t", i);
        for (int c = 0; c < em->nComps; c++) {
            fprintf(fp, "%.3E\t", p[c]);
        }
        fprintf(fp, "%.3E\t\n", em->scales[i]);
        free(p);
    }
    fclose(fp);
}

void initMuFourComps(VectorDouble ***mu, int coverage, int *nMixtures) {
    double errCov = coverage > 20 ? (double) coverage / 20 : 1;
    double dupCov = (double) coverage / 2.0;
    double hapCov = (double) coverage;
    double colCov = (double) coverage * 2.0;

    // Erroneous
    mu[0][0]->data[0] = errCov;

    // Duplicated
    mu[1][0]->data[0] = dupCov;

    // Haploid Hom
    mu[2][0]->data[0] = hapCov;

    // Collapsed
    for (int i = 0; i < nMixtures[3]; i++) {
        mu[3][i]->data[0] = colCov * ((double) i * 0.5 + 1);
    }
    fprintf(stderr, "Set collapsed initials\n");
}

HMM *makeAndInitModel(int *coverages, int nClasses, int nComps, int nEmit, int *nMixtures, double maxHighMapqRatio,
                      double *regionFreqRatios, char *numMatrixFile, ModelType modelType) {

    VectorDouble ****mu = malloc(nClasses * sizeof(VectorDouble * **));
    for (int r = 0; r < nClasses; r++) {
        mu[r] = malloc(nComps * sizeof(VectorDouble * *));
        for (int c = 0; c < nComps; c++) {
            mu[r][c] = VectorDouble_constructArray1D(nMixtures[c], nEmit);
        }
    }
    for (int r = 0; r < nClasses; r++) {
        initMuFourComps(mu[r], coverages[r], nMixtures);
    }
    VectorDouble ***muFactors = VectorDouble_constructArray2D(nComps, 10, nEmit);

    //Erroneous
    muFactors[0][0]->data[0] = -1;
    // Duplicated
    muFactors[1][0]->data[0] = 0.5;

    // Haploid
    /// Hom
    muFactors[2][0]->data[0] = 1.0;

    // Collapsed
    for (int i = 0; i < nMixtures[3]; i++) {
        fprintf(stderr, "%d\n", i);
        VectorDouble_setValue(muFactors[3][i], i + 2);
    }
    fprintf(stderr, "Set mu factors\n");

    MatrixDouble ***covFactors = MatrixDouble_constructArray2D(nComps, 10, nEmit, nEmit);

    MatrixDouble_setValue(covFactors[0][0], -1);
    MatrixDouble_setValue(covFactors[1][0], 0.5);

    MatrixDouble_setValue(covFactors[2][0], 1);

    // Collapsed
    for (int i = 0; i < nMixtures[3]; i++) {
        MatrixDouble_setValue(covFactors[3][i], i + 2);
    }
    fprintf(stderr, "Set cov factors\n");

    MatrixDouble *pseudoCountMatrix = MatrixDouble_parseFromFile(numMatrixFile, nComps, nComps);

    fprintf(stderr, "read pseudoCountMatrix\n");
    fprintf(stderr, "pseudoCountMatrix = \n%s\n\n", MatrixDouble_toString(pseudoCountMatrix));
    // Numerator pseudo-counts for updating transition probs
    MatrixDouble **transNum = MatrixDouble_constructArray1D(nClasses, nComps, nComps);

    // Denominator pseudo-counts for updating transition probs
    MatrixDouble **transDenom = MatrixDouble_constructArray1D(nClasses, nComps, nComps);

    double rowSum = 0;
    for (int r = 0; r < nClasses; r++) {
        for (int i = 0; i < nComps; i++) {
            rowSum = 0;
            for (int j = 0; j < nComps; j++) {
                rowSum += pseudoCountMatrix->data[i][j]; //* regionFreqRatios[r];
                transNum[r]->data[i][j] = pseudoCountMatrix->data[i][j]; // * regionFreqRatios[r];
            }
            for (int j = 0; j < nComps; j++) {
                transDenom[r]->data[i][j] = rowSum;
            }
        }
    }

    // Numerator pseudo-counts for updating transition probs
    // MatrixDouble** transDenom = MatrixDouble_constructArray1D(nClasses, nComps, nComps);

    HMM *model = HMM_construct(nClasses, nComps, nEmit, nMixtures, mu, muFactors, covFactors, maxHighMapqRatio,
                               transNum, transDenom, modelType);

    //model->emit[1][1]->cov->data[1][1] = model->emit[1][2]->mu->data[0] / 8.0;
    //model->emit[1][2]->cov->data[1][1] = 4.0 * model->emit[1][2]->mu->data[0];
    //model->emit[1][3]->cov->data[1][1] = 8.0 * model->emit[1][2]->mu->data[0];
    for (int r = 0; r < nClasses; r++) {
        for (int c = 0; c < nComps; c++) {
            VectorDouble_destructArray1D(mu[r][c], nMixtures[c]);
        }
        free(mu[r]);
    }
    free(mu);
    return model;
}

EM **Batch_buildEmArray(Batch *batch, HMM *model) {
    EM **emArray = (EM **) malloc(batch->nThreadChunks * sizeof(EM * ));
    Chunk *chunk;
    for (int t = 0; t < batch->nThreadChunks; t++) {
        chunk = batch->threadChunks[t];
        emArray[t] = EM_construct(chunk->seqEmit, chunk->seqClass, chunk->seqLen, model);
    }
    return emArray;
}

void *trainModelSaveStats(void *args_) {
    Arguments *args = (Arguments *) args_;
    HMM *model = args->model;
    EM *em = args->em;
    char *name = args->name;
    char *dir = args->path;
    pthread_mutex_t *mutexPtr = args->mutexPtr;
    char path[200];
    //for(int itr=1; itr <= nItr; itr++){
    //fprintf(stderr,"\t\tIteration %d\n", itr);
    //fprintf(stderr,"\t\t\tRun forward\n");
    runForward(model, em);
    //fprintf(stderr,"\t\t\tRun borward\n");
    runBackward(model, em);

    sprintf(path, "%s/forward_%s.txt", dir, name);
    saveEMStats(FORWARD, path, em);
    sprintf(path, "%s/backward_%s.txt", dir, name);
    saveEMStats(BACKWARD, path, em);
    //sprintf(path, "%s/posterior_%d.txt", dir, itr);
    //saveEMStats(POSTERIOR, path, em);

    //fprintf(stderr, "\t\t\tUpdate sufficient Stats ...\n");
    pthread_mutex_lock(mutexPtr);
    updateSufficientStats(model, em);
    pthread_mutex_unlock(mutexPtr);

    //estimateParameters(model);
    //resetSufficientStats(model);

    //sprintf(path, "%s/transition_%d.txt", dir, itr);
    //saveHMMStats(TRANSITION, path, model);
    //sprintf(path, "%s/emission_%d.txt", dir, itr);
    //saveHMMStats(EMISSION, path, model);
    //}
}

void *infer(void *args_) {
    Arguments *args = (Arguments *) args_;
    HMM *model = args->model;
    EM *em = args->em;
    char path[200];
    runForward(model, em);
    runBackward(model, em);

    //sprintf(path, "%s/forward_%d.txt", dir, itr);
    //saveEMStats(FORWARD, path, em);
    //sprintf(path, "%s/backward_%d.txt", dir, itr);
    //saveEMStats(BACKWARD, path, em);
    //sprintf(path, "%s/posterior_%d.txt", dir, itr);
    //saveEMStats(POSTERIOR, path, em);
    //}
}


void Batch_trainModelSaveStats(Batch *batch, int batchIdx, HMM *model, char *dir) {
    EM **emArray = Batch_buildEmArray(batch, model);
    char dir1[100];
    char dir2[100];
    char name[200];
    EM *em;
    pthread_t *tids = malloc(batch->nThreadChunks * sizeof(pthread_t));
    Arguments **args = malloc(batch->nThreadChunks * sizeof(Arguments *));
    for (int t = 0; t < batch->nThreadChunks; t++) {
        em = emArray[t];
        /*sprintf(dir1, "%s/batch_%d", dir, batchIdx);
        sprintf(dir2, "%s/batch_%d/chunk_%d", dir, batchIdx, t);
        if (stat(dir1, &st) == -1) {
            mkdir(dir1, 0700);
        }
        if (stat(dir2, &st) == -1) {
                        mkdir(dir2, 0700);
                }*/
        Chunk *chunk = batch->threadChunks[t];
        fprintf(stderr, "Batch %d, Chunk %d, %s: [%d-%d] \n", batchIdx, t, chunk->ctg, chunk->s, chunk->e);
        sprintf(name, "ba%d_ch%d", batchIdx, t);
        args[t] = Arguments_construct(model, em, name, dir, &mutex);
        fprintf(stderr, "Args is created\n");
        pthread_create(&tids[t], NULL, trainModelSaveStats, (void *) args[t]);
        fprintf(stderr, "Thread %d is running\n", t);
    }
    for (int t = 0; t < batch->nThreadChunks; t++) {
        assert(pthread_join(tids[t], NULL) == 0);
        fprintf(stderr, "Thread %d is finished\n", t);
        free(args[t]);
        EM_destruct(emArray[t]);
    }
    free(args);
    free(emArray);
    free(tids);
}

uint8_t getCompIdx(EM *em, int loc) {
    double *p = getPosterior(em, loc);
    uint8_t idx = 0;
    for (int c = 1; c < em->nComps; c++) {
        idx = p[idx] < p[c] ? c : idx;
    }
    free(p);
    return idx;
}

void Batch_inferSaveOutput(Batch *batch, int batchIdx, HMM *model, FILE *outputFile, double minColScore, int minColLen,
                           double maxDupScore, int minDupLen) {
    EM **emArray = Batch_buildEmArray(batch, model);
    char dir1[100];
    char dir2[100];
    EM *em;
    pthread_t *tids = malloc(batch->nThreadChunks * sizeof(pthread_t));
    for (int t = 0; t < batch->nThreadChunks; t++) {
        em = emArray[t];
        Chunk *chunk = batch->threadChunks[t];
        fprintf(stderr, "Batch %d, Chunk %d, %s: [%d-%d] \n", batchIdx, t, chunk->ctg, chunk->s, chunk->e);
        Arguments *args = Arguments_construct(model, em, "infer", NULL, NULL);
        pthread_create(&tids[t], NULL, infer, (void *) args);
        fprintf(stderr, "Thread %d is running\n", t);
    }
    int s = -1;
    int e = -1;
    int compIdx = -1;
    int preCompIdx = -1;
    int windowLen; // windowLen may be different for small contigs; look at chunk->windowLen instead of batch->windowLen
    double cov;
    double hapMu;
    double score;
    double sumScore = 0.0;
    double avgScore = 0.0;
    for (int t = 0; t < batch->nThreadChunks; t++) {
        assert(pthread_join(tids[t], NULL) == 0);
        fprintf(stderr, "Thread %d is finished\n", t);
        fprintf(stderr, "Writing the output...\n");
        Chunk *chunk = batch->threadChunks[t];
        windowLen = chunk->windowLen;
        s = chunk->s;
        e = chunk->s;
        preCompIdx = -1;
        for (int i = 0; i < chunk->seqLen; i++) {
            compIdx = getCompIdx(emArray[t], i);
            cov = (double) chunk->seqEmit[i]->data[0];
            if (model->modelType == GAUSSIAN) {
                Gaussian *gaussian = model->emit[chunk->seqClass[i]][2];
                hapMu = (double) gaussian->mu[0]->data[0]; // TODO: assuming index 2 is always haploid

            } else if (model->modelType == NEGATIVE_BINOMIAL) {
                NegativeBinomial *nb = model->emit[chunk->seqClass[i]][2];
                hapMu = (double) nb->mu[0]->data[0]; // TODO: assuming index 2 is always haploid
            }
            score = cov / hapMu;
            sumScore += score * chunk->windowLen;
            // if component changed write the block
            if (preCompIdx != -1 && preCompIdx != compIdx) {
                avgScore = sumScore / (e - s);
                if (preCompIdx == 1 && avgScore > maxDupScore && e - s < minDupLen) {
                    preCompIdx = 2;
                } else if (preCompIdx == 3 && avgScore < minColScore && e - s < minColLen) {
                    preCompIdx = 2;
                }
                fprintf(outputFile, "%s\t%d\t%d\t%s\t0\t+\t%d\t%d\t%s\t%.2f\n",
                        chunk->ctg, s, e, COMP_NAMES[preCompIdx],
                        s, e, COMP_COLORS[preCompIdx], avgScore);
                s = e;
                sumScore = 0;
            }
            e = i == chunk->seqLen - 1 ? chunk->e + 1 : e + windowLen; // maybe the last window is not complete
            preCompIdx = compIdx;
        }
        avgScore = sumScore / (e - s);
        if (preCompIdx == 1 && avgScore > maxDupScore && e - s < minDupLen) {
            preCompIdx = 2;
        } else if (preCompIdx == 3 && avgScore < minColScore && e - s < minColLen) {
            preCompIdx = 2;
        }
        fprintf(outputFile, "%s\t%d\t%d\t%s\t0\t+\t%d\t%d\t%s\t%0.2f\n",
                chunk->ctg, s, e, COMP_NAMES[preCompIdx],
                s, e, COMP_COLORS[preCompIdx], avgScore); // there will be an unwritten component finally
        fprintf(stderr, "Done writing!\n");
        EM_destruct(emArray[t]);
    }
    free(emArray);
    free(tids);
}

void *readChunkAndUpdateStats(void *arg_) {
    work_arg_t *arg = arg_;
    Batch *batch = arg->batch;
    HMM *model = arg->model;
    // Prepare the batch for reading the chunk
    batch->templateChunkIdx = arg->templateChunkIdx;
    batch->nThreadChunks = 0;
    // Read the chunk
    Batch_readNextChunk(batch);
    Chunk *chunk = batch->threadChunks[0];
    fprintf(stderr, "Chunk %d: %d, %d, %d, ..., %d, %d, %d\n", batch->templateChunkIdx,
            chunk->seqEmit[0]->data[0],
            chunk->seqEmit[1]->data[0],
            chunk->seqEmit[2]->data[0],
            chunk->seqEmit[chunk->seqLen - 3]->data[0],
            chunk->seqEmit[chunk->seqLen - 2]->data[0],
            chunk->seqEmit[chunk->seqLen - 1]->data[0]);
    // Make an EM object
    EM *em = EM_construct(chunk->seqEmit, chunk->seqClass, chunk->seqLen, model);
    // Run forward and backward
    fprintf(stderr, "Run forward\n");
    runForward(model, em);
    fprintf(stderr, "Run backward\n");
    runBackward(model, em);
    fprintf(stderr, "Update sufficient stats\n");
    /*char path[200];
    sprintf(path, "%s/forward.%s.txt", arg->dir, arg->name);
        saveEMStats(FORWARD, path, em);
        sprintf(path, "%s/backward.%s.txt", arg->dir, arg->name);
        saveEMStats(BACKWARD, path, em);*/
    // Update statistics for estimating parameters
    //pthread_mutex_lock(model->mutexPtr);
    if (model->modelType == NEGATIVE_BINOMIAL) {
        NegativeBinomial_updateSufficientStats(model, em);
    } else if (model->modelType == GAUSSIAN) {
        updateSufficientStats(model, em);
    }
    //pthread_mutex_unlock(model->mutexPtr);
    fprintf(stderr, "Update sufficient stats finished\n");
    // Free the EM object
    EM_destruct(em);
}

void *
runOneRound(HMM *model, char *covPath, int chunkLen, int nThreads, int emissionDim, int windowLen, char *dir, int itr) {
    // Prepare all batches for reading chunks
    // Each batch here contain only one chunk
    Batch *batch = Batch_construct(covPath, 0, 0, 0, 0);
    Batch **batches = malloc(nThreads * sizeof(Batch * ));
    for (int t = 0; t < nThreads; t++) {
        batches[t] = Batch_construct(covPath, chunkLen, 1, emissionDim, windowLen);
    }

    // Create a thread pool
    // Each thread recieves only one batch
    tpool_t *tm = tpool_create(nThreads, batches);
    for (int i = 0; i < stList_length(batch->templateChunks); i++) {
        work_arg_t *work_arg = malloc(sizeof(work_arg_t));
        work_arg->templateChunkIdx = i;
        work_arg->model = model;
        strcpy(work_arg->dir, dir);
        sprintf(work_arg->name, "itr_%d.chunk_%d", itr, i);
        tpool_add_work(tm, readChunkAndUpdateStats, work_arg);
    }

    tpool_wait(tm);
    tpool_destroy(tm);



    // Free all batches
    for (int t = 0; t < nThreads; t++) {
        Batch_destruct(batches[t]);
    }
    free(batches);
    Batch_destruct(batch);
}

static struct option long_options[] =
        {
                {"inputCov",         required_argument, NULL, 'i'},
                {"threads",          required_argument, NULL, 't'},
                {"chunkLen",         required_argument, NULL, 'l'},
                {"iterations",       required_argument, NULL, 'n'},
                {"windowLen",        required_argument, NULL, 'w'},
                {"coverage",         required_argument, NULL, 'c'},
                {"regions",          required_argument, NULL, 'r'},
                {"trackName",        required_argument, NULL, 'm'},
                {"outputDir",        required_argument, NULL, 'o'},
                {"regionFactors",    required_argument, NULL, 'f'},
                {"minColScore",      required_argument, NULL, 's'},
                {"minColLen",        required_argument, NULL, 'e'},
                {"maxDupScore",      required_argument, NULL, 'd'},
                {"minDupLen",        required_argument, NULL, 'a'},
                {"maxHighMapqRatio", required_argument, NULL, 'p'},
                {"regionRatios",     required_argument, NULL, 'g'},
                {"transPseudoPath",  required_argument, NULL, 'u'},
                {"model",            required_argument, NULL, 'M'},
                {NULL,               0,                 NULL, 0}
        };


int main(int argc, char *argv[]) {
    int c;
    char covPath[200];
    char covIndexPath[200];
    char outputDir[200];
    char trackName[200];
    int nThreads;
    int chunkLen;
    int windowLen = 1000;
    int nIteration = 5;
    int meanCoverage = -1;
    int nClasses = -1;
    double minColScore = 1.8;
    int minColLen = 40e3;
    double maxDupScore = 0.4;
    int minDupLen = 40e3;
    char regionFactors[200];
    double maxHighMapqRatio = 0.25;
    char regionFreqRatios[200];
    char transPseudoPath[200];
    char *program;
    ModelType modelType;
    (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
    while (~(c = getopt_long(argc, argv, "i:t:l:n:w:c:r:f:g:m:o:s:e:d:a:p:u:M:h", long_options, NULL))) {
        switch (c) {
            case 'i':
                strcpy(covPath, optarg);
                break;
            case 't':
                nThreads = atoi(optarg);
                break;
            case 'o':
                strcpy(outputDir, optarg);
                break;
            case 'l':
                chunkLen = atoi(optarg);
                break;
            case 'w':
                windowLen = atoi(optarg);
                break;
            case 'n':
                nIteration = atoi(optarg);
                break;
            case 'm':
                strcpy(trackName, optarg);
                break;
            case 'c':
                meanCoverage = atoi(optarg);
                break;
            case 'r':
                nClasses = atoi(optarg);
                break;
            case 's':
                minColScore = atof(optarg);
                break;
            case 'e':
                minColLen = atoi(optarg);
                break;
            case 'd':
                maxDupScore = atof(optarg);
                break;
            case 'a':
                minDupLen = atoi(optarg);
                break;
            case 'f':
                strcpy(regionFactors, optarg);
                break;
            case 'p':
                maxHighMapqRatio = atof(optarg);
                break;
            case 'g':
                strcpy(regionFreqRatios, optarg);
                break;
            case 'u':
                strcpy(transPseudoPath, optarg);
                break;
            case 'M':
                if (strcmp(optarg, "gaussian") == 0) {
                    modelType = GAUSSIAN;
                } else if (strcmp(optarg, "nb") == 0) {
                    modelType = NEGATIVE_BINOMIAL;
                } else {
                    fprintf(stderr, "Error: --model should be either 'gaussian' or 'nb'!");
                    exit(EXIT_FAILURE);
                }
                break;
            default:
                if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
            help:
                fprintf(stderr, "\nUsage: %s\n", program);
                fprintf(stderr, "Options:\n");
                fprintf(stderr, "         --inputCov, -i         path to .cov file\n");
                fprintf(stderr,
                        "         --model, -M         use gaussian or negative binomial emission model (values can be either 'gaussian' or 'nb' [default:'nb'])\n");
                fprintf(stderr, "         --threads, -t         number of threads\n");
                fprintf(stderr, "         --chunkLen, -l         chunk length\n");
                fprintf(stderr, "         --iterations, -n         number of iterations\n");
                fprintf(stderr, "         --windowLen, -w         window length (default = 100)\n");
                fprintf(stderr, "         --coverage, -c         mean coverage\n");
                fprintf(stderr, "         --regions, -r         number of region classes [only 2 or 3]\n");
                fprintf(stderr, "         --trackName, -m         track name\n");
                fprintf(stderr, "         --outputDir, -o         output dir\n");
                fprintf(stderr, "         --minColScore, -s       minimum score of short collapsed blocks\n");
                fprintf(stderr, "         --minColLen, -e       minimum length of low score collapsed blocks\n");
                fprintf(stderr, "         --maxDupScore, -d       maximum score of short duplicated blocks\n");
                fprintf(stderr, "         --minDupLen, -a       minimum length of high score duplicated blocks\n");
                fprintf(stderr,
                        "         --regionFactors, -f       factors to scale the initial mean values for different regions\n");
                fprintf(stderr,
                        "         --maxHighMapqRatio, -p       maximum ratio of high mapq coverage for duplicated component\n");
                fprintf(stderr, "         --transPseudoPath, -u       path to transition pseudo counts\n");
                fprintf(stderr,
                        "         --regionFreqRatios, -g	region frequency ratios for adjusting transition pseudo counts\n");

                return 1;
        }
    }
    fprintf(stderr, "CONSTRUCTED MODEL\n");
    int nComps = 4;
    int nEmit = 2;
    int *coverages = malloc(nClasses * sizeof(int));
    Splitter *splitter = Splitter_construct(regionFactors, ',');
    char *token;
    double *numericalFactors = malloc(nClasses * sizeof(double));
    int i = 0;
    while ((token = Splitter_getToken(splitter)) != NULL) {
        numericalFactors[i] = atof(token);
        coverages[i] = meanCoverage * numericalFactors[i];
        i++;
    }
    Splitter_destruct(splitter);
    int *nMixtures = malloc(nComps * sizeof(int));
    nMixtures[0] = 1;
    nMixtures[1] = 1;
    nMixtures[2] = 1;
    nMixtures[3] = 10;
    fprintf(stderr, "CONSTRUCTED MODEL\n");

    // Read region freq ratios to adjust transition pseudo counts
    splitter = Splitter_construct(regionFreqRatios, ',');
    double *regionFreqRatiosDouble = malloc(nClasses * sizeof(double));
    i = 0;
    while ((token = Splitter_getToken(splitter)) != NULL) {
        regionFreqRatiosDouble[i] = atof(token);
        i++;
    }
    Splitter_destruct(splitter);

    HMM *model = makeAndInitModel(coverages, nClasses, nComps, nEmit, nMixtures, maxHighMapqRatio,
                                  regionFreqRatiosDouble, transPseudoPath, modelType);

    fprintf(stderr, "CONSTRUCTED MODEL\n");
    for (int r = 0; r < nClasses; r++) {
        char *transStr = MatrixDouble_toString(model->trans[r]);
        fprintf(stderr, "r=%d, trans=\n%s\n", r, transStr);
        for (int c = 0; c < nComps; c++) {
            if (model->modelType == GAUSSIAN) {
                Gaussian *gaussian = model->emit[r][c];
                for (int m = 0; m < gaussian->n; m++) {
                    char *muStr = VectorDouble_toString(gaussian->mu[m]);
                    char *covStr = MatrixDouble_toString(gaussian->cov[m]);
                    fprintf(stderr, "r=%d, c=%d, m=%d\n%s\n%s\n\n", r, c, m, muStr, covStr);
                }
            } else if (model->modelType == NEGATIVE_BINOMIAL) {
                NegativeBinomial *nb = model->emit[r][c];
                for (int m = 0; m < nb->n; m++) {
                    char *muStr = VectorDouble_toString(nb->mu[m]);
                    char *covStr = MatrixDouble_toString(nb->cov[m]);
                    fprintf(stderr, "r=%d, c=%d, m=%d\n%s\n%s\n\n", r, c, m, muStr, covStr);
                }
            }
        }
    }

    char outputPath[200];
    for (int itr = 0; itr < nIteration; itr++) {
        /*batch->templateChunkIdx = 0;
     batch->nThreadChunks = 0;
        fprintf(stderr, "Batch is built\n");
        batchIdx = 0;
     while(Batch_readThreadChunks(batch)){
         fprintf(stderr, "Run EM for batch %d\n", batchIdx);
                 Batch_trainModelSaveStats(batch, batchIdx, model, outputDir);
         batchIdx += 1;
        }*/
        fprintf(stderr, "ROUND STARTED\n");
        runOneRound(model, covPath, chunkLen, nThreads, nEmit, windowLen, outputDir, itr);
        fprintf(stderr, "ROUND FINISHED\n");

        if (model->modelType == GAUSSIAN) {
            fprintf(stderr, "Estimating parameters for Gaussian model\n");
            estimateParameters(model);
            fprintf(stderr, "Reseting sufficient stats\n");
            resetSufficientStats(model);
        } else if (model->modelType == NEGATIVE_BINOMIAL) {
            fprintf(stderr, "Estimating parameters for Negative Binomial model\n");
            NegativeBinomial_estimateParameters(model);
            fprintf(stderr, "Reseting sufficient stats\n");
            NegativeBinomial_resetSufficientStats(model);
        }

        fprintf(stderr, "Saving HMM stats\n");
        sprintf(outputPath, "%s/transition_%d.txt", outputDir, itr);
        fprintf(stderr, "%s\n", outputPath);
        saveHMMStats(TRANSITION, outputPath, model);
        sprintf(outputPath, "%s/emission_%d.txt", outputDir, itr);
        saveHMMStats(EMISSION, outputPath, model);
    }

    //Run inference
    sprintf(outputPath, "%s/%s.flagger.bed", outputDir, trackName);
    FILE *fp = fopen(outputPath, "w+");
    fprintf(fp, "track name=%s visibility=1 itemRgb=\"On\"\n", trackName);

    Batch *batch = Batch_construct(covPath, chunkLen, nThreads, nEmit, windowLen);
    int batchIdx = 0;
    batch->templateChunkIdx = 0;
    batch->nThreadChunks = 0;
    while (Batch_readThreadChunks(batch)) {
        fprintf(stderr, "[Inference] Running EM for batch %d\n", batchIdx);
        Batch_inferSaveOutput(batch, batchIdx, model, fp, minColScore, minColLen, maxDupScore, minDupLen);
        batchIdx += 1;
    }
    fflush(fp);
    fclose(fp);
    Batch_destruct(batch);
    HMM_destruct(model);
}
