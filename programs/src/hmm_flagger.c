#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include "hmm_utils.h"
#include "hmm.h"
#include "data_types.h"
#include "common.h"
#include "chunk.h"
#include "summary_table.h"

ChunksCreator *getChunksCreator(char *inputPath,
                                int chunkCanonicalLen,
                                int windowLen,
                                int threads,
                                stList *contigList) {
    ChunksCreator *chunksCreator = NULL;

    char *inputExtension = extractFileExtension(inputPath);
    if (strcmp(inputExtension, "cov") == 0 ||
        strcmp(inputExtension, "cov.gz") == 0) {
        char *faiPath = NULL;
        fprintf(stderr, "[%s] The given input file is not binary so chunks will be constructed from cov file.\n",
                get_timestamp());
        chunksCreator = ChunksCreator_constructFromCov(inputPath, faiPath, chunkCanonicalLen, threads, windowLen);
        if (ChunksCreator_parseChunks(chunksCreator) != 0) {
            fprintf(stderr, "[%s] Error: creating chunks from cov file failed.\n", get_timestamp());
            exit(EXIT_FAILURE);
        }
        fprintf(stderr, "[%s] Chunks are constructed from cov file.\n", get_timestamp());
    } else if (strcmp(inputExtension, "bin") == 0) {
        chunksCreator = ChunksCreator_constructEmpty();
        ChunksCreator_parseChunksFromBinaryFile(chunksCreator, inputPath);
        fprintf(stderr, "[%s] Binary chunks are parsed. Based on the header chunkCanonicalLen=%d and windowLen=%d .\n",
                get_timestamp(),
                chunksCreator->chunkCanonicalLen,
                chunksCreator->windowLen);
        if (chunksCreator->chunkCanonicalLen != chunkCanonicalLen || chunksCreator->windowLen != windowLen) {
            fprintf(stderr,
                    "[%s] Warning: chunk/window lengths (based on bin file) do not match the values given through program parameters. The values read from the binary file will be considered.\n",
                    get_timestamp());
        }
    }

    // if there is a list keep only the given contigs
    if (contigList != NULL) {
        fprintf(stderr,
                "[%s] Including only the chunks whose contigs match the given contig list.\n",
                get_timestamp());
        ChunksCreator_subsetChunksToContigs(chunksCreator, contigList);
    }
    free(inputExtension);

    return chunksCreator;
}

int getBestNumberOfCollapsedComps(ChunksCreator *chunksCreator) {
    CoverageHeader *header = chunksCreator->header;
    int maxCoverage = ChunksCreator_getMaximumCoverageValue(chunksCreator);
    int minRegionModeCoverage = Int_getMinValue1DArray(header->regionCoverages, header->numberOfRegions);
    int numberOfCollapsedComps = maxCoverage / minRegionModeCoverage + 1;
    return numberOfCollapsedComps;
}

double getRandomNumber(double start, double end) {
    srand(time(NULL));
    return (double) rand() / (double) (RAND_MAX / (end - start)) + start;
}


void writeParameterStats(HMM *model, char *outputDir, char *suffix) {
    char path[2000];
    /*
    if (writePosterior) {
        fprintf(stderr, "writing posterior tsv...\n");
        sprintf(path, "%s/posterior_prediction_%s.tsv", outputDir, suffix);
        FILE *posteriorTsvFile = fopen(path, "w+");
        EM_printPosteriorInTsvFormat(em, posteriorTsvFile);
        fclose(posteriorTsvFile);
    }*/
    sprintf(path, "%s/transition_%s.tsv", outputDir, suffix);
    fprintf(stderr, "[%s] Writing transition tsv...\n", get_timestamp());
    FILE *transitionTsvFile = fopen(path, "w+");
    HMM_printTransitionMatrixInTsvFormat(model, transitionTsvFile);
    fclose(transitionTsvFile);

    sprintf(path, "%s/emission_%s.tsv", outputDir, suffix);
    fprintf(stderr, "[%s] Writing emission tsv ...\n", get_timestamp());
    FILE *emissionTsvFile = fopen(path, "w+");
    HMM_printEmissionParametersInTsvFormat(model, emissionTsvFile);
    fclose(emissionTsvFile);
}

void writeBenchmarkingStats(ChunksCreator *chunksCreator,
                            char *outputDir,
                            char *suffix,
                            stList *labelNamesWithUnknown,
                            char *binArrayFilePath,
                            double overlapRatioThreshold,
                            int threads) {
    char outputPath[2000];
    sprintf(outputPath, "%s/prediction_summary_%s.tsv", outputDir, suffix);

    // define object and functions for iterating blocks
    BlockIteratorType blockIteratorType = ITERATOR_BY_CHUNK;
    void *iterator = (void *) ChunkIterator_construct(chunksCreator);
    CoverageHeader *header = chunksCreator->header;

    // create and write all summary tables with final stats
    SummaryTableList_createAndWriteAllTables(iterator,
                                             blockIteratorType,
                                             header,
                                             outputPath,
                                             binArrayFilePath,
                                             labelNamesWithUnknown,
                                             overlapRatioThreshold,
                                             threads);

    // free iterator
    ChunkIterator_destruct((ChunkIterator *) iterator);
}


HMM *createModel(ModelType modelType,
                 int numberOfCollapsedComps,
                 CoverageHeader *header,
                 MatrixDouble *alphaMatrix,
                 double maxHighMapqRatio,
                 double initialDeviation) {
    // calculate initial mean coverages for all region classes
    int numberOfStates = 4; // Err, Dup, Hap, and Col
    int numberOfRegions = header->numberOfRegions;

    // set number of mixture components for each state
    int *numberOfCompsPerState = Int_construct1DArray(numberOfStates);
    numberOfCompsPerState[STATE_ERR] = 1;
    numberOfCompsPerState[STATE_DUP] = 1;
    numberOfCompsPerState[STATE_HAP] = 1;
    numberOfCompsPerState[STATE_COL] = numberOfCollapsedComps;


    int maxNumberOfComps = Int_getMaxValue1DArray(numberOfCompsPerState, numberOfStates);

    // select the coverage of the first region as the baseline
    double medianCoverage = header->regionCoverages[0];
    // calculate region scales against the baseline (first region)
    double *regionScales = Double_construct1DArray(header->numberOfRegions);
    for (int i = 0; i < header->numberOfRegions; i++) {
        regionScales[i] = header->regionCoverages[i] / medianCoverage;
    }

    double **means = Double_construct2DArray(numberOfStates, maxNumberOfComps);
    // initialize mean values and deviate them from the real values randomly
    // by at most a factor of 0.2 which can be either upward or downward
    // This is to make sure the EM algorithm can find the real values even if
    // out inital values are not exact
    if (0.5 < initialDeviation) {
        fprintf(stderr, "[%s] Error: Initial random deviation for the model parameters cannot be greater than 0.5. \n",
                get_timestamp());
        exit(EXIT_FAILURE);
    }
    // ERR_COMP_BINDING_COEF is 0.1 (defined in hmm_utils.h)
    means[STATE_ERR][0] =
            medianCoverage * ERR_COMP_BINDING_COEF * getRandomNumber(1.0 - initialDeviation, 1.0 + initialDeviation);
    means[STATE_DUP][0] = medianCoverage * 0.5 * getRandomNumber(1.0 - initialDeviation, 1.0 + initialDeviation);
    means[STATE_HAP][0] = medianCoverage * 1.0 * getRandomNumber(1.0 - initialDeviation, 1.0 + initialDeviation);
    for (int i = 0; i < numberOfCompsPerState[STATE_COL]; i++) {
        means[STATE_COL][i] =
                means[STATE_HAP][0] * (i + 2) * getRandomNumber(1.0 - initialDeviation, 1.0 + initialDeviation);
    }

    double minHighlyClippedRatio = 1.0; // Msj is not supported now

    HMM *model = HMM_construct(numberOfStates,
                               numberOfRegions,
                               numberOfCompsPerState,
                               means,
                               regionScales,
                               maxHighMapqRatio,
                               minHighlyClippedRatio,
                               NULL,
                               modelType,
                               alphaMatrix,
                               true);
    return model;
}

void runHMMFlagger(ChunksCreator *chunksCreator,
                   HMM **modelPtr,
                   int numberOfIterations,
                   double convergenceTol,
                   char *outputDir,
                   int threads,
                   bool writeParameterStatsPerIteration,
                   bool writeBenchmarkingStatsPerIteration,
                   stList *labelNamesWithUnknown,
                   char *binArrayFilePath,
                   double overlapRatioThreshold,
                   bool acceleration) {

    HMM *model = *modelPtr;

    char loglikelihoodPath[2000];
    sprintf(loglikelihoodPath, "%s/loglikelihood.tsv", outputDir);
    fprintf(stderr, "[%s] Opening file for writing loglikelihood value per EM iteration...\n", get_timestamp());
    FILE *loglikelihoodTsvFile = fopen(loglikelihoodPath, "w+");
    // write header for loglikelihood tsv
    fprintf(loglikelihoodTsvFile, "#Iteration\tEffective_Iteration\tLoglikelihood\n");

    char suffix[200];
    stList *chunks = chunksCreator->chunks;
    int numberOfChunks = stList_length(chunks);


    fprintf(stderr, "[%s] Creating EM arrays (for %d chunks) ...\n",
            get_timestamp(),
            numberOfChunks);

    stList *emPerChunk = stList_construct3(0, EM_destruct);
    for (int chunkIndex = 0; chunkIndex < numberOfChunks; chunkIndex++) {
        Chunk *chunk = stList_get(chunks, chunkIndex);
        EM *em = EM_construct(chunk->coverageInfoSeq, chunk->coverageInfoSeqLen, model);
        stList_append(emPerChunk, em);
    }

    // write initial parameter values in a tsv file
    writeParameterStats(model, outputDir, "initial");

    int iter = 1;
    bool converged = false;
    while (iter <= numberOfIterations && converged == false) {
        fprintf(stderr, "[%s] [Iteration %s = %d] Running EM jobs for %d chunks (with %d threads) ...\n",
                get_timestamp(),
                acceleration ? "accelerated" : "",
                iter,
                numberOfChunks,
                threads);
        EM_runOneIterationForList(emPerChunk, model, threads);


        fprintf(stderr, "[%s] [Iteration %s = %d] EM jobs are all finished.\n",
                get_timestamp(),
                acceleration ? "accelerated" : "",
                iter);

        // chunks now definitely contain prediction labels
        chunksCreator->header->isPredictionAvailable = true;
        chunksCreator->header->numberOfLabels = 4;

        // save loglikelihood
        fprintf(loglikelihoodTsvFile, "%d\t%d\t%.4f\n", iter - 1, acceleration ? 3 * (iter - 1) : iter - 1,
                model->loglikelihood);

        // write benchmarking stats
        if (writeBenchmarkingStatsPerIteration || iter == 1) {
            if (iter == 1) {
                strcpy(suffix, "initial");
            } else {
                // (iter - 1) because this prediction is for the parameters updated after the previous iteration
                if (acceleration) {
                    sprintf(suffix, "iteration_accelerated_%d", iter - 1);
                } else {
                    sprintf(suffix, "iteration_%d", iter - 1);
                }
            }
            writeBenchmarkingStats(chunksCreator,
                                   outputDir,
                                   suffix,
                                   labelNamesWithUnknown,
                                   binArrayFilePath,
                                   overlapRatioThreshold,
                                   threads);
        }

        // if acceleration is true it will run two more EM update per iteration
        if (acceleration) {
            fprintf(stderr, "[%s] [Iteration %s = %d] Running SQUAREM acceleration.\n", get_timestamp(),
                    acceleration ? "accelerated" : "",
                    iter);
            SquareAccelerator *accelerator = SquareAccelerator_construct();
            // set model 0
            SquareAccelerator_setModel0(accelerator, model);
            // compute and set model 1
            HMM_estimateParameters(model, convergenceTol);
            SquareAccelerator_setModel1(accelerator, model);
            // compute and set model 2
            HMM_resetEstimators(model);
            EM_runOneIterationForList(emPerChunk, model, threads);
            HMM_estimateParameters(model, convergenceTol);
            SquareAccelerator_setModel2(accelerator, model);
            // get model prime and run an additional round of EM on it
            HMM *modelPrime = SquareAccelerator_getModelPrime(accelerator, emPerChunk, threads);
            // run EM on model prime
            HMM_resetEstimators(modelPrime);
            EM_runOneIterationForList(emPerChunk, modelPrime, threads);
            // destruct old model and replace with a copy of the final model after acceleration
            HMM_destruct(model);
            model = HMM_copy(modelPrime);
            *modelPtr = model;
            // update model object in em objects
            for (int chunkIndex = 0; chunkIndex < stList_length(emPerChunk); chunkIndex++) {
                EM *em = stList_get(emPerChunk, chunkIndex);
                EM_renewParametersAndEstimatorsFromModel(em, model);
            }
            // destruct accelerator
            SquareAccelerator_destruct(accelerator);
            fprintf(stderr, "[%s] [Iteration %s = %d] Finished SQUAREM acceleration.\n", get_timestamp(),
                    acceleration ? "accelerated" : "",
                    iter);
        }

        // update parameters
        converged = HMM_estimateParameters(model, convergenceTol);
        fprintf(stderr, "[%s] [Iteration %s = %d] Parameters are estimated and updated.\n",
                get_timestamp(),
                acceleration ? "accelerated" : "",
                iter);

        HMM_resetEstimators(model);
        fprintf(stderr, "[%s] [Iteration %s = %d] Model parameter estimators are reset.\n",
                get_timestamp(),
                acceleration ? "accelerated" : "",
                iter);

        if (writeParameterStatsPerIteration) {
            fprintf(stderr,
                    "[%s] [Iteration %s = %d] Writing parameter values into tsv file.\n",
                    acceleration ? "accelerated" : "",
                    get_timestamp(),
                    iter);
            if (acceleration) {
                sprintf(suffix, "iteration_accelerated_%d", iter);
            } else {
                sprintf(suffix, "iteration_%d", iter);
            }
            writeParameterStats(model, outputDir, suffix);
        }
        iter += 1;
    }

    if (converged == true) {
        fprintf(stderr, "[%s] Parameters converged after %d iterations (tol=%.2e)\n",
                get_timestamp(),
                iter - 1,
                convergenceTol);
    } else {
        fprintf(stderr,
                "[%s] Parameter estimation stopped (not yet converged based on the given tolerance) after %d iterations (tol=%.2e)\n",
                get_timestamp(),
                iter - 1,
                convergenceTol);
    }

    fprintf(stderr, "[%s] [Final Inference] Running EM jobs for %d chunks (with %d threads) ...\n",
            get_timestamp(),
            numberOfChunks,
            threads);
    EM_runOneIterationForList(emPerChunk, model, threads);
    fprintf(stderr, "[%s] [Final Inference] EM jobs are all finished.\n", get_timestamp());

    fprintf(loglikelihoodTsvFile, "%d\t%d\t%.4f\n", iter - 1, acceleration ? 3 * (iter - 1) : iter - 1,
            model->loglikelihood);


    sprintf(suffix, "final");
    fprintf(stderr,
            "[%s] [Final Inference] Writing final parameter and benchmarking stats into tsv file.\n",
            get_timestamp());
    writeParameterStats(model, outputDir, suffix);
    writeBenchmarkingStats(chunksCreator,
                           outputDir,
                           suffix,
                           labelNamesWithUnknown,
                           binArrayFilePath,
                           overlapRatioThreshold,
                           threads);

    stList_destruct(emPerChunk);
    fclose(loglikelihoodTsvFile);
}

// input can be NULL
MatrixDouble *getAlphaMatrix(char *alphaTsvPath) {
    if (alphaTsvPath == NULL) {
        MatrixDouble *alpha = MatrixDouble_construct0(NUMBER_OF_STATES - 1, NUMBER_OF_STATES - 1);
        MatrixDouble_setValue(alpha, 0.0);
        return alpha;
    }

    bool skipFirstLine = false;
    // -1 because of ignoring MSJ
    MatrixDouble *alpha = MatrixDouble_parseFromFile(alphaTsvPath, NUMBER_OF_STATES - 1, NUMBER_OF_STATES - 1,
                                                     skipFirstLine);
    // check all alpha values are between 0 and 1
    for (int i = 0; i < alpha->dim1; i++) {
        for (int j = 0; j < alpha->dim2; j++) {
            if (1.0 < alpha->data[i][j] || alpha->data[i][j] < 0.0) {
                fprintf(stderr, "[%s] Error: There is at least one alpha value in '%s' not between 0 and 1. \n%s\n",
                        get_timestamp(),
                        MatrixDouble_toString(alpha)
                );
                exit(EXIT_FAILURE);
            }
        }
    }
    return alpha;
}


static struct option long_options[] =
        {
                {"input",                              required_argument, NULL, 'i'},
                {"iterations",                         required_argument, NULL, 'n'},
                {"convergenceTol",                     required_argument, NULL, 't'},
                {"modelType",                          required_argument, NULL, 'm'},
                {"maxHighMapqRatio",                   required_argument, NULL, 'q'},
                {"chunkLen",                           required_argument, NULL, 'C'},
                {"windowLen",                          required_argument, NULL, 'W'},
                {"contigsList",                        required_argument, NULL, 'c'},
                {"threads",                            required_argument, NULL, '@'},
                {"collapsedComps",                     required_argument, NULL, 'p'},
                {"alphaTsv",                           required_argument, NULL, 'A'},
                {"binArrayFile",                       required_argument, NULL, 'a'},
                {"writeParameterStatsPerIteration",    no_argument,       NULL, 'w'},
                {"writeBenchmarkingStatsPerIteration", no_argument,       NULL, 'k'},
                {"outputDir",                          required_argument, NULL, 'o'},
                {"overlapRatioThreshold",              required_argument, NULL, 'v'},
                {"labelNames",                         required_argument, NULL, 'l'},
                {"initialRandomDev",                   required_argument, NULL, 'D'},
                {"trackName",                          required_argument, NULL, 'N'},
                {"dumpBin",                            no_argument,       NULL, 'B'},
                {"accelerate",                         no_argument,       NULL, 's'},
                {"minimumLengths",                     required_argument, NULL, 'M'},
                {NULL,                                 0,                 NULL, 0}
        };


int main(int argc, char *argv[]) {
    int c;
    char *trackName = copyString("final_flagger");
    char *inputPath = NULL;
    char *alphaTsvPath = NULL;
    char *contigListPath = NULL;
    char *binArrayFilePath = NULL;
    stList *labelNamesWithUnknown = NULL;
    stList *contigList = NULL;
    int numberOfIterations = 100;
    double convergenceTol = 0.001;
    double maxHighMapqRatio = 0.25;
    int numberOfCollapsedComps = -1;
    char *outputDir = NULL;
    bool writeParameterStatsPerIteration = false;
    bool writeBenchmarkingStatsPerIteration = false;
    ModelType modelType = MODEL_GAUSSIAN;
    double overlapRatioThreshold = 0.4;
    double initialRandomDeviation = 0.0;
    int chunkCanonicalLen = 20000000; //20Mb
    int windowLen = 100;
    int threads = 4;
    bool dumpBin = false;
    bool acceleration = false;
    int *minLenPerState = malloc(4 * sizeof(int));
    minLenPerState[0] = 0;
    minLenPerState[1] = 0;
    minLenPerState[2] = 0;
    minLenPerState[3] = 0;
    char *program;
    (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
    while (~(c = getopt_long(argc, argv, "i:n:t:m:q:C:W:c:@:p:A:a:wko:v:l:D:BN:M:s", long_options, NULL))) {
        switch (c) {
            case 'i':
                inputPath = optarg;
                break;
            case 'n':
                numberOfIterations = atoi(optarg);
                break;
            case 'B':
                dumpBin = true;
                break;
            case 'N':
                trackName = optarg;
                break;
            case 't':
                convergenceTol = atof(optarg);
                break;
            case 'm':
                modelType = getModelTypeFromString(optarg);
                break;
            case 'a':
                binArrayFilePath = optarg;
                break;
            case 'c':
                contigListPath = optarg;
                break;
            case 'p':
                numberOfCollapsedComps = atoi(optarg);
                break;
            case '@':
                threads = atoi(optarg);
                break;
            case 'A':
                alphaTsvPath = optarg;
                break;
            case 'C':
                chunkCanonicalLen = atoi(optarg);
                break;
            case 'W':
                windowLen = atoi(optarg);
                break;
            case 'w':
                writeParameterStatsPerIteration = true;
                break;
            case 'k':
                writeBenchmarkingStatsPerIteration = true;
                break;
            case 'o':
                outputDir = optarg;
                break;
            case 'D':
                initialRandomDeviation = atof(optarg);
                break;
            case 'l':
                labelNamesWithUnknown = Splitter_getStringList(optarg, ',');
                stList_append(labelNamesWithUnknown, copyString("Unk"));
                break;
            case 'v':
                overlapRatioThreshold = atof(optarg);
                break;
            case 'q':
                maxHighMapqRatio = atof(optarg);
                break;
            case 's':
                acceleration = true;
                break;
            case 'M':
                int arraySize = 0;
                int *minLenPerStateTemp = Splitter_getIntArray(optarg, ',', &arraySize);
                if (arraySize != 3) {
                    fprintf(stderr,
                            "[%s] Error: --minimumLengths should contain only 3 tab-delimited positive integers.\n",
                            get_timestamp());
                    exit(EXIT_FAILURE);
                }
                minLenPerState[0] = minLenPerStateTemp[0];
                minLenPerState[1] = minLenPerStateTemp[1];
                minLenPerState[3] = minLenPerStateTemp[2];
                free(minLenPerStateTemp);
                break;
            default:
                if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
            help:
                fprintf(stderr, "\nUsage: %s  -i <INPUT_FILE> -o <OUTPUT_DIR> \n", program);
                fprintf(stderr,
                        "Options:\n");
                fprintf(stderr,
                        "         --input, -i\n"
                        "                           Path to the input file (cov/cov.gz or binary output of create_bin_chunks) \n");
                fprintf(stderr,
                        "         --outputDir, -o\n"
                        "                           Directory for saving output files.\n");
                fprintf(stderr,
                        "         --modelType, -m\n"
                        "                           Model type can be either 'gaussian', 'negative_binomial', or \n"
                        "                           'trunc_exp_gaussian' [Default = 'gaussian']\n");
                fprintf(stderr,
                        "         --trackName, -N\n"
                        "                           The track name that will appear in the final BED.[Default = 'final_flagger']\n");
                fprintf(stderr,
                        "         --chunkLen, -C\n"
                        "                           Chunk length. Each chunk is the length of the genome (in bases) for which \n"
                        "                           forward/backward algorithm will be performed once in each iteration. Splitting\n"
                        "                           genome into chunks is mainly for enabling multi-threading and running \n"
                        "                           this algorithm for multiple parts of the genome simultaneously. \n"
                        "                           [Default = 20000000 (20Mb)]\n");
                fprintf(stderr,
                        "         --windowLen, -W\n"
                        "                           Window length. Coverage information will be averaged along each window \n"
                        "                           (non-overlapping) and the average value of each window will be regarded \n"
                        "                           as one observation for the EM algorithm. Having larger windows will lead \n"
                        "                           to lower resolution but faster runtime. (default = 100)\n");
                fprintf(stderr,
                        "         --labelNames, -l\n"
                        "                           (Optional) (For naming rows in the final benchmarking tsv file) A \n"
                        "                           comma-delimited string of label names (for example \n"
                        "                           'Err,Dup,Hap,Col'). It should match the number of labels in the \n"
                        "                           header of the input file.[default: none]\n");
                fprintf(stderr,
                        "         --iterations, -n\n"
                        "                           Maximum number of iterations [Default = 100]\n");
                fprintf(stderr,
                        "         --contigsList, -c\n"
                        "                           (Optional) Path to a file with a list of contig names to include (one \n"
                        "                           contig name per line) [Optional]\n");
                fprintf(stderr,
                        "         --convergenceTol, -t\n"
                        "                           Convergence tolerance. The EM iteration will stop once the difference \n"
                        "                           between all model parameter values in two consecutive iterations is \n"
                        "                           less than this value. [Default = 0.001]\n");
                fprintf(stderr,
                        "         --maxHighMapqRatio, -q\n"
                        "                           Maximum ratio of high mapq coverage for duplicated component\n");
                fprintf(stderr,
                        "         --alphaTsv, -A\n"
                        "                           (Optional) The dependency factors of the current emission density\n"
                        "                           to the previous emission. This parameter is a tsv file with 4 rows\n"
                        "                           and 4 columns with no header line. All numbers should be between \n"
                        "                           0 and 1. \n"
                        "                           [Default = all alpha factors set to 0]\n");
                fprintf(stderr,
                        "         --collapsedComps, -p\n"
                        "                           (Optional) Force the number of components of the collapsed state\n"
                        "                           to this positive integer number (Only for testing purposes) \n"
                        "                           [Default = it will automatically detect the best number of \n"
                        "                           components by taking the maximum observed coverage.]\n");
                fprintf(stderr,
                        "         --writeParameterStatsPerIteration, -w\n"
                        "                           (Optional) Write emission, transition statistics per \n"
                        "                           each iteration (For only investigating how the model \n"
                        "                           is improved over EM iterations) [Default = disabled].\n");
                fprintf(stderr,
                        "         --writeBenchmarkingStatsPerIteration, -k\n"
                        "                           (Optional) Write benchmarking (precision/recall/f1score)\n"
                        "                           statistics per each iteration (For only investigating how \n"
                        "                           the model is improved over EM iterations) [Default = disabled].\n");

                fprintf(stderr,
                        "         --binArrayFile, -b\n"
                        "                           (Optional) A tsv file (tab-delimited) that contains bin arrays \n"
                        "                           for stratifying results by event size. Bin intervals can have overlap.\n"
                        "                           It should contain three columns. \n"
                        "                           1st column is the closed start of the bin and the 2nd \n"
                        "                           column is the open end. The 3rd column has a name for each bin. \n"
                        "                           For example one row can be '0\t100\t[0-100)'\n"
                        "                           If no file is passed it will consider one large bin as the default value.\n"
                        "                           (Default = [0,1e9) with the name 'ALL_SIZES')\n");
                fprintf(stderr,
                        "         --overlapRatioThreshold, -v\n"
                        "                           Minimum overlap ratio in calculating overlap-based metrics for \n"
                        "                           considering a hit between a ref label (for example truth label for\n"
                        "                           recall) and query label (for example prediction label for recall) \n"
                        "                           [default: 0.4]\n");
                fprintf(stderr,
                        "         --initialRandomDev, -D\n"
                        "                           Randomly deviate the initial mean values for EM algorithm.\n"
                        "                           This is only for experimenting how much HMM-Flagger is tolerant to\n"
                        "                           starting with approximate values. It should be greater than or equal\n"
                        "                           to 0 and less than 0.5 . [default: 0.0]\n");
                fprintf(stderr,
                        "         --threads, -@\n"
                        "                           Number of threads [default: 4]\n");
                fprintf(stderr,
                        "         --dumpBin, -B\n"
                        "                           Dump chunks in binary format in the output dir (it will make \n"
                        "                           later runs faster by skipping chunk creation part if using the\n"
                        "                           same bin file) [default: disabled]\n");
                fprintf(stderr,
                        "         --accelerate, -s\n"
                        "                           Run EM in accelerated mode. Each iteration will run 3 rounds \n"
                        "                           of EM internally so each iteration will take longer in the \n"
                        "                           acceleration mode but it will have a boosted update for \n"
                        "                           parameters. The whole run should converge faster \n"
                        "                           [default: disabled]\n");
                fprintf(stderr,
                        "         --minimumLengths, -M\n"
                        "                           Comma-delimited list of minimum lengths for converting \n"
                        "                           non-Hap short blocks into Hap blocks. Given numbers should be \n"
                        "                           related to states Err, Dup and Col respectively. \n"
                        "                           [default: '0,0,0']\n");
                return 1;
        }
    }

    double realtimeStart = System_getRealTimePoint();

    if (inputPath == NULL) {
        fprintf(stderr, "[%s] Error: Input path cannot be NULL.\n", get_timestamp());
        exit(EXIT_FAILURE);
    }
    if (modelType == MODEL_UNDEFINED) {
        fprintf(stderr,
                "[%s] Error: Model type is not defined. Specify the model type with --model (-m) argument.\n",
                get_timestamp());
        exit(EXIT_FAILURE);
    }
    if (convergenceTol <= 0.0 || convergenceTol > 1.0) {
        fprintf(stderr, "[%s] Error: convergence tol = %2.f should be between 0 and 1.\n",
                get_timestamp(),
                convergenceTol);
        exit(EXIT_FAILURE);
    }
    if (outputDir == NULL) {
        fprintf(stderr, "[%s] Error: --outputDir, -o should be specified.\n", get_timestamp());
        exit(EXIT_FAILURE);
    } else if (folder_exists(outputDir) == false) {
        fprintf(stderr, "[%s] Error: Output directory %s does not exist!\n", get_timestamp(), outputDir);
        exit(EXIT_FAILURE);
    }
    if (contigListPath != NULL) {
        contigList = Splitter_parseLinesIntoList(contigListPath);
        fprintf(stderr,
                "[%s] Parsed the list of contig names to include for this analysis. (Number of contigs = %d)",
                get_timestamp(),
                stList_length(contigList));
    }

    // check input file extensions
    char *inputExtension = extractFileExtension(inputPath);
    if (strcmp(inputExtension, "cov") != 0 &&
        strcmp(inputExtension, "cov.gz") != 0 &&
        strcmp(inputExtension, "bin") != 0) {
        fprintf(stderr,
                "[%s] Error: input file should either cov/cov.gz or a binary file made with create_bin_chunks.\n",
                get_timestamp());
        exit(EXIT_FAILURE);
    }
    free(inputExtension);

    // 1. get chunks and subset to contigs if given
    fprintf(stderr, "[%s] Parsing/Creating coverage chunks. \n", get_timestamp());

    ChunksCreator *chunksCreator = getChunksCreator(inputPath,
                                                    chunkCanonicalLen,
                                                    windowLen,
                                                    threads,
                                                    contigList);

    if (dumpBin) {
        char binPath[1000];
        sprintf(binPath, "%s/chunks.c_%d.w_%d.bin", outputDir, chunksCreator->chunkCanonicalLen,
                chunksCreator->windowLen);
        fprintf(stderr, "[%s] Writing bin file into %s . \n", get_timestamp(), binPath);
        ChunksCreator_writeChunksIntoBinaryFile(chunksCreator, binPath);
    }

    int numberOfChunks = ChunksCreator_getTotalNumberOfChunks(chunksCreator);
    int64_t totalLengthOfChunks = ChunksCreator_getTotalLength(chunksCreator);
    fprintf(stderr, "[%s] %d chunks are parsed covering total length of %ld bases. \n",
            get_timestamp(),
            numberOfChunks,
            totalLengthOfChunks);

    // 2. get maximum coverage and set number of collapsed comps
    if (numberOfCollapsedComps == -1) {
        fprintf(stderr, "[%s] Determining the number of components for the 'collapsed' state. \n", get_timestamp());
        numberOfCollapsedComps = getBestNumberOfCollapsedComps(chunksCreator);
        // make sure the number of components is not too large or too small
        numberOfCollapsedComps = numberOfCollapsedComps < 2 ? 2 : numberOfCollapsedComps;
        numberOfCollapsedComps = numberOfCollapsedComps > 10 ? 10 : numberOfCollapsedComps;
        fprintf(stderr,
                "[%s] The number of collapsed components (n=%d) is determined and adjusted automatically by taking the maximum observed coverage. \n",
                get_timestamp(),
                numberOfCollapsedComps);
    } else {
        fprintf(stderr, "[%s] The number of components for the 'collapsed' state is set by the program argument %d. \n",
                get_timestamp(),
                numberOfCollapsedComps);
    }

    // 3. create a model
    fprintf(stderr, "[%s] Creating HMM model. \n", get_timestamp());

    MatrixDouble *alphaMatrix = getAlphaMatrix(alphaTsvPath);
    HMM *model = createModel(modelType,
                             numberOfCollapsedComps,
                             chunksCreator->header,
                             alphaMatrix,
                             maxHighMapqRatio,
                             initialRandomDeviation);

    // 4. run EM for estimating parameters
    fprintf(stderr, "[%s] Running EM for estimating parameters. \n", get_timestamp());

    runHMMFlagger(chunksCreator,
                  &model,
                  numberOfIterations,
                  convergenceTol,
                  outputDir,
                  threads,
                  writeParameterStatsPerIteration,
                  writeBenchmarkingStatsPerIteration,
                  labelNamesWithUnknown,
                  binArrayFilePath,
                  overlapRatioThreshold,
                  acceleration);


    // 5. write final BED
    fprintf(stderr, "[%s] Writing final BED file. \n", get_timestamp());

    char outputBEDPath[2000];
    sprintf(outputBEDPath, "%s/final_flagger_prediction.bed", outputDir);
    ChunksCreator_writePredictionIntoFinalBED(chunksCreator, outputBEDPath, trackName, minLenPerState);

    ChunksCreator_destruct(chunksCreator);
    HMM_destruct(model);

    fprintf(stderr, "[%s] Done! \n", get_timestamp());

    double realtime = System_getRealTimePoint() - realtimeStart;
    double cputime = System_getCpuTime();
    double rssgb = System_getPeakRSSInGB();
    double usage = System_getCpuUsage(cputime, realtime);
    // copied from https://github.com/chhylp123/hifiasm/blob/70fd9a0b1fea45e442eb5f331922ea91ef4f71ae/main.cpp#L73
    fprintf(stderr, "Real time:  %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB; CPU usage: %.1f\%\n", realtime, cputime,
            rssgb, usage * 100.0);
}
