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
        fprintf(stderr, "[%s] Chunks are constructed from cov file.\n", get_timestamp());
    } else if (strcmp(inputExtension, "bin") == 0) {
        chunksCreator = ChunksCreator_constructEmpty();
        ChunkCreator_parseChunksFromBinaryFile(chunksCreator, inputPath);
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
    fprintf(stderr, "writing transition tsv...\n");
    FILE *transitionTsvFile = fopen(path, "w+");
    HMM_printTransitionMatrixInTsvFormat(model, transitionTsvFile);
    fclose(transitionTsvFile);

    sprintf(path, "%s/emission_%s.tsv", outputDir, suffix);
    fprintf(stderr, "writing emission tsv ...\n");
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
                   HMM *model,
                   int numberOfIterations,
                   double convergenceTol,
                   char *outputDir,
                   int threads,
                   bool writeParameterStatsPerIteration,
                   bool writeBenchmarkingStatsPerIteration,
                   stList *labelNamesWithUnknown,
                   char *binArrayFilePath,
                   double overlapRatioThreshold) {

    char suffix[200];
    stList *chunks = chunksCreator->chunks;
    int numberOfChunks = stList_length(chunks);


    fprintf(stderr, "[%s] Creating EM arrays (for %d chunks) ...\n",
            get_timestamp(),
            numberOfChunks);

    stList *emPerChunk = stList_construct3(0, EM_destruct());
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
        fprintf(stderr, "[%s] [Iteration = %d] Running EM jobs for %d chunks (with %d threads) ...\n",
                get_timestamp(),
                iter,
                numberOfChunks,
                threads);
        EM_runOneIterationForList(emPerChunk, threads);

        fprintf(stderr, "[%s] [Iteration = %d] EM jobs are all finished.\n",
                get_timestamp(),
                iter);

        // chunks now definitely contain prediction labels
        chunksCreator->header->isPredictionAvailable = true;
        chunksCreator->header->numberOfLabels = 4;

        // write benchmarking stats
        if (writeBenchmarkingStatsPerIteration || iter == 1) {
            if (iter == 1) {
                strcpy(suffix, "initial");
            } else {
                // (iter - 1) because this prediction is for the parameters updated after the previous iteration
                sprintf(suffix, "iteration_%d", iter - 1);
            }
            writeBenchmarkingStats(chunksCreator,
                                   outputDir,
                                   suffix,
                                   labelNamesWithUnknown,
                                   binArrayFilePath,
                                   overlapRatioThreshold,
                                   threads);
        }

        // update parameters
        converged = HMM_estimateParameters(model, convergenceTol);
        fprintf(stderr, "[%s] [Iteration = %d] Parameters are estimated and updated.\n",
                get_timestamp(),
                iter);

        HMM_resetEstimators(model);
        fprintf(stderr, "[%s] [Iteration = %d] Model parameter estimators are reset.\n",
                get_timestamp(),
                iter);

        if (writeParameterStatsPerIteration) {
            fprintf(stderr,
                    "[%s] [Iteration = %d] Writing parameter values into tsv file.\n",
                    get_timestamp(),
                    iter);
            sprintf(suffix, "iteration_%d", iter);
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
    EM_runOneIterationForList(emPerChunk, threads);
    fprintf(stderr, "[%s] [Final Inference] EM jobs are all finished.\n", get_timestamp());


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
}

// input can be NULL
MatrixDouble *getAlphaMatrix(char *alphaInputString) {

    // -1 because of ignoring MSJ
    MatrixDouble *alpha = MatrixDouble_construct0(NUMBER_OF_STATES - 1, NUMBER_OF_STATES - 1);

    MatrixDouble_setValue(alpha, 0.0);

    //overwrite alpha if it's given as a comma-separated string
    if (alphaInputString != NULL && alphaInputString[0] != '\0') {
        int alphaSize = 0;
        double *alphaInputArray = Splitter_getDoubleArray(alphaInputString, ',', &alphaSize);
        if (alphaSize < 5) {
            fprintf(stderr, "[%s] Error: alpha string '%s' should be comma-separated and has at least 5 numbers \n",
                    get_timestamp(),
                    alphaInputString);
            exit(EXIT_FAILURE);
        }
        // check all alpha values are between 0 and 1
        for (int i = 0; i < alphaSize; i++) {
            if (1.0 < alphaInputArray[i] || alphaInputArray[i] < 0.0) {
                fprintf(stderr, "[%s] Error: There is at least one alpha value in '%s' not between 0 and 1. \n",
                        get_timestamp(),
                        alphaInputString);
                exit(EXIT_FAILURE);
            }
        }
        alpha->data[ALPHA_ERR][ALPHA_ERR] = alphaInputArray[0];
        alpha->data[ALPHA_DUP][ALPHA_DUP] = alphaInputArray[1];
        alpha->data[ALPHA_HAP][ALPHA_HAP] = alphaInputArray[2];
        alpha->data[ALPHA_COL][ALPHA_COL] = alphaInputArray[3];
        alpha->data[ALPHA_TRN][ALPHA_TRN] = alphaInputArray[4];
        free(alphaInputArray);
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
                {"alpha",                              required_argument, NULL, 'A'},
                {"binArrayFile",                       required_argument, NULL, 'a'},
                {"writeParameterStatsPerIteration",    no_argument,       NULL, 'w'},
                {"writeBenchmarkingStatsPerIteration", no_argument,       NULL, 'k'},
                {"outputDir",                          required_argument, NULL, 'o'},
                {"overlapRatioThreshold",              required_argument, NULL, 'v'},
                {"labelNames",                         required_argument, NULL, 'l'},
                {"initialRandomDev",                   required_argument, NULL, 'D'},
                {"trackName",                          required_argument, NULL, 'N'},
                {NULL,                                 0,                 NULL, 0}
        };


int main(int argc, char *argv[]) {
    int c;
    char *trackName = copyString("final_flagger");
    char *inputPath = NULL;
    char *alphaString = NULL;
    char *contigListPath = NULL;
    stList *labelNamesWithUnknown = NULL;
    stList *contigList = NULL;
    int numberOfIterations = 100;
    double convergenceTol = 0.001;
    double maxHighMapqRatio = 0.25;
    int numberOfCollapsedComps = -1;
    char *outputDir = NULL;
    bool writeParameterStatsPerIteration = false;
    bool writeBenchmarkingStatsPerIteration = false;
    ModelType modelType = MODEL_UNDEFINED;
    double overlapRatioThreshold = 0.4;
    double initialRandomDeviation = 0.0;
    int chunkCanonicalLen = 20000000; //20Mb
    int windowLen = 100;
    int threads = 4;
    char *program;
    (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
    while (~(c = getopt_long(argc, argv, "i:n:t:m:q:C:W:c:@:p:A:a:wko:v:l:D:N:", long_options, NULL))) {
        switch (c) {
            case 'i':
                inputPath = optarg;
                break;
            case 'n':
                numberOfIterations = atoi(optarg);
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
                alphaString = optarg;
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
            default:
                if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
            help:
                fprintf(stderr,
                        "\nUsage: %s\n", program);
                fprintf(stderr,
                        "Options:\n");
                fprintf(stderr,
                        "         --input, -i	                Path to the input file (cov/cov.gz or binary output of create_bin_chunks) \n");
                fprintf(stderr,
                        "         --outputDir, -o             Directory for saving output files.\n");
                fprintf(stderr,
                        "         --modelType, -m	            Model type can be either 'gaussian', 'negative_binomial', or \n"
                        "                                     'trunc_exp_gaussian' [Default = not defined]\n");
                fprintf(stderr,
                        "         --trackName, -N	            The track name that will appear in the final BED.[Default = 'final_flagger']\n");
                fprintf(stderr,
                        "         --chunkLen, -C              Chunk length. Each chunk is the length of the genome (in bases) for which \n"
                        "                                     forward/backward algorithm will be performed once in each iteration. It is \n"
                        "                                     mainly for using multi threads while running this algorithm for multiple \n"
                        "                                     parts of the genome simultaneously. [Default = 20000000 (20Mb)]\n");
                fprintf(stderr,
                        "         --windowLen, -W             Window length. Coverage information will be averaged along each window \n"
                        "                                     (non-overlapping) and the average value of each window will be regarded \n"
                        "                                     as one observation for the EM algorithm. Having larger windows will lead \n"
                        "                                     to lower resolution but faster runtime. (default = 100)\n");
                fprintf(stderr,
                        "         --labelNames, -l           (Optional) (For naming rows in the final benchmarking tsv file) A \n"
                        "                                     comma-delimited string of label names (for example \n"
                        "                                     'Err,Dup,Hap,Col'). It should match the number of labels in the \n"
                        "                                     header of the input file.[default: none]\n");
                fprintf(stderr,
                        "         --iterations, -n		    Maximum number of iterations [Default = 100]\n");
                fprintf(stderr,
                        "         --contigsList, -c		    (Optional) Path to a file with a list of contig names to include (one \n"
                        "                                     contig name per line) [Optional]\n");
                fprintf(stderr,
                        "         --convergenceTol, -t		Convergence tolerance. The EM iteration will stop once the difference \n"
                        "                                     between all model parameter values in two consecutive iterations is \n"
                        "                                     less than this value. [Default = 0.001]\n");
                fprintf(stderr,
                        "         --maxHighMapqRatio, -q      Maximum ratio of high mapq coverage for duplicated component\n");
                fprintf(stderr,
                        "         --alpha, -A                 Alpha is the dependency factor of the current emission density\n"
                        "                                     to the previous emission. It should be a comma-separated string\n"
                        "                                     of 5 numbers for these states respectively err,dup,hap,col,trans.\n"
                        "                                     (trans is for transitioning from one state to a different one)\n");
                fprintf(stderr,
                        "         --collapsedComps, -p		(Optional) Force the number of components of the collapsed state\n"
                        "                                     to this positive integer number (Only for testing purposes) \n"
                        "                                     [Default = it will automatically detect the best number of \n"
                        "                                     components by taking the maximum observed coverage.]\n");
                fprintf(stderr,
                        "         --writeParameterStatsPerIteration, -w	    (Optional) Write emission, transition statistics per \n"
                        "                                                     each iteration (For only investigating how the model \n"
                        "                                                     is improved over EM iterations) [Default = disabled].\n");
                fprintf(stderr,
                        "         --writeBenchmarkingStatsPerIteration, -k	(Optional) Write benchmarking (precision/recall/f1score)\n"
                        "                                                     statistics per each iteration (For only investigating how \n"
                        "                                                     the model is improved over EM iterations) [Default = disabled].\n");

                fprintf(stderr,
                        "         -b,--binArrayFile                           (Optional) A tsv file (tab-delimited) that contains bin arrays \n"
                        "                                                     for stratifying results by event size. It should contain three \n"
                        "                                                     columns. 1st column is the closed start of the bin and the 2nd \n"
                        "                                                     column is the open end. The 3rd column has a name for each bin. \n"
                        "                                                     For example one row can be '0\t100\t[0-100). If no file is passed\n"
                        "                                                     it will consider one large bin as the default value. "
                        "                                                     (Default = [0,1e9) with the name 'ALL_SIZES')\n");
                fprintf(stderr,
                        "         -v, --overlapRatioThreshold                 Minimum overlap ratio in calculating overlap-based metrics for \n"
                        "                                                     considering a hit between a ref label (for example truth label for\n"
                        "                                                     recall) and query label (for example prediction label for recall) \n"
                        "                                                     [default: 0.4]\n");
                return 1;
        }
    }

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

    // 2. get maximum coverage and set number of collapsed comps
    if (numberOfCollapsedComps == -1) {
        fprintf(stderr, "[%s] Determining the number of components for the 'collapsed' state. \n", get_timestamp());
        numberOfCollapsedComps = getBestNumberOfCollapsedComps(chunksCreator);
        // make sure the number of components is not too large or too small
        numberOfCollapsedComps = numberOfCollapsedComps < 2 ? 2 : numberOfCollapsedComps;
        numberOfCollapsedComps = numberOfCollapsedComps > 20 ? 20 : numberOfCollapsedComps;
        fprintf(stderr,
                "[%s] The number of collapsed components (n=%d) is determined and adjusted automatically by taking the maximum observed coverage. \n",
                get_timestamp(),
                numberOfCollapsedComps);
    }else{
        fprintf(stderr, "[%s] The number of components for the 'collapsed' state is set by the program argument %d. \n",
                get_timestamp(),
                numberOfCollapsedComps);
    }

    // 3. create a model
    fprintf(stderr, "[%s] Creating HMM model. \n", get_timestamp());

    MatrixDouble *alphaMatrix = getAlphaMatrix(alphaString);
    HMM *model = createModel(modelType,
                             numberOfCollapsedComps,
                             chunsCreator->header,
                             alphaMatrix,
                             maxHighMapqRatio,
                             initialRandomDeviation)

    // 4. run EM for estimating parameters
    fprintf(stderr, "[%s] Running EM for estimating parameters. \n", get_timestamp());

    runHMMFlagger(chunksCreator,
                  model,
                  numberOfIterations,
                  convergenceTol,
                  outputDir,
                  threads,
                  writeParameterStatsPerIteration,
                  writeBenchmarkingStatsPerIteration,
                  labelNamesWithUnknown,
                  binArrayFilePath,
                  overlapRatioThreshold);


    // 5. write final BED
    fprintf(stderr, "[%s] Writing final BED file. \n", get_timestamp());

    char outputBEDPath[2000];
    sprintf(outputBEDPath, "%s/final_flagger_prediction.bed", outputDir);
    ChunksCreator_writePredictionIntoFinalBED(chunksCreator, outputBEDPath, trackName);

    ChunksCreator_destruct(chunksCreator);
    HMM_destruct(model);

    fprintf(stderr, "[%s] Done! \n", get_timestamp());
}
