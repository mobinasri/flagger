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



size_t CoverageInfo_getNumberOfLinesTestData(char* pathToTestData) {
    FILE* f = fopen(pathToTestData, "r");
    size_t len = 0;
    char* line = NULL;
    char* token;
    ssize_t read;
    size_t numberOfLines = 0;
    while ((read = getline(&line, &len, f)) != -1){
        numberOfLines += 1;
    }
    fclose(f);
    return numberOfLines;
}

CoverageInfo **CoverageInfo_parseTestData(char* pathToTestData, int *seqLen){
    int numberOfLines = CoverageInfo_getNumberOfLinesTestData(pathToTestData);
    CoverageInfo **coverageInfoSeq = (CoverageInfo **) malloc(numberOfLines * sizeof(CoverageInfo *));
    int numberOfColumns = 5;
    // skip header
    MatrixDouble* testTable = MatrixDouble_parseFromFile(pathToTestData, numberOfLines - 1, numberOfColumns, true);
    for(int row=0; row < numberOfLines - 1; row++){
        uint8_t coverage =  (uint8_t) testTable->data[row][0];
        uint8_t coverage_high_mapq =  (uint8_t) testTable->data[row][1];
        uint8_t coverage_high_clip =  (uint8_t) testTable->data[row][2];
        int32_t annotation_flag =  (int32_t) testTable->data[row][3];
        CoverageInfo *coverageInfo = CoverageInfo_construct(annotation_flag,
                                              coverage,
                                              coverage_high_mapq,
                                              coverage_high_clip);
        coverageInfoSeq[row] = coverageInfo;
    }
    *seqLen = numberOfLines - 1;
    return coverageInfoSeq;
}


uint8_t * parseStateIndicesFromTestData(char* pathToTestData){
    int numberOfLines = CoverageInfo_getNumberOfLinesTestData(pathToTestData);
    uint8_t *stateIndices = (uint8_t *) malloc(numberOfLines * sizeof(uint8_t));
    int numberOfColumns = 5;
    // skip header
    MatrixDouble* testTable = MatrixDouble_parseFromFile(pathToTestData, numberOfLines - 1, numberOfColumns, true);
    for(int row=0; row < numberOfLines - 1; row++){
        uint8_t stateIndex =  (uint8_t) testTable->data[row][4];
        stateIndices[row] = stateIndex;
    }
    return stateIndices;
}

double getRandomNumber(double start, double end){
	srand(time(NULL));
	return (double) rand() / (double)(RAND_MAX/(end - start)) + start;
}

void writeStats(EM* em, char* outputDir, char *suffix, bool writePosterior){
    char path[100];

    if (writePosterior){
        fprintf(stderr, "writing posterior tsv...\n");
        sprintf(path, "%s/posterior_prediction_%s.tsv", outputDir, suffix);
        FILE* posteriorTsvFile = fopen(path,"w+");
        EM_printPosteriorInTsvFormat(em, posteriorTsvFile);
        fclose(posteriorTsvFile);
    }

    sprintf(path, "%s/transition_%s.tsv", outputDir, suffix);
    fprintf(stderr, "writing transition tsv...\n");
    FILE* transitionTsvFile = fopen(path,"w+");
    HMM_printTransitionMatrixInTsvFormat(em->model, transitionTsvFile);
    fclose(transitionTsvFile);

    sprintf(path, "%s/emission_%s.tsv", outputDir, suffix);
    fprintf(stderr, "writing emission tsv ...\n");
    FILE* emissionTsvFile = fopen(path,"w+");
    HMM_printEmissionParametersInTsvFormat(em->model, emissionTsvFile);
    fclose(emissionTsvFile);
}


HMM *createModel(ModelType modelType, double medianCoverage, int collapsedComps, double *regionScales, int numberOfRegions){
    // calculate initial mean coverages for all region classes
    int numberOfStates = 4; // Err, Dup, Hap, and Col

    // set number of mixture components for each state
    int *numberOfCompsPerState = malloc(numberOfStates * sizeof(int));
    numberOfCompsPerState[STATE_ERR] = 1;
    numberOfCompsPerState[STATE_DUP] = 1;
    numberOfCompsPerState[STATE_HAP] = 1;
    numberOfCompsPerState[STATE_COL] = collapsedComps;


    int maxNumberOfComps = maxIntArray(numberOfCompsPerState, numberOfStates);
    double **means = Double_construct2DArray(numberOfStates, maxNumberOfComps);
    // initialize mean values and deviate them from the real values randomly
    // by at most a factor of 0.2 which can be either upward or downward
    // This is to make sure the EM algorithm can find the real values even if
    // out inital values are not exact
    means[STATE_ERR][0] = medianCoverage * 0.1 * getRandomNumber(0.8, 1.2);
    means[STATE_DUP][0] = medianCoverage * 0.5 * getRandomNumber(0.8, 1.2);
    means[STATE_HAP][0] = medianCoverage * 1.0 * getRandomNumber(0.8, 1.2);
    for (int i = 0; i < numberOfCompsPerState[STATE_COL]; i++) {
        means[STATE_COL][i] = means[STATE_HAP][0] * (i + 2) * getRandomNumber(0.8, 1.2);
    }


    // alpha is zero for all transitions
    MatrixDouble *alpha = MatrixDouble_construct0(numberOfStates + 1, numberOfStates + 1);

    double maxHighMapqRatio = 1.0; // Dup state is not dependent on high mapq coverage in this test
    double minHighlyClippedRatio = 1.0; // Msj is not present in this test so it doesn't matter
                                        //
    fprintf(stderr, "creating model ...\n");
    HMM *model = HMM_construct(numberOfStates,
                               numberOfRegions,
                               numberOfCompsPerState,
                               means,
                               regionScales,
                               maxHighMapqRatio,
                               minHighlyClippedRatio,
                               NULL,
                               modelType,
                               alpha,
                               true);
    return model;

}

int runTest(char *pathToTestData, HMM *model, int numberOfIterations, double convergenceTol, char* outputDir, bool writeStatsPerIteration) {
    int seqLen = 0;
    fprintf(stderr, "parsing test file ...\n");
    CoverageInfo **coverageInfoSeq = CoverageInfo_parseTestData(pathToTestData, &seqLen);
    uint8_t *truthStateIndices = parseStateIndicesFromTestData(pathToTestData); 

    char suffix[200];
    fprintf(stderr, "creating em ...\n");
    EM *em = EM_construct(coverageInfoSeq, seqLen, model);
    writeStats(em, outputDir, "initial", false); // will not write posterior tsv
    int iter = 0;
    bool converged = false;
    while(iter < numberOfIterations && converged == false) {
	fprintf(stderr, "Running forward algorithm (round = %d) ...\n", iter);
        EM_runForward(em);
	fprintf(stderr, "Running backward algorithm (round = %d) ...\n", iter);
        EM_runBackward(em);
	fprintf(stderr, "Updating estimators (round = %d) ...\n", iter);
        EM_updateEstimators(em);
	fprintf(stderr, "Estimating and updating parameters (round = %d) ...\n", iter);
        converged = EM_estimateParameters(em, convergenceTol);
	fprintf(stderr, "Reseting estimators (round = %d) ...\n", iter);
        EM_resetEstimators(em);
	iter += 1;

	if (writeStatsPerIteration){
	    sprintf(suffix, "iter_%d", iter);
	    writeStats(em, outputDir, suffix, false); // will not write posterior tsv
	}
    }
    fprintf(stderr, "Parameters converged after %d iterations (tol=%.2e)\n", iter, convergenceTol);
    fprintf(stderr, "Running forward algorithm (round = inference) ...\n");
    EM_runForward(em);
    fprintf(stderr, "Running backward algorithm(round = inference) ...\n");
    EM_runBackward(em);
    int correct = 0;
    for(int pos=0; pos < seqLen ; pos++){
	    int prediction = EM_getMostProbableState(em, pos);
	    int truth = truthStateIndices[pos];
	    correct += prediction == truth ? 1 : 0;
    }
    fprintf(stderr, "Accuracy = %.2e\n", (double) correct / seqLen);

    writeStats(em, outputDir, "final", true); // write posterior tsv

    EM_destruct(em);
}


static struct option long_options[] =
{
    {"testData",         required_argument, NULL, 'i'},
    {"numberOfIterations",         required_argument, NULL, 'n'},
    {"convergenceTol",         required_argument, NULL, 't'},
    {"modelType",         required_argument, NULL, 'm'},
    {"coverage",         required_argument, NULL, 'c'},
    {"collapsedComps",         required_argument, NULL, 'p'},
    {"writeStatsPerIteration",         no_argument, NULL, 'w'},
    {"outputDir",         required_argument, NULL, 'o'},
    {"regionScales", required_argument, NULL, 's'},
    {"numberOfRegions", required_argument, NULL, 'r'},
    {NULL,               0,                 NULL, 0}
};


int main(int argc, char *argv[]) {
    int c;
    char testDataPath[1000];
    int numberOfIterations = 100;
    double convergenceTol = 0.001;
    double medianCoverage = 0.0;
    int collapsedComps = 4;
    char* outputDir = NULL;
    bool writeStatsPerIteration = false;
    ModelType modelType = MODEL_UNDEFINED;
    double *regionScales;
    char *regionScalesStr="1.0";
    int numberOfRegions=1;
    char *program;
    (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
    while (~(c = getopt_long(argc, argv, "i:n:t:m:c:wn:p:r:s:o:h", long_options, NULL))) {
        switch (c) {
            case 'i':
                strcpy(testDataPath, optarg);
                break;
	    case 'n':
		numberOfIterations = atoi(optarg);
		break;
	    case 't':
		convergenceTol = atof(optarg);
		break;
	    case 'm':
		modelType = getModelTypeFromString(optarg);
		break;
	    case 'c':
		medianCoverage = atof(optarg);
		break;
	    case 'p':
		collapsedComps = atoi(optarg);
		break;
	    case 'w':
		writeStatsPerIteration = true;
		break;
	    case 'o':
		outputDir = optarg;
		break;
	    case 'r':
		numberOfRegions = atoi(optarg);
		break;
	    case 's':
		regionScalesStr = optarg;
		break;
            default:
                if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
            help:
                fprintf(stderr,
                        "\nUsage: %s\n", program);
                fprintf(stderr,
                        "Options:\n");
                fprintf(stderr,
                        "         --testData, -i		path to the test file\n");
		fprintf(stderr,
                        "         --outputDir, -o		directory for saving output files.\n");
		fprintf(stderr,
                        "         --iterations, -n		maximum number of iterations [Default = 100]\n");
		fprintf(stderr,
                        "         --convergenceTol, -t		convergence tolerance [Default = 0.001]\n");
		fprintf(stderr,
                        "         --model, -m			model type can be either 'gaussian', 'negative_binomial', or 'trunc_exp_gaussian' [Default = not defined]\n");
		fprintf(stderr,
			"         --coverage, -c		median coverage for initializing the parameters for the EM algorithm. Note that the inital values will be slightly deviated from the given value by a random factor to make sure the EM algorithm can handle imprecise initial values.\n");
		fprintf(stderr,
			"         --collapsedComps, -p		number of components of the collapsed state [Default = 4]\n");
		fprintf(stderr,
			"         --writeStatsPerIteration, -w	write emission, transition and posterior statisitics per each iteration.\n");
		fprintf(stderr,
                        "         --numberOfRegions, -r		number of regions [Default = 1]\n");
                fprintf(stderr,
                        "         --regionScales, -s		a comma-delimited list of scaling factors for setting initial means of different regions [Default = '1.0'] (As an example with --numberOfRegions 3 --regionScales can be set to '1.0,0.5,1.25')\n");
                return 1;
        }
    }
    if (modelType == MODEL_UNDEFINED){
	    fprintf(stderr, "(Error) Model type is not defined. Specify the model type with --model (-m) argument.\n");
	    exit(EXIT_FAILURE);
    }
    if (medianCoverage <= 0.0){
	    fprintf(stderr, "(Error) --coverage, -c should be speficied. It can only be a positive number.\n");
	    exit(EXIT_FAILURE);
    }
    if (numberOfRegions < 1){
	    fprintf(stderr, "(Error) --numberOfRegions, -r should be a positive number.\n");
            exit(EXIT_FAILURE);
    }
    if (outputDir == NULL){
	    fprintf(stderr, "(Error) --outputDir, -o should be speficied.\n");
            exit(EXIT_FAILURE);
    }
    else if (folder_exists(outputDir) == false){
	   fprintf(stderr, "(Error) Output directory %s does not exist!\n", outputDir);
           exit(EXIT_FAILURE);
    }
    Splitter *splitter = Splitter_construct(regionScalesStr, ',');
    char *token;
    regionScales = malloc(numberOfRegions * sizeof(double));
    int i = 0;
    while ((token = Splitter_getToken(splitter)) != NULL) {
	    regionScales[i] = atof(token);
	    i++;
    }
    Splitter_destruct(splitter);
    fprintf(stdout, "Running test... \n");
    HMM *model = createModel(modelType, medianCoverage, collapsedComps, regionScales, numberOfRegions);
    runTest(testDataPath, model, numberOfIterations, convergenceTol, outputDir, writeStatsPerIteration);
    HMM_destruct(model);
    fprintf(stdout, "Test Finished \n");
}
