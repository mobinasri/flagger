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
    int numberOfColumns = 4;
    MatrixDouble* testTable = MatrixDouble_parseFromFile(pathToTestData, numberOfLines, numberOfColumns);
    for(int row=0; row < numberOfLines; row++){
        uint8_t coverage =  (uint8_t) testTable->date[row][0];
        uint8_t coverage_high_mapq =  (uint8_t) testTable->date[row][1];
        uint8_t coverage_high_clip =  (uint8_t) testTable->date[row][2];
        int32_t annotation_flag =  (int32_t) testTable->date[row][3];
        CoverageInfo *coverageInfo = CoverageInfo_construct(annotation_flag,
                                              coverage,
                                              coverage_high_mapq,
                                              coverage_high_clip);
        coverageInfoSeq[row] = coverageInfo;
    }
    *seqLen = numberOfLines;
    return coverageInfoSeq;
}

int test1(char * pathToTestData) {
    // calculate initial mean coverages for all region classes
    int numberOfRegions = 1;
    int numberOfStates = 4; // Err, Dup, Hap, and Col

    // set number of mixture components for each state
    int *numberOfCompsPerState = malloc(numberOfStates * sizeof(int));
    numberOfCompsPerState[STATE_ERR] = 1;
    numberOfCompsPerState[STATE_DUP] = 1;
    numberOfCompsPerState[STATE_HAP] = 1;
    numberOfCompsPerState[STATE_COL] = 10;


    int maxNumberOfComps = maxIntArray(numberOfCompsPerState, numberOfStates);
    double **means = Double_construct2DArray(numberOfStates, maxNumberOfComps);
    means[STATE_ERR][0] = 2.0;
    means[STATE_DUP][0] = 10.0;
    means[STATE_HAP][0] = 20.0;
    for (int i = 0; i < numberOfCompsPerState[STATE_COL]; i++) {
        means[STATE_COL][i] = means[STATE_HAP][0] * (i + 2);
    }

    double *meanScalePerRegion = malloc(numberOfRegions * sizeof(double));
    meanScalePerRegion[0] = 1.0;

    // alpha is zero for all transitions
    MatrixDouble *alpha = MatrixDouble_construct0(numberOfStates + 1, numberOfStates + 1);

    double maxHighMapqRatio = 1.0; // Dup state is not dependent on high mapq coverage in this test
    double minHighlyClippedRatio = 1.0; // Msj is not present in this test so it doesn't matter

    HMM *model = HMM_construct(numberOfStates,
                               numberOfRegions,
                               numberOfCompsPerState,
                               means,
                               meanScalePerRegion,
                               maxHighMapqRatio,
                               minHighlyClippedRatio,
                               NULL,
                               MODEL_NEGATIVE_BINOMIAL,
                               alpha);

    int seqLen = 0;
    CoverageInfo **coverageInfoSeq = CoverageInfo_parseTestData(pathToTestData, &seqLen);
    EM *em = EM_construct(coverageInfoSeq, seqLen, model);
    for(int iter=0; iter < 10; iter++) {
        EM_runForward(em);
        EM_runBackward(em);
        EM_updateEstimators(em);
        EM_estimateParameters(em);
    }

    transitionTsvFile = fopen("transition.tsv","w+");
    HMM_printTransitionMatrixInTsvFormat(model, transitionTsvFile);
    fclose(transitionTsvFile);

    emissionTsvFile = fopen("emission.tsv","w+");
    HMM_printEmissionParametersInTsvFormat(model, emissionTsvFile);
    fclose(emissionTsvFile);

    EM_destruct(em);
    HMM_destruct(model);
}


int main(int argc, char *argv[]) {
    int c;
    char testDataPath[1000];
    char *program;
    (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
    while (~(c = getopt_long(argc, argv, "i:h", long_options, NULL))) {
        switch (c) {
            case 'i':
                strcpy(testDataPath, optarg);
                break;
            default:
                if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
            help:
                fprintf(stderr,
                        "\nUsage: %s\n", program);
                fprintf(stderr,
                        "Options:\n");
                fprintf(stderr,
                        "         --testData, -i         path to the test file\n");
                return 1;
        }
    }
    fprintf(stdout, "Running test 1 ... \n");
    test1(testDataPath);
    fprintf(stdout, "Test 1 Finished \n");
}