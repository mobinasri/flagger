#include <getopt.h>
#include <time.h>
#include "bgzf.h"
#include "sonLib.h"
#include <assert.h>
#include <math.h>
#include <float.h>
#include "common.h"
#include <time.h>
#include <string.h>
#include "ptBlock.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include "stdlib.h"
#include "sonLib.h"
#include "chunk.h"
#include "track_reader.h"
#include "cov_fast_reader.h"


static struct option long_options[] =
        {
                {"input",          required_argument, NULL, 'i'},
                {"fai",            required_argument, NULL, 'f'},
                {"predictionBed",  required_argument, NULL, 'p'},
                {"truthBed",       required_argument, NULL, 't'},
                {"numberOfLabels", required_argument, NULL, 'n'},
                {"output",         required_argument, NULL, 'o'},
                {"threads",         required_argument, NULL, '@'},
                {NULL,             0,                 NULL, 0}
        };


int main(int argc, char *argv[]) {
    int c;
    char *inputPath = NULL;
    char *outputPath = NULL;
    char *faiPath = NULL;
    char *truthPath = NULL;
    char *predictionPath = NULL;
    int numberOfLabels = 4;
    int threads = 4;
    char *program;
    (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
    while (~(c = getopt_long(argc, argv, "i:o:f:w:t:n:p:@:h", long_options, NULL))) {
        switch (c) {
            case 'i':
                inputPath = optarg;
                break;
            case 'o':
                outputPath = optarg;
                break;
            case 'f':
                faiPath = optarg;
                break;
            case 'p':
                predictionPath = optarg;
                break;
            case 't':
                truthPath = optarg;
                break;
            case 'n':
                numberOfLabels = atoi(optarg);
                break;
            case '@':
                threads = atoi(optarg);
                break;
            default:
                if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
            help:
                fprintf(stderr, "\nUsage: %s  -i <INPUT_FILE> -f <FAI> -o <OUTPUT_FILE> \n", program);
                fprintf(stderr, "Options:\n");
                fprintf(stderr,
                        "         --input, -i                  input path (can have formats '.cov', '.cov.gz', '.bed' or '.bed.gz')\n");
                fprintf(stderr, "         --fai, -f                    fai path\n");
                fprintf(stderr,
                        "         --output, -o                 output path (can have formats '.cov', '.cov.gz', '.bed' or '.bed.gz')\n");
                fprintf(stderr,
                        "         --truthBed, -t              path to a truth bed file. 4th column in the bed file should contain the truth integer label (0<= label <= --numberOfLabels). Labels with a value of -1 will be considered as not defined. If no truth bed is provided (but a prediction bed exists) the related column will be set all to -1 since 8th column is reserved for truth labels.\n");
                fprintf(stderr,
                        "         --predictionBed, -p         path to a truth bed file. 4th column in the bed file should contain the prediction integer label (0<= label <= --numberOfLabels). Labels with a value of -1 will be considered as not defined. The prediction labels will appear in the 9th column of the output coverage file.\n");
                fprintf(stderr, "         --numberOfLabels, -n         number of labels [Default = 4]\n");
                fprintf(stderr, "         --threads, -@       number of threads [Default = 4]\n");

                return 1;
        }
    }

    double realtimeStart = System_getRealTimePoint();

    if (truthPath == NULL && predictionPath == NULL) {
        fprintf(stderr, "[%s] Error: At least one of --truthBed or --predictionBed should be provided!\n",
                get_timestamp());
        exit(EXIT_FAILURE);
    }
    char *inputExtension = extractFileExtension(inputPath);
    char *outputExtension = extractFileExtension(outputPath);
    if (strcmp(outputExtension, "cov") != 0 &&
        strcmp(outputExtension, "cov.gz") != 0 &&
        strcmp(outputExtension, "bed") != 0 &&
        strcmp(outputExtension, "bed.gz") != 0) {
        fprintf(stderr, "[%s] Error: output file should either cov/cov.gz/bed/bed.gz  !\n", get_timestamp());
        free(inputExtension);
        free(outputExtension);
        exit(EXIT_FAILURE);
    }
    if (strcmp(inputExtension, "cov") != 0 &&
        strcmp(inputExtension, "cov.gz") != 0 &&
        strcmp(inputExtension, "bed") != 0 &&
        strcmp(inputExtension, "bed.gz") != 0) {
        fprintf(stderr, "[%s] Error: input file should either cov/cov.gz/bed/bed.gz  !\n", get_timestamp());
        free(inputExtension);
        free(outputExtension);
        exit(EXIT_FAILURE);
    }

    fprintf(stderr, "[%s] Parsing %s to create a table of contig lengths.\n", get_timestamp(), faiPath);
    stHash *ctgToLen = ptBlock_get_contig_length_stHash_from_fai(faiPath);

    fprintf(stderr, "[%s] Parsing %s.\n", get_timestamp(), inputPath);
    int chunkLen = 40e6;
    CovFastReader *covFastReader =  CovFastReader_construct(inputPath, chunkLen, threads);
    stHash *blockTable = CovFastReader_getBlockTablePerContig(covFastReader);
    // parse information from header
    CoverageHeader *header = CoverageHeader_construct(inputPath);

    //ptBlock_write_blocks_per_contig(blockTable, "test_.bed", "all", ctgToLen, header);

    fprintf(stderr, "[%s] Parsed blocks : tot_len=%ld, number=%ld\n", get_timestamp(),
            ptBlock_get_total_length_by_rf(blockTable),
            ptBlock_get_total_number(blockTable));

    // parse truth/prediction bed files and extend blocks with the related label blocks
    bool isLabelTruth;
    // add truth labels
    if (truthPath != NULL) {
        isLabelTruth = true;
        stHash *blockTableTruth = ptBlock_parse_inference_label_blocks(truthPath, isLabelTruth);
        // add truth labels to the block table
        ptBlock_extend_block_tables(blockTable, blockTableTruth);
    }
    // add prediction labels
    if (predictionPath != NULL) {
        isLabelTruth = false;
        stHash *blockTablePrediction = ptBlock_parse_inference_label_blocks(predictionPath, isLabelTruth);
        // add prediction labels to the block table
        ptBlock_extend_block_tables(blockTable, blockTablePrediction);
    }

    //sort
    ptBlock_sort_stHash_by_rfs(blockTable);

    //merge and create the final block table
    stHash *finalBlockTable = ptBlock_merge_blocks_per_contig_by_rf_v2(blockTable);

    fprintf(stderr, "[%s] Merged blocks after adding label bed tracks : tot_len=%ld, number=%ld\n", get_timestamp(),
            ptBlock_get_total_length_by_rf(finalBlockTable),
            ptBlock_get_total_number(finalBlockTable));



    if (numberOfLabels != header->numberOfLabels) {
        fprintf(stderr,
                "[%s] Warning: The number of labels in the header of the input file (%d) does not match --numberOfLabels %d. It will be overwritten to %d.\n",
                get_timestamp(), header->numberOfLabels, numberOfLabels, numberOfLabels);
    }
    // create a new header with update attributes
    bool isTruthAvailable = truthPath != NULL || header->isTruthAvailable;
    bool isPredictionAvailable = predictionPath != NULL || header->isPredictionAvailable;

    CoverageHeader *newHeader = CoverageHeader_constructByAttributes(header->annotationNames,
                                                                     header->regionCoverages,
                                                                     header->numberOfRegions,
                                                                     numberOfLabels,
                                                                     isTruthAvailable,
                                                                     isPredictionAvailable);

    fprintf(stderr, "[%s] Writing %s.\n", get_timestamp(), outputPath);
    // write header and tracks into output file
    // all means write all available columns
    ptBlock_write_blocks_per_contig(finalBlockTable, outputPath, "all", ctgToLen, newHeader);

    // free memory
    CoverageHeader_destruct(header);
    CoverageHeader_destruct(newHeader);
    CovFastReader_destruct(covFastReader);
    //stHash_destruct(blockTable);
    stHash_destruct(finalBlockTable);
    stHash_destruct(ctgToLen);
    free(outputExtension);
    free(inputExtension);

    fprintf(stderr, "[%s] Done!\n", get_timestamp());

    // log used time/resources
    double realtime = System_getRealTimePoint() - realtimeStart;
    double cputime = System_getCpuTime();
    double rssgb = System_getPeakRSSInGB();
    double usage = System_getCpuUsage(cputime, realtime);
    // copied from https://github.com/chhylp123/hifiasm/blob/70fd9a0b1fea45e442eb5f331922ea91ef4f71ae/main.cpp#L73
    fprintf(stderr, "Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB; CPU usage: %.1f\%\n", realtime, cputime,
            rssgb, usage * 100.0);

}
