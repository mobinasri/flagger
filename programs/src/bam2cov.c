//
// Created by mobin on 9/22/23.
//

#include <getopt.h>
#include "sam.h"
#include "faidx.h"
#include <time.h>
#include "bgzf.h"
#include <regex.h>
#include "sonLib.h"
#include <assert.h>
#include <math.h>
#include <float.h>
#include "cigar_it.h"
#include "common.h"
#include "vcf.h"
#include "ptBlock.h"
#include "ptAlignment.h"
#include <time.h>
#include <string.h>
#include "ptBlock.h"
#include "ptAlignment.h"
#include "tpool.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include "cJSON.h"
#include "stdlib.h"
#include "bias_detector.h"


static struct option long_options[] =
        {
                {"bam",                     required_argument, NULL, 'i'},
                {"annotationJson",          required_argument, NULL, 'j'},
                {"mapqThreshold",           required_argument, NULL, 'm'},
                {"clipRatioThreshold",      required_argument, NULL, 'c'},
                {"threads",                 required_argument, NULL, 't'},
                {"restrictBiasAnnotations", required_argument, NULL, 'r'},
                {"baselineAnnotation",      required_argument, NULL, 'b'},
                {"covDiffThreshold",        required_argument, NULL, 'd'},
                {"minBiasCoverage",         required_argument, NULL, 'g'},
                {"minBiasLength",           required_argument, NULL, 'l'},
                {"output",                  required_argument, NULL, 'o'},
                {"format",                  required_argument, NULL, 'f'},
                {"runBiasDetection",        no_argument,       NULL, 'u'},
                {NULL,                      0,                 NULL, 0}
        };


int main(int argc, char *argv[]) {
    int c;
    int mapqThreshold = 20;
    double clipRatioThreshold = 0.1;
    int threads = 4;
    char *bamPath;
    char *outPath;
    char *jsonPath;
    char *program;
    char *restrictBiasAnnotationsPath = NULL;
    char *baselineAnnotationName = copyString("no_annotation");
    double covDiffNormalizedThreshold = 0.15;
    int minBiasCoverage = 4;
    int minBiasLength = 100e3;
    bool runBiasDetection = false;
    char *format = copyString("all");
    (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);

    while (~(c = getopt_long(argc, argv, "i:t:j:m:r:f:o:g:c:b:d:g:a:uh", long_options, NULL))) {
        switch (c) {
            case 'i':
                bamPath = optarg;
                break;
            case 't':
                threads = atoi(optarg);
                break;
            case 'o':
                outPath = optarg;
                break;
            case 'j':
                jsonPath = optarg;
                break;
            case 'm':
                mapqThreshold = atoi(optarg);
                break;
            case 'c':
                clipRatioThreshold = atof(optarg);
                break;
            case 'f':
                format = optarg;
                break;
            case 'r':
                restrictBiasAnnotationsPath = optarg;
                break;
            case 'b':
                free(baselineAnnotationName);
                baselineAnnotationName = optarg;
                break;
            case 'd':
                covDiffNormalizedThreshold = atof(optarg);
                break;
            case 'g':
                minBiasCoverage = atoi(optarg);
                break;
            case 'l':
                minBiasLength = atoi(optarg);
                break;
            case 'u':
                runBiasDetection = true;
                break;
            default:
                if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
            help:
                fprintf(stderr, "\nUsage: %s  -i <BAM_FILE> -j <JSON_FILE> -t <THREADS> -o <OUT_BED_OR_COV_FILE> \n",
                        program);
                fprintf(stderr, "Options:\n");
                fprintf(stderr, "         -i, --bam                            Input bam file (should be indexed)\n");
                fprintf(stderr,
                        "         -j, --annotationJson                 Json file for the annotation bed files. At least one BED file should be in this json and it can be a BED file covering the whole genome/assembly [maximum 58 bed files can be given and the keys can be any number between 1-58 for example {\"1\":\"/path/to/1.bed\", \"2\":\"/path/to/2.bed\"}]\n");
                fprintf(stderr,
                        "         -m, --mapqThreshold                  Minimum mapq for the measuring the coverage of the alignments with high mapq [Default = 20]\n");
                fprintf(stderr,
                        "         -c, --clipRatioThreshold             Minimum clipping ratio for the measuring the coverage of the highly clipped alignments [Default = 0.1]\n");
                fprintf(stderr,
                        "         -f, --format                         If this parameter is enabled and output file extension is either cov or cov.gz then the output will be formatted based on the value of this parameter. options: [\"all\", \"only_total\", \"only_high_mapq\"][Default = \"all\"]\n");
                fprintf(stderr, "         -t, --threads                        Number of threads [default: 4]\n");
                fprintf(stderr,
                        "         -o, --output                         Output path [output file can be either cov/cov.gz/bed/bed.gz]\n");
                fprintf(stderr,
                        "         -r, --restrictBiasAnnotations        Path to a text file that contains one annotation name per line. Only these annotations will be considered for detecting coverage biases. If no file is provided it will run bias detection on all annotations in the json file.\n");
                fprintf(stderr,
                        "         -b, --baselineAnnotation             [coverage bias detection] name of the baseline annotation. If no baseline name is provided it will consider all blocks with no given annotation as baseline.\n");
                fprintf(stderr,
                        "         -d, --covDiffThreshold               Threshold for reporting an annotation as biased or not  (It is being applied on the coverage deviation normalized by the baseline coverage) [default:0.15]\n");
                fprintf(stderr,
                        "         -g, --minBiasCoverage                [coverage bias detection] the most frequent coverage will selected from among the coverages greater than or equal to the value of this parameter (It can be useful for ignoring the regions with no coverage) [default:5]\n");
                fprintf(stderr,
                        "         -l, --minBiasLength                  [coverage bias detection] the minimum length of an annotation to be considered biased [default:100000 (100k)]\n");
                fprintf(stderr,
                        "         -u, --runBiasDetection               Run coverage bias detection. It will update the number of regions. [Default: disabled]\n");
                return 1;
        }
    }
    double realtimeStart = System_getRealTimePoint();

    //merge and create the final block table
    stHash *blockTable = ptBlock_multi_threaded_coverage_extraction_with_zero_coverage_and_annotation(bamPath,
                                                                                                      jsonPath,
                                                                                                      threads,
                                                                                                      mapqThreshold,
                                                                                                      clipRatioThreshold);


    // make a list of all annotation names
    // index 0 is reserved for blocks with no annotation
    const char *annotationZeroName = "no_annotation";
    stList *annotationNames = parse_annotation_names_and_save_in_stList(jsonPath, annotationZeroName);
    fprintf(stderr,
            "[%s] Number of annotations added to the blocks table (+1 for 'no_annotation' which is a reserved annotation name for blocks with no annotation): %d\n",
            get_timestamp(), stList_length(annotationNames));

    char *extension = extractFileExtension(outPath);

    if (runBiasDetection && jsonPath == NULL) {
        fprintf(stderr,
                "[%s] (Warning) Json file is not provided so for bias detection there is only one annotation ('no_annotation')!\n",
                get_timestamp());
    }


    // an array for mapping annotation to region index
    int *annotationToRegionMap = NULL;
    // the corresponding coverage values will updated later
    int *coveragePerRegion = NULL;
    int numberOfRegions = 0;

    if (runBiasDetection) {

        // make a list of annotations names that can be biased
        stList *annotationNamesToCheck = NULL;
        if (restrictBiasAnnotationsPath != NULL) {
            annotationNamesToCheck = Splitter_parseLinesIntoList(restrictBiasAnnotationsPath);
        } else { // copy all names
            fprintf(stderr,
                    "[%s] --restrictBiasAnnotations is not provided so all annotations will be checked for coverage bias!\n",
                    get_timestamp());
            annotationNamesToCheck = stList_construct3(0, free);
            for (int i = 0; i < stList_length(annotationNames); i++) {
                char *annotationName = stList_get(annotationNames, i);
                stList_append(annotationNamesToCheck, copyString(annotationName));
            }
        }
        BiasDetector *biasDetector = BiasDetector_construct(annotationNames, baselineAnnotationName, minBiasCoverage,
                                                            minBiasLength, covDiffNormalizedThreshold);
        BiasDetector_setStatisticsPerAnnotation(biasDetector, blockTable);

        // make a file name for saving bias detection table
        char *prefix = copyString(outPath);
        prefix[strlen(outPath) - strlen(extension) - 1] = '\0';
        char *tsvPathToWriteTable = malloc(strlen(prefix) + 100);
        sprintf(tsvPathToWriteTable, "%s.bias_detection_table.tsv", prefix);

        // run bias detection
        fprintf(stderr, "[%s] Running bias detection. Results table will be saved in %s\n", get_timestamp(),
                tsvPathToWriteTable);
        BiasDetector_runBiasDetection(biasDetector, annotationNamesToCheck, tsvPathToWriteTable);
        // update mapping from annotation to region
        // and region median coverages
        annotationToRegionMap = Int_copy1DArray(biasDetector->annotationToRegionMap, stList_length(annotationNames));
        coveragePerRegion = Int_copy1DArray(biasDetector->coveragePerRegion, biasDetector->numberOfRegions);
        numberOfRegions = biasDetector->numberOfRegions;

        stList_destruct(annotationNamesToCheck);
        BiasDetector_destruct(biasDetector);
        free(prefix);
        free(tsvPathToWriteTable);
    } else {
        fprintf(stderr,
                "[%s] Bias detection is disabled. Running bias module only for getting the whole genome median coverage.\n",
                get_timestamp());

        BiasDetector *biasDetector = BiasDetector_construct(annotationNames, annotationZeroName, 1, 1,
                                                            covDiffNormalizedThreshold);
        BiasDetector_setStatisticsPerAnnotation(biasDetector, blockTable);

        annotationToRegionMap = Int_construct1DArray(stList_length(annotationNames));
        coveragePerRegion = Int_construct1DArray(1);
        numberOfRegions = 1;

        // index 0 is for 'no_annotation'
        coveragePerRegion[0] = biasDetector->mostFrequentCoveragePerAnnotation[0];
        BiasDetector_destruct(biasDetector);
    }

    // set region indices based on the mapping generated above
    ptBlock_set_region_indices_by_mapping(blockTable, annotationToRegionMap, stList_length(annotationNames));

    // get the table from contig name to contig length
    // it is only used once the output type is either .cov or .cov.gz
    stHash *ctgToLen = ptBlock_get_contig_length_stHash_from_bam(bamPath);

    // create a list of header lines
    int numberOfLabels = 0;
    bool isTruthAvailable = false;
    bool isPredictionAvailable = false;
    stList *headerLines = ptBlock_create_headers(annotationNames,
                                                 coveragePerRegion,
                                                 numberOfRegions,
                                                 numberOfLabels,
                                                 isTruthAvailable,
                                                 isPredictionAvailable);
    // write header and tracks into output file
    ptBlock_write_blocks_per_contig(blockTable, outPath, format, ctgToLen, headerLines);

    free(extension);
    stList_destruct(headerLines);
    stHash_destruct(blockTable);
    stHash_destruct(ctgToLen);
    stList_destruct(annotationNames);
    free(annotationToRegionMap);
    free(coveragePerRegion);
    fprintf(stderr, "[%s] Done.\n", get_timestamp());

    double realtime = System_getRealTimePoint() - realtimeStart;
    double cputime = System_getCpuTime();
    double rssgb = System_getPeakRSSInGB();
    double usage = System_getCpuUsage(cputime, realtime);
    // copied from https://github.com/chhylp123/hifiasm/blob/70fd9a0b1fea45e442eb5f331922ea91ef4f71ae/main.cpp#L73
    fprintf(stderr, "Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB; CPU usage: %.1f\%\n", realtime, cputime,
            rssgb, usage * 100.0);
}
