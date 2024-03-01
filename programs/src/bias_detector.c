//
// Created by mobin on 1/25/24.
//

#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include "sonLib.h"
#include "common.h"
#include "hmm_utils.h"
#include "tpool.h"
#include <pthread.h>
#include <time.h>



CountData** createCountDataPerAnnotation(stHash *final_block_table, int numberOfAnnotations){
    CountData** countDataPerAnnotation = CountData_construct1DArray(MAX_COVERAGE_VALUE, numberOfAnnotations);
    stHashIterator *it = stHash_getIterator(final_block_table);
    char *contig_name;
    while ((contig_name = stHash_getNext(it)) != NULL) {
        stList* blocks = stHash_search(final_block_table, contig_name);
        for (int i = 0; i < stList_length(blocks); i++) {
            ptBlock *block = stList_get(blocks, i);
            CoverageInfo *covInfo = (CoverageInfo *) block->data;
            int annotationIndex = getFirstIndexWithNonZeroBitFromRight(covInfo->annotation_flag);
	    if (annotationIndex < 0) continue;
            int count = block->rfe - block->rfs + 1;
            int value = covInfo->coverage;
            CountData_increment(countDataPerAnnotation[annotationIndex], value, (double) count);
        }
    }
    return countDataPerAnnotation;
}

int* getMostFrequentCoveragesPerAnnotation(CountData** countDataPerAnnotation, int numberOfAnnotations, int lowest_coverage){
    int *mostFrequentCoverages = malloc(numberOfAnnotations * sizeof(int));
    for (int annotationIndex = 0; annotationIndex < numberOfAnnotations; annotationIndex++) {
        int mostFreqCoverage = CountData_getMostFrequentValue(countDataPerAnnotation[annotationIndex], lowest_coverage, MAX_COVERAGE_VALUE);
        mostFrequentCoverages[annotationIndex] = mostFreqCoverage;
    }
    return mostFrequentCoverages;
}

int* getMaxCountsPerAnnotation(CountData** countDataPerAnnotation, int numberOfAnnotations, int lowest_coverage){
    int *maxCounts = malloc(numberOfAnnotations * sizeof(int));
    for (int annotationIndex = 0; annotationIndex < numberOfAnnotations; annotationIndex++) {
        int maxCount = CountData_getMaxCount(countDataPerAnnotation[annotationIndex], lowest_coverage, MAX_COVERAGE_VALUE);
        maxCounts[annotationIndex] = maxCount;
    }
    return maxCounts;
}


int main(int argc, char *argv[]) {
    int c;
    int min_mapq = 20;
    double min_clipping_ratio = 0.1;
    int threads=4;
    char *bam_path;
    char *json_path;
    char *baseline_annot_name;
    double cov_diff_normalized_threshold = 0.15;
    int lowest_coverage = 5;
    int min_count = 10000; // 10Kb
    char *program;
    (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
    while (~(c = getopt(argc, argv, "i:t:j:b:d:e:m:h"))) {
        switch (c) {
            case 'i':
                bam_path = optarg;
                break;
            case 't':
                threads = atoi(optarg);
                break;
            case 'j':
                json_path = optarg;
                break;
            case 'b':
                baseline_annot_name = optarg;
                break;
            case 'd':
		cov_diff_normalized_threshold = atof(optarg);
		break;
            case 'e':
		lowest_coverage = atoi(optarg);
		break;
	    case 'm':
		min_count = atoi(optarg);
		break;
            default:
                if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
            help:
                fprintf(stderr, "\nUsage: %s  -i <BAM_FILE> -t <THREADS> -j <ANNOTATION_JSON> \n", program);
                fprintf(stderr, "Options:\n");
                fprintf(stderr, "         -i         input bam file (should be indexed)\n");
                fprintf(stderr, "         -j         JSON file for the annotation bed files [maximum 32 files can be given and the keys can any string {\"hsat1\":\"/path/to/hsat1.bed\", \"bsat\":\"/path/to/bsat.bed\"}]\n");
                fprintf(stderr, "         -t         number of threads [default: 4]\n");
                fprintf(stderr, "         -b         name of the baseline annotation\n");
		fprintf(stderr, "         -d         threshold for reporting an annotation as biased or not  (It is being applied on the coverage deviation normalized by the baseline coverage) [default:0.15]\n");
		fprintf(stderr, "         -e         the most frequent coverage will selected from among the coverages greater than or equal to the value of this parameter [default:5]\n");
		fprintf(stderr, "         -m         minimum value for the mode count of an annotation that can be reported as biased [default:10000]\n");
                return 1;
        }
    }
    //merge and create the final block table
    stHash * final_block_table = ptBlock_multi_threaded_coverage_extraction_with_zero_coverage_and_annotation(
            bam_path,
            json_path,
            threads,
            min_mapq,
            min_clipping_ratio
            );
    
    stList *annotationPaths = parse_annotation_paths_and_save_in_stList(json_path);
    stList *annotationNames = parse_annotation_names_and_save_in_stList(json_path);
    int numberOfAnnotations = stList_length(annotationNames);
    int baselineIndex = get_annotation_index(annotationNames, baseline_annot_name);

    CountData** countDataPerAnnotation = createCountDataPerAnnotation(final_block_table, numberOfAnnotations);
    int* mostFrequentCoverages = getMostFrequentCoveragesPerAnnotation(countDataPerAnnotation, numberOfAnnotations, lowest_coverage);
    int* maxCounts = getMaxCountsPerAnnotation(countDataPerAnnotation, numberOfAnnotations, lowest_coverage);
    int baselineCoverage = mostFrequentCoverages[baselineIndex];


    fprintf(stderr, "Printing coverage table in stdout ...\n");    
    fprintf(stdout, "annotation\tstatus\tmost_freq_cov\tcov_diff_normalized\tpath\n");
    for(int annotIndex=0; annotIndex < numberOfAnnotations; annotIndex++){
        char* annotationName = stList_get(annotationNames, annotIndex);
        char* annotationPath = stList_get(annotationPaths, annotIndex);
        double coverageDiffNormalized = ((double) mostFrequentCoverages[annotIndex]  - baselineCoverage) / baselineCoverage;
	double coverageDiffNormalizedAbs = coverageDiffNormalized < 0 ? -1 * coverageDiffNormalized : coverageDiffNormalized;
	int maxCount = maxCounts[annotIndex];
        fprintf(stdout, "%s\t%s\t%d\t%+.3f\t%s\n",
                annotationName,
		(cov_diff_normalized_threshold < coverageDiffNormalizedAbs) && (maxCount > min_count) ? "biased" : "not_biased",
                mostFrequentCoverages[annotIndex],
                coverageDiffNormalized,
		annotationPath);
    }

    CountData_destruct1DArray(countDataPerAnnotation, numberOfAnnotations);
    free(mostFrequentCoverages);
    stList_destruct(annotationNames);
    stList_destruct(annotationPaths);
    stHash_destruct(final_block_table);
    fprintf(stderr, "[%s] Done.\n", get_timestamp());
}

