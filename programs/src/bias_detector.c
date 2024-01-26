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
    while ((contig_name = stHash_getNext(it)) != NULL) {
        blocks = stHash_search(final_block_table, contig_name);
        for (int i = 0; i < stList_length(blocks); i++) {
            ptBlock *block = stList_get(blocks, i);
            CoverageInfo *covInfo = (CoverageInfo *) block->data;
            int annotationIndex = getFirstIndexWithNonZeroBitFromRight(covInfo->annotation_flag);
            int count = block->rfe - block->rfs + 1;
            int value = covInfo->coverage;
            CountData_increment(countDataPerAnnotation[annotationIndex], value, (double) count);
        }
    }
    return countDataPerAnnotation;
}

int* getMostFrequentCoveragesPerAnnotation(CountData** countDataPerAnnotation, int numberOfAnnotations){
    int *mostFrequentCoverages = malloc(numberOfAnnotations * sizeof(int));
    for (int annotationIndex = 0; annotationIndex < numberOfAnnotations; annotationIndex++) {
        int mostFreqCoverage = CountData_getMostFrequentValue(countDataPerAnnotation[annotationIndex]);
        mostFrequentCoverages[annotationIndex] = mostFreqCoverage;
    }
    return mostFrequentCoverages;
}

int main(int argc, char *argv[]) {
    int c;
    int min_mapq = 20;
    double min_clipping_ratio = 0.1;
    int threads=4;
    char *bam_path;
    char *json_path;
    char *baseline_annot_name;
    char *program;
    (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
    while (~(c = getopt(argc, argv, "i:t:j:h"))) {
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
            default:
                if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
            help:
                fprintf(stderr, "\nUsage: %s  -i <BAM_FILE> -t <THREADS> -j <ANNOTATION_JSON> \n", program);
                fprintf(stderr, "Options:\n");
                fprintf(stderr, "         -i         input bam file (should be indexed)\n");
                fprintf(stderr, "         -j         JSON file for the annotation bed files [maximum 32 files can be given and the keys can any string {\"hsat1\":\"/path/to/hsat1.bed\", \"bsat\":\"/path/to/bsat.bed\"}]\n");
                fprintf(stderr, "         -t         number of threads [default: 4]\n");
                fprintf(stderr, "         -b         name of the baseline annotation\n");
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

    stList *annotation_names = parse_annotation_names_and_save_in_stList(json_path);
    int numberOfAnnotations = stList_length(annotation_names);
    int baselineIndex = get_annotation_index(annotation_names, baseline_annot_name);

    CountData** countDataPerAnnotation = createCountDataPerAnnotation(final_block_table, numberOfAnnotations);
    int* mostFrequentCoverages = getMostFrequentCoveragesPerAnnotation(countDataPerAnnotation, numberOfAnnotations);
    int baselineCoverage = mostFrequentCoverages[baselineIndex];

    for(int annotIndex=0; annotIndex < numberOfAnnotations; annotIndex++){
        char* annotation = stList_get(annotation_names, annotIndex);
        double coverageDiffRatio = abs((double) mostFrequentCoverages[annotIndex]  - baselineCoverage) / baselineCoverage;
        fprintf(stdout, "[%s]\tmod_cov=%d\tcov_diff_ratio = %.3f\n",
                annotation,
                mostFrequentCoverages[annotIndex],
                coverageDiffRatio);
    }

    CountData_destruct1DArray(countDataPerAnnotation, numberOfAnnotations);
    free(mostFrequentCoverages);
    stList_destruct(annotation_names);
    stHash_destruct(final_block_table);
    fprintf(stderr, "[%s] Done.\n", get_timestamp());
}

