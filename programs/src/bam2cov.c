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



int main(int argc, char *argv[]) {
    int c;
    int min_mapq = 20;
    double min_clipping_ratio = 0.1;
    int threads=4;
    char *bam_path;
    char *out_path;
    char *json_path;
    char *out_type;
    char *program;
    (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
    while (~(c = getopt(argc, argv, "i:t:j:m:r:O:o:h"))) {
        switch (c) {
            case 'i':
                bam_path = optarg;
                break;
            case 't':
                threads = atoi(optarg);
                break;
            case 'o':
                out_path = optarg;
                break;
            case 'j':
                json_path = optarg;
                break;
            case 'm':
                min_mapq = atoi(optarg);
                break;
            case 'r':
                min_clipping_ratio = atof(optarg);
                break;
            case 'O':
                out_type = optarg;
                break;
            default:
                if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
            help:
                fprintf(stderr, "\nUsage: %s  -i <BAM_FILE> -t <THREADS> -o <OUT_BED_FILE> \n", program);
                fprintf(stderr, "Options:\n");
                fprintf(stderr, "         -i         input bam file (should be indexed)\n");
                fprintf(stderr, "         -j         JSON file for the annotation bed files [maximum 32 files can be given and the keys can be any number between 1-32 for example {\"1\":\"/path/to/1.bed\", \"2\":\"/path/to/2.bed\"}]\n");
                fprintf(stderr, "         -m         minimum mapq for the measuring the coverage of the alignments with high mapq [Default = 20]\n");
                fprintf(stderr, "         -r         minimum clipping ratio for the measuring the coverage of the highly clipped alignments [Default = 0.1]\n");
                fprintf(stderr, "         -t         number of threads [default: 4]\n");
                fprintf(stderr, "         -O         output type [\"b\" for bed, \"bz\" for gzipped bed, \"c\" for cov, \"cz\" for gzipped cov]\n");
                fprintf(stderr, "         -o         output path \n");
                return 1;
        }
    }

    //merge and create the final block table
    stHash * final_block_table = ptBlock_multi_threaded_coverage_extraction_with_annotation(bam_path,
                                                                                            json_path,
                                                                                            threads,
                                                                                            min_mapq,
                                                                                            min_clipping_ratio);

    // get the table from contig name to contig length
    // it is only used once the output type is either .cov or .cov.gz
    stHash* ctg_to_len = ptBlock_get_contig_length_stHash_from_bam(bam_path);

    fprintf(stderr, "[%s] Started writing to %s.\n", get_timestamp(), out_path);
    if (strcmp(out_type, "c") == 0) { // uncompressed cov
        FILE* fp = fopen(out_path, "w");
        if (fp == NULL) {
            fprintf(stderr, "[%s] Error: Failed to open file %s.\n", get_timestamp(), out_path);
        }
        ptBlock_print_blocks_stHash_in_cov(final_block_table,
                                           get_string_cov_info_data_format_2,
                                           fp,
                                           false,
                                           ctg_to_len);
        fclose(fp);
    }else if(strcmp(out_type, "cz") == 0){ // gzip-compressed cov
        gzFile fp = gzopen(out_path, "w6h");
        if (fp == NULL) {
            fprintf(stderr, "[%s] Error: Failed to open file %s.\n", get_timestamp(), out_path);
        }
        ptBlock_print_blocks_stHash_in_cov(final_block_table,
                                           get_string_cov_info_data_format_2,
                                           &fp,
                                           true,
                                           ctg_to_len);
        gzclose(fp);
    }else if (strcmp(out_type, "b") == 0){ // uncompressed BED
        FILE* fp = fopen(out_path, "w");
        if (fp == NULL) {
            fprintf(stderr, "[%s] Error: Failed to open file %s.\n", get_timestamp(), out_path);
        }
        ptBlock_print_blocks_stHash_in_bed(final_block_table,
                                           get_string_cov_info_data_format_1,
                                           fp,
                                           false);
        fclose(fp);
    }else if (strcmp(out_type, "bz") == 0){ // gzip-compressed BED
        gzFile fp = gzopen(out_path, "w6h");
        if (fp == NULL) {
            fprintf(stderr, "[%s] Error: Failed to open file %s.\n", get_timestamp(), out_path);
        }
        ptBlock_print_blocks_stHash_in_bed(final_block_table,
                                           get_string_cov_info_data_format_1,
                                           &fp,
                                           true);
        gzclose(fp);
    }else{
        fprintf(stderr, "[%s] Error: Output type should only be from this list ['c','cz','b','bz']\n", get_timestamp());
        exit(EXIT_FAILURE);
    }
    stHash_destruct(final_block_table);
    stHash_destruct(ctg_to_len);
    fprintf(stderr, "[%s] Done.\n", get_timestamp());
}
