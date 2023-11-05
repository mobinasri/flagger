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


stList* parse_all_annotations_and_save_in_stList(char* json_path){
    int buffer_size = 0;
    char* json_buffer = read_whole_file(json_path, &buffer_size, "r");
    fwrite(json_buffer,1,buffer_size,stderr);
    cJSON *annotation_json = cJSON_ParseWithLength(json_buffer, buffer_size);
    if (annotation_json == NULL)
    {
        const char *error_ptr = cJSON_GetErrorPtr();
        if (error_ptr != NULL)
        {
            fprintf(stderr, "Error before: %s\n", error_ptr);
        }
        return NULL;
    }

    int annotation_count = cJSON_GetArraySize(annotation_json);
    stList* block_table_list = stList_construct3(annotation_count, stHash_destruct);

    // iterate over key-values in json
    // each key is an index
    // each value is a path to a bed file
    cJSON *element = NULL;
    cJSON_ArrayForEach(element, annotation_json)
    {
	if (cJSON_IsString(element)){
            char* bed_path = cJSON_GetStringValue(element);
            stHash *annotation_block_table = ptBlock_parse_bed(bed_path);
	    fprintf(stderr, "[%s] Parsed  annotation %s:%s\n", get_timestamp(), element->string, cJSON_GetStringValue(element));
            int index = atoi(element->string) - 1; // given index is 1-based
            stList_set(block_table_list, index, annotation_block_table);
        }
    }
    cJSON_Delete(annotation_json);
    fprintf(stderr, "[%s] Number of parsed annotations = %d\n", get_timestamp(), stList_length(block_table_list));
    //ptBlock_print_blocks_stHash_in_bed(stList_get(block_table_list, 0), false, stderr, false);
    return block_table_list;

}

// annotation block tables need to have a CoverageInfo object as their "data"
// this CoverageInfo will have zero values for the coverage related attributes
// the annotation flag is set based on the order of the annotation blocks in the 
// given stList as the input
// The data augmentation is happening in place
void add_coverage_info_to_all_annotation_block_tables(stList *block_table_list){
    for(int i=0; i < stList_length(block_table_list); i++){
        stHash * block_per_contig = stList_get(block_table_list, i);
        // Each annotation has an associated flag represented by a bit-vector
        // the size of the bit-vector is 32, so it can be saved in an int32_t variable
        // for example for i=0 -> flag = 1 and for i=6 -> flag= 64
        int32_t annotation_flag = 1 << i;
        CoverageInfo * cov_info = CoverageInfo_construct(annotation_flag, 0, 0, 0);
        // The cov_info object created above will be copied and added as "data" to all blocks
        // for the current annotation. This process is happening in place
        ptBlock_add_data_to_all_blocks_stHash(block_per_contig,
                                              cov_info,
                                              destruct_cov_info_data,
                                              copy_cov_info_data,
                                              extend_cov_info_data);
    }
}

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

    // parse annotation bed files
    stList* annotation_block_table_list = parse_all_annotations_and_save_in_stList(json_path);

    // add coverage info objects as data to all annotation blocks
    // each coverage info will contain only the related annotation flag with 0 coverage
    add_coverage_info_to_all_annotation_block_tables(annotation_block_table_list);

    // parse alignments and make a block table that contains the necessary coverage values per block
    // the coverage values will be related to the total alignments, alignments with high mapq and
    // each block in the output table is a maximal contiguous block with no change in the depth of coverage
    stHash* coverage_block_table = ptBlock_multi_threaded_coverage_extraction(bam_path,
                                                                              threads,
                                                                              min_mapq,
                                                                              min_clipping_ratio);
    // print len/number stats for the coverage block table
    fprintf(stderr, "[%s] Created block table with coverage data : tot_len=%ld, number=%ld\n", get_timestamp(),
		    ptBlock_get_total_length_by_rf(coverage_block_table), 
		    ptBlock_get_total_number(coverage_block_table));

    // cover the whole genome with blocks that have zero coverage
    // this is useful to save the blocks with no coverage
    stHash* whole_genome_block_table = ptBlock_get_whole_genome_blocks_per_contig(bam_path);
    // print len/number stats for the whole genome block table
    fprintf(stderr, "[%s] Created block table for whole genome  : tot_len=%ld, number=%ld\n", get_timestamp(),
		    ptBlock_get_total_length_by_rf(whole_genome_block_table), 
		    ptBlock_get_total_number(whole_genome_block_table));

    CoverageInfo * cov_info = CoverageInfo_construct(0, 0, 0, 0);
    // The cov_info object created above will be copied and added as "data" to all whole genome blocks
    // This process is happening in place.
    ptBlock_add_data_to_all_blocks_stHash(whole_genome_block_table,
                                          cov_info,
                                          destruct_cov_info_data,
                                          copy_cov_info_data,
                                          extend_cov_info_data);

    // print len/number stats for the whole genome block table
    fprintf(stderr, "[%s] Created block table for whole genome  : tot_len=%ld, number=%ld\n", get_timestamp(),
		    ptBlock_get_total_length_by_rf(whole_genome_block_table), 
		    ptBlock_get_total_number(whole_genome_block_table));


    // add annotation and coverage blocks to the whole genome blocks in place
    ptBlock_extend_block_tables(whole_genome_block_table, coverage_block_table);
    fprintf(stderr, "[%s] Added 0-coverage whole genome blocks to coverage block tables : tot_len=%ld, number=%ld\n", get_timestamp(), 
		    ptBlock_get_total_length_by_rf(whole_genome_block_table), 
		    ptBlock_get_total_number(whole_genome_block_table));
    for(int i=0; i < stList_length(annotation_block_table_list); i++){
        ptBlock_extend_block_tables(whole_genome_block_table, stList_get(annotation_block_table_list, i));
    }
    fprintf(stderr, "[%s] Added annotation blocks to coverage block tables: \n", get_timestamp(), 
		    ptBlock_get_total_length_by_rf(whole_genome_block_table),
		    ptBlock_get_total_number(whole_genome_block_table));

    fprintf(stderr, "[%s] Started sorting and merging blocks\n", get_timestamp());
    //sort
    ptBlock_sort_stHash_by_rfs(whole_genome_block_table);
    fprintf(stderr, "[%s] Merged blocks : tot_len=%ld, number=%ld\n", get_timestamp(),
		    ptBlock_get_total_length_by_rf(whole_genome_block_table), 
		    ptBlock_get_total_number(whole_genome_block_table));

    //merge and create the final block table
    stHash * final_block_table = ptBlock_merge_blocks_per_contig_by_rf_v2(whole_genome_block_table);

    fprintf(stderr, "[%s] Created final block table : tot_len=%ld, number=%ld\n", get_timestamp(),
		    ptBlock_get_total_length_by_rf(final_block_table),
		    ptBlock_get_total_number(final_block_table));


    // free unmerged blocks
    stHash_destruct(coverage_block_table);
    stHash_destruct(whole_genome_block_table);
    stList_destruct(annotation_block_table_list);

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
