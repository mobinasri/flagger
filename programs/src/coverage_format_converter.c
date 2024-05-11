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

int main(int argc, char *argv[]) {
    int c;
    char *inputPath = NULL;
    char *outputPath = NULL;
    char *faiPath = NULL;
    int windowLen = 1000;
    int nThreads = 4;
    char *trackName = "coverage_wig";
    char *program;
    (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
    while (~(c = getopt(argc, argv, "i:o:f:w:t:n:h"))) {
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
	    case 'w':
		windowLen = atoi(optarg);
		break;
	    case 't':
		nThreads = atoi(optarg);
		break;
	    case 'n':
		trackName = optarg;
		break;
	    default:
                if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
            help:
                fprintf(stderr, "\nUsage: %s  -i <INPUT_FILE> -f <FAI> -o <OUTPUT_FILE> \n", program);
                fprintf(stderr, "Options:\n");
		fprintf(stderr, "         -i         input path (can have formats '.cov', '.cov.gz', '.bed' or '.bed.gz')\n");
		fprintf(stderr, "         -f         fai path\n");
                fprintf(stderr, "         -o         output path (can have formats '.cov', '.cov.gz', '.bed', '.bed.gz', or 'bedgraph')\n");
		fprintf(stderr, "         -w         window length (only be used for generating wig file) [Default = 1000]\n");
		fprintf(stderr, "         -t         number of threads for taking window average (only be used for generating wig file) [Default = 4]\n");
		fprintf(stderr, "         -n         track name (only be used for generating wig file) [Default = 'coverage_wig']\n");
                return 1;
        }
    }
    char *inputExtension = extractFileExtension(inputPath);
    char *outputExtension = extractFileExtension(outputPath);
    if (strcmp(outputExtension, "cov") != 0 &&
	strcmp(outputExtension, "cov.gz") != 0 &&
	strcmp(outputExtension, "bed") != 0 &&
	strcmp(outputExtension, "bed.gz") != 0 &&
	strcmp(outputExtension, "bedgraph") != 0){
	    fprintf(stderr, "[%s] Error: output file should either cov/cov.gz/bed/bed.gz/bedgraph  !\n", get_timestamp());
	    free(inputExtension);
	    free(outputExtension);
	    exit(EXIT_FAILURE);
    }
    if (strcmp(inputExtension, "cov") != 0 &&
        strcmp(inputExtension, "cov.gz") != 0 &&
        strcmp(inputExtension, "bed") != 0 &&
        strcmp(inputExtension, "bed.gz") != 0){
            fprintf(stderr, "[%s] Error: input file should either cov/cov.gz/bed/bed.gz  !\n", get_timestamp());
	    free(inputExtension);
            free(outputExtension);
            exit(EXIT_FAILURE);
    }

    fprintf(stderr, "[%s] Parsing %s.\n", get_timestamp(), faiPath);
    stHash *ctg_to_len = ptBlock_get_contig_length_stHash_from_fai(faiPath);

    char trackName2[100];
    char color[100];
    if (strcmp(outputExtension, "bedgraph") == 0){
	    char outputPath2[1000];
	    char prefix[1000];
	    memcpy(prefix, outputPath, strlen(outputPath)-9);
	    prefix[strlen(outputPath)-9] = '\0';

            int chunkCanonicalLen = ptBlock_get_max_contig_length(ctg_to_len);
            ChunksCreator *chunksCreator = ChunksCreator_constructFromCov(inputPath, faiPath, chunkCanonicalLen, nThreads, windowLen);
	    ChunksCreator_parseChunks(chunksCreator);
	    ChunksCreator_sortChunks(chunksCreator);

	    fprintf(stderr, "[%s] Writing %s.\n", get_timestamp(), outputPath);
	    // since chunkCanonicalLen is set to maximum contig size
	    // then each contig is parsed in a single chunk
	    sprintf(trackName2, "%s_total", trackName);
	    sprintf(outputPath2, "%s.bedgraph", prefix);
	    ChunksCreator_writeChunksIntoBedGraph(chunksCreator, outputPath2, trackName2, CoverageInfo_getCoverage, "w", "51,51,255");
	    
	    sprintf(trackName2, "%s_high_mapq", trackName);
	    sprintf(outputPath2, "%s.high_mapq.bedgraph", prefix);
            ChunksCreator_writeChunksIntoBedGraph(chunksCreator, outputPath2, trackName2, CoverageInfo_getCoverageHighMapq, "w", "0,0,153");
	    
	    sprintf(trackName2, "%s_high_clip", trackName);
	    sprintf(outputPath2, "%s.high_clip.bedgraph", prefix);
            ChunksCreator_writeChunksIntoBedGraph(chunksCreator, outputPath2, trackName2, CoverageInfo_getCoverageHighClip, "w", "204,0,204");

	    sprintf(trackName2, "%s_region_idx", trackName);
	    sprintf(outputPath2, "%s.region_idx.bedgraph", prefix);
            ChunksCreator_writeChunksIntoBedGraph(chunksCreator, outputPath2, trackName2, CoverageInfo_getRegionIndex, "w", "64,64,64");
	    
	    ChunksCreator_destruct(chunksCreator);
	    stHash_destruct(ctg_to_len);
	    free(outputExtension);
	    fprintf(stderr, "[%s] Done!\n", get_timestamp());
	    return 0;
    }

    fprintf(stderr, "[%s] Parsing %s.\n", get_timestamp(), inputPath);
    stHash *coverage_blocks_per_contig = ptBlock_parse_coverage_info_blocks(inputPath);

    fprintf(stderr, "[%s] Writing %s.\n", get_timestamp(), outputPath);
    if(strcmp(outputExtension, "cov") == 0){
	FILE* fp = fopen(outputPath, "w");
        if (fp == NULL) {
            fprintf(stderr, "[%s] Error: Failed to open file %s.\n", get_timestamp(), outputPath);
        }
	ptBlock_print_blocks_stHash_in_cov(coverage_blocks_per_contig,
                                           get_string_cov_info_data_format_2,
                                           fp,
                                           false,
                                           ctg_to_len);
	fclose(fp);
    }
    else if(strcmp(outputExtension, "cov.gz") == 0){
        gzFile fp = gzopen(outputPath, "w6h");
        if (fp == Z_NULL) {
            fprintf(stderr, "[%s] Error: Failed to open file %s.\n", get_timestamp(), outputPath);
        }
        ptBlock_print_blocks_stHash_in_cov(coverage_blocks_per_contig,
                                           get_string_cov_info_data_format_2,
                                           &fp,
                                           true,
                                           ctg_to_len);
	gzclose(fp);

    }
    else if(strcmp(outputExtension, "bed") == 0){
        FILE* fp = fopen(outputPath, "w");
        if (fp == NULL) {
            fprintf(stderr, "[%s] Error: Failed to open file %s.\n", get_timestamp(), outputPath);
        }
        ptBlock_print_blocks_stHash_in_bed(coverage_blocks_per_contig,
                                           get_string_cov_info_data_format_2,
                                           fp,
                                           false);
        fclose(fp);
    }
    else if(strcmp(outputExtension, "bed.gz") == 0){
        gzFile fp = gzopen(outputPath, "w6h");
        if (fp == Z_NULL) {
            fprintf(stderr, "[%s] Error: Failed to open file %s.\n", get_timestamp(), outputPath);
        }
        ptBlock_print_blocks_stHash_in_bed(coverage_blocks_per_contig,
                                           get_string_cov_info_data_format_2,
                                           &fp,
                                           true);
        gzclose(fp);
    }
    stHash_destruct(coverage_blocks_per_contig);
    stHash_destruct(ctg_to_len);
    free(outputExtension);
    fprintf(stderr, "[%s] Done!\n", get_timestamp());
}
