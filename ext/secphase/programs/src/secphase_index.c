#include <getopt.h>
#include "sam.h"
#include "bgzf.h"
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
#include "edlib.h"
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include "ptBlock.h"
#include "ptVariant.h"
#include "ptAlignment.h"
#include "ptMarker.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

static struct option long_options[] =
        {
                {"inputBam",          required_argument, NULL, 'i'},
                {"stepSize",          required_argument, NULL, 's'},
                {NULL,                0,                 NULL, 0}
        };


int main(int argc, char *argv[]) {
    int c;
    int step_size = 10;
    char inputPath[200];
    char *program;
    (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
    while (~(c = getopt_long(argc, argv, "i:h", long_options, NULL))) {
        switch (c) {
            case 'i':
                strcpy(inputPath, optarg);
                break;
            case 's':
                step_size = atoi(optarg);
                break;
            default:
                if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
            help:
                fprintf(stderr, "\nUsage: %s  -i <INPUT_BAM> \n", program);
                fprintf(stderr, "Options:\n");
                fprintf(stderr, "         --inputBam, -i         Input BAM file\n");
                fprintf(stderr, "         --stepSize, -s         Step size for indexing\n");
                return 1;
        }
    }
    int64_t step_log_index = 1;
    int64_t step_log_size = 100000;

    // open input sam/bam file for parsing alignment records
    samFile *fp = sam_open(inputPath, "r");
    sam_hdr_t *sam_hdr = sam_hdr_read(fp);
    bam1_t *b = bam_init1();
    char read_name[100];
    char read_name_new[100];
    memset(read_name, '\0', 100);
    memset(read_name_new, '\0', 100);
    int64_t count_parsed_reads = 0;
    int64_t step_index = 1;
    int64_t addresses_size = 100000;
    int64_t addresses_number = 0;
    int64_t* addresses = malloc(addresses_size * sizeof(int64_t));
    // save the address of the first alignment
    addresses[addresses_number] = bgzf_tell(fp->fp.bgzf);
    addresses_number += 1;
    int bytes_read;
    while (true) {
        bytes_read = sam_read1(fp, sam_hdr, b);
        if (bytes_read > -1) {
            strcpy(read_name_new, bam_get_qname(b));
            if (read_name[0] == '\0') {
                strcpy(read_name, read_name_new);
            }
        }
	    if (bytes_read <= -1) break; // file is finished so break
        // If read name has changed or file is finished
        if ((strcmp(read_name_new, read_name) != 0) || (bytes_read <= -1)) {
            count_parsed_reads += 1;
        }
        if (count_parsed_reads == step_index * step_size) {
            if (addresses_size == addresses_number) addresses_size *= 2;
            addresses = realloc(addresses, addresses_size * sizeof(int64_t));
            addresses[addresses_number] = bgzf_tell(fp->fp.bgzf);
            addresses_number += 1;
            step_index += 1;
        }
        if (count_parsed_reads == step_log_index * step_log_size) {
            fprintf(stderr, "[%s] # Parsed reads = %d\n", get_timestamp(), count_parsed_reads);
            step_log_index += 1;
        }
    }
    // make sure to save the last address in the file
    if (count_parsed_reads >=  (step_index - 1) * step_size) {
        if (addresses_size == addresses_number) addresses_size += 1;
        addresses = realloc(addresses, addresses_size * sizeof(int64_t));
        addresses[addresses_number] = bgzf_tell(fp->fp.bgzf);
        addresses_number += 1;
    }
    fprintf(stderr, "[%s] # Total parsed reads = %d\n", get_timestamp(), count_parsed_reads);
    char index_path[200];
    snprintf(index_path, 200, "%s.secphase.index", inputPath);
    fprintf(stderr, "[%s] Writing index file %s\n", get_timestamp(), index_path);
    FILE *index_fp = fopen(index_path,"wb");
    // the first 8 bytes is the number of addresses in this index file
    fwrite(&addresses_number, sizeof(int64_t), 1, index_fp);
    // the remaining bytes are addresses one after another (each address occupies 8 bytes)
    fwrite(addresses, sizeof(int64_t), addresses_number, index_fp);
    sam_hdr_destroy(sam_hdr);
    sam_close(fp);
    bam_destroy1(b);
    fclose(index_fp);
    fprintf(stderr, "[%s] Done!\n", get_timestamp());
}
