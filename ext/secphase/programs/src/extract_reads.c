#include <getopt.h>
#include "sam.h"
#include <assert.h>
#include "faidx.h"
#include <time.h>
#include "bgzf.h"
#include <regex.h>
#include "sonLib.h"

stHash* get_read_name_table(char* read_names_path){
	FILE* fp = fopen(read_names_path, "r");
	size_t read;
	size_t len;
	char* line = NULL;
	char read_name[100];
	char* key;
	char* token;
	stHash* read_name_table = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, NULL);
	while((read = getline(&line, &len, fp)) != -1){
		token = strtok(line, "\n");
        	strcpy(read_name, token);
		key = malloc(strlen(read_name) + 1);
		strcpy(key, read_name);
		stHash_insert(read_name_table, key, key);
	}
	return read_name_table;
}

int main(int argc, char *argv[]){
        int c;
        char* input_path;
	char* output_path;
	char* read_names_path;
        char *program;
        (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
        while (~(c=getopt(argc, argv, "i:o:r:h"))) {
                switch (c) {
                        case 'i':
                                input_path = optarg;
                                break;
                        case 'o':
                                output_path = optarg;
                                break;
			case 'r':
                                read_names_path = optarg;
                                break;
                        default:
                                if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
                        help:
                                fprintf(stderr, "\nUsage: %s  -i <INPUT_BAM> -o <OUTPUT_BAM> -r <READ_NAMES>\n", program);
                                fprintf(stderr, "Options:\n");
                                fprintf(stderr, "         -i         input bam file\n");
                                fprintf(stderr, "         -o         output dir\n");
				fprintf(stderr, "         -r         a single column file containing readnames\n");
                                return 1;
                }
        }
	//open input bam file and read the header
        samFile* fp = sam_open(input_path, "r");
        sam_hdr_t* sam_hdr = sam_hdr_read(fp);
	//open output bam file and write the header
	samFile* fo = sam_open(output_path, "w");
	sam_hdr_write(fo, sam_hdr);
	//get the table of read names
	stHash* read_name_table = get_read_name_table(read_names_path);
        char read_name[100];
        memset(read_name, '\0', 100);
	bam1_t* b = bam_init1();
	while(sam_read1(fp, sam_hdr, b) > -1) {
        	strcpy(read_name, bam_get_qname(b));
		if (stHash_search(read_name_table, read_name)){
			if (sam_write1(fo, sam_hdr, b) == -1) {
				fprintf(stderr, "Couldn't write %s\n", read_name);
                	}
		}
	}
	sam_close(fo);
	sam_close(fp);
        bam_destroy1(b);
        sam_hdr_destroy(sam_hdr);
}
