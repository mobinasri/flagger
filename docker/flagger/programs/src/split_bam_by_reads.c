#include <getopt.h>
#include "sam.h"
//#include "bgzf.h"
#include <assert.h>

int main(int argc, char *argv[]){
        int c;
        char* input_path;
	char* output_dir;
	int n;
        char *program;
        (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
        while (~(c=getopt(argc, argv, "i:o:n:h"))) {
                switch (c) {
                        case 'i':
                                input_path = optarg;
                                break;
                        case 'o':
                                output_dir = optarg;
                                break;
			case 'n':
                                n = atoi(optarg);
                                break;
                        default:
                                if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
                        help:
                                fprintf(stderr, "\nUsage: %s  -i <INPUT_BAM> -o <OUTPUT_DIR> -n <N>\n", program);
                                fprintf(stderr, "Options:\n");
                                fprintf(stderr, "         -i         input bam file\n");
                                fprintf(stderr, "         -o         output dir\n");
				fprintf(stderr, "         -n         number of reads per bam\n");
                                return 1;
                }
        }
        samFile* fp = sam_open(input_path, "r");
        sam_hdr_t* sam_hdr = sam_hdr_read(fp);
        char read_name[100];
        char read_name_new[100];
        memset(read_name, '\0', 100);
        memset(read_name_new, '\0', 100);
	int reads_count = 0;
	int split_idx = 0;
	char output_path[100];
        sprintf(output_path, "%s/%d.bam", output_dir, split_idx);
        samFile* fo = sam_open(output_path, "wb");
	sam_hdr_write(fo, sam_hdr);
	bam1_t* b = bam_init1();
	while(sam_read1(fp, sam_hdr, b) > -1) {
        	strcpy(read_name_new, bam_get_qname(b));
		if (strcmp(read_name_new, read_name) != 0) {
			strcpy(read_name, read_name_new);
			reads_count++;
			if (reads_count > n){
				sam_close(fo);
				split_idx++;
                                sprintf(output_path, "%s/%d.bam", output_dir, split_idx);
                                fo = sam_open(output_path, "wb");
				sam_hdr_write(fo, sam_hdr);
				reads_count = 1;
			}
		}
		if (sam_write1(fo, sam_hdr, b) == -1) {
			fprintf(stderr, "Couldn't write %s\n", read_name);
                }

	}
	sam_close(fo);
	sam_close(fp);
        bam_destroy1(b);
        sam_hdr_destroy(sam_hdr);
}
