#include "sam.h"
#include <assert.h>
#include <getopt.h>
#include "bgzf.h"


int main(int argc, char *argv[]){
        int c;
        char* inputPath;
	char* outputPath;
        int threshold = 30;
        char *program;
        (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
        while (~(c=getopt(argc, argv, "i:t:o:h"))) {
                switch (c) {
                        case 'i':
                                inputPath = optarg;
                                break;
                        case 'o':
                                outputPath = optarg;
                                break;
                        case 't':
                                threshold = atoi(optarg);
                                break;
                        default:
                                if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
                        help:
                                fprintf(stderr, "\nUsage: %s  -i <INPUT_BAM> -o <OUTPUT_BAM> -t <THRESHOLD> \n", program);
                                fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -i         input bam file\n");
                                fprintf(stderr, "         -t         mapq threshold\n");
                                fprintf(stderr, "         -o         output bam file\n");
				return 1;
		}
	}
	samFile* fp = sam_open(inputPath, "r");
        samFile* fo = sam_open(outputPath, "wb");
        sam_hdr_t* sam_hdr = sam_hdr_read(fp);
        sam_hdr_write(fo, sam_hdr);
        bam1_t* b = bam_init1();
	int qual;
	while(sam_read1(fp, sam_hdr, b) > -1){
                if ((b->core.flag & BAM_FSECONDARY) > 0) continue; // skip secondary alignments
                if (b->core.qual < threshold) b->core.qual = threshold; // modify mapq
                if (sam_write1(fo, sam_hdr, b) == -1) {
                	fprintf(stderr, "Couldn't write %s\n", bam_get_qname(b));
		}
	}
	sam_hdr_destroy(sam_hdr);
	sam_close(fo);
	sam_close(fp);
}

