#include "sam.h"
#include <time.h>
#include "bgzf.h"
#include <regex.h>
#include "sonLib.h"
#include <assert.h>
#include "faidx.h"
#include <getopt.h>
#include "cigar_it.h"
#include "thread_pool.h"

stHash* getSnpTable(char* vcfPath){
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    char* token;
    char contigName[50];
    char contigNamePrev[50];
    contigNamePrev[0] = '\0';
    int loc;
    stList* snpList;
    stHash* snpTable = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, (void (*)(void *)) stList_destruct);
    fp = fopen(vcfPath, "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);
    while ((read = getline(&line, &len, fp)) != -1) {
	if(line[0] == '#') continue;
        line[strlen(line)-1] = '\0';
	token = strtok(line, "\t");
        strcpy(contigName, token);
	// get location
        token = strtok(NULL, "\t");
        loc = atoi(token) - 1; // 0-based
	if (strcmp(contigNamePrev, contigName) != 0){ // Assume that vcf is sorted
		strcpy(contigNamePrev, contigName);
		snpList = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
		printf("%s\n", contigName);
		char* tmp = malloc(50);
		strcpy(tmp, contigName);
		stHash_insert(snpTable, tmp, snpList);
	}
	snpList = stHash_search(snpTable, contigName);
	stList_append(snpList, stIntTuple_construct1(loc));
    }
    return snpTable;
}


void filterReads(stHash* snpTable, char* inputPath, char* outputPath, char* filteredPath, int nthreads, int clusterMargin, double ratioThreshold, int minQual){
	samFile* fp = sam_open(inputPath, "r");
	samFile* fo = sam_open(outputPath, "wb");
	samFile* ff = sam_open(filteredPath, "wb");
	// Make a multi threading pool
        htsThreadPool p = {NULL, 0};
        if (nthreads > 0) {
                p.pool = hts_tpool_init(nthreads);
                if (!p.pool) {
                        fprintf(stderr, "Error creating thread pool\n");
                }
		else{ // Add I/O streams to the threading pool
			hts_set_opt(fp, HTS_OPT_THREAD_POOL, &p);
			hts_set_opt(fo, HTS_OPT_THREAD_POOL, &p);
			hts_set_opt(ff, HTS_OPT_THREAD_POOL, &p);
		}
        }
        sam_hdr_t* sam_hdr = sam_hdr_read(fp);
	sam_hdr_write(fo, sam_hdr);
	sam_hdr_write(ff, sam_hdr);
        bam1_t* b = bam_init1();
	ptCigarIt* cigarIt;
	int len;
	char* readName;
	int tid;
	char* contigName;
	char contigNamePrev[50];
	contigNamePrev[0] = '\0';
	stList* snpList = NULL;
	int startSnpIndex=0; int snpIndex=0;
	int loc;
	int locPrevious;
	stIntTuple* locTuple;
	bool writeFlag = true;
	int clusterMismatchCount = 0;
        int clusterTotalCount = 0;
	double clusterMismatchRatio = 0.0;
	uint8_t* qual;
	fprintf(stderr, "[LOG:] readname\tcontig\tstart\tcluster_mismatch_count\tcluster_total_count\tcluster_mismatch_ratio\n");
	while( sam_read1(fp, sam_hdr, b) > -1){
		if((b->core.flag & BAM_FUNMAP) > 0) continue;
		readName = bam_get_qname(b);
		qual = bam_get_qual(b);
                tid = b->core.tid;
                contigName = sam_hdr_tid2name(sam_hdr, tid);
		//update snps list if contig has changed
		if (strcmp(contigName, contigNamePrev) != 0){
			snpList = stHash_search(snpTable, contigName);
			startSnpIndex = 0;
			strcpy(contigNamePrev, contigName);
		}
		if ((b->core.flag & BAM_FSECONDARY) > 0) continue;
		if (snpList == NULL || stList_length(snpList) == 0) { // This contig has no snps
			if (sam_write1(fo, sam_hdr, b) == -1) {
                                fprintf(stderr, "Couldn't write %s\n", readName);
                        }
			continue;
		}
		locTuple = stList_get(snpList, startSnpIndex);
		loc = stIntTuple_get(locTuple, 0);
		while (loc < b->core.pos && startSnpIndex < stList_length(snpList) - 1){
			startSnpIndex++;
			locTuple = stList_get(snpList, startSnpIndex);
                	loc = stIntTuple_get(locTuple, 0);
			//printf("%d %d %d %d\n",loc, b->core.pos, startSnpIndex , stList_length(snpList));
		}
		//printf("%d\n",loc);
		if (loc < b->core.pos){
			//printf("writing\n");
			if (sam_write1(fo, sam_hdr, b) == -1) {
                                fprintf(stderr, "Couldn't write %s\n", readName);
                        }
			continue;
		}

		cigarIt = ptCigarIt_construct(b, true, true);
		snpIndex = startSnpIndex;
		writeFlag = true;
		clusterMismatchCount = 0;
		clusterTotalCount = 0;
		locPrevious = -1 * clusterMargin;
		clusterMismatchRatio = 0;
		if (b->core.pos <= loc) {
			clusterTotalCount = 1;
		}
		while(ptCigarIt_next(cigarIt)){
			if (cigarIt->op == BAM_CDIFF) {
				len = cigarIt->rfe - cigarIt->rfs + 1;
                                for(int i=0; i < len; i++){
					while(loc < (cigarIt->rfs + i) && snpIndex < stList_length(snpList) - 1){
						snpIndex++;
                        			locTuple = stList_get(snpList, snpIndex);
						locPrevious = loc;
                        			loc = stIntTuple_get(locTuple, 0);
						if (locPrevious + clusterMargin < loc && 0 < clusterTotalCount) { // one cluster is passed
							clusterMismatchRatio = (double) clusterMismatchCount / clusterTotalCount;
							if (ratioThreshold <= clusterMismatchRatio) {
								writeFlag = false;
								break;
							}
							clusterTotalCount = 0;
							clusterMismatchCount = 0;
						}
						clusterTotalCount += 1;
					}
					if (writeFlag == false) break;
					if (cigarIt->rfs + i == loc && minQual <= qual[cigarIt->sqs + i]){
						clusterMismatchCount += 1;
					}
					if (loc < (cigarIt->rfs + i)){
                       				break;
                			}
					//printf("##%d - %d\n", cigarIt->rfs, cigarIt->rfe);
					//printf("#%d\n", cigarIt->len);
				}
				if (writeFlag == false) break;
			}
		}

		if(0 < clusterTotalCount) clusterMismatchRatio = (double) clusterMismatchCount / clusterTotalCount;
		if(clusterMismatchRatio < ratioThreshold){
			if (sam_write1(fo, sam_hdr, b) == -1) {
                        	fprintf(stderr, "Couldn't write %s\n", readName);
			}
		}
		else {
			if (sam_write1(ff, sam_hdr, b) == -1) {
                                fprintf(stderr, "Couldn't write %s\n", readName);
                        }
			fprintf(stderr, "[LOG:] %s\t%s\t%d\t%d\t%d\t%.2f\n", readName, contigName, b->core.pos, clusterMismatchCount, clusterTotalCount, clusterMismatchRatio);
		}
		ptCigarIt_destruct(cigarIt);
	}
	printf("All Done!\n");
	sam_hdr_destroy(sam_hdr);
	bam_destroy1(b);
	sam_close(fo);
	sam_close(fp);
	sam_close(ff);
	if (p.pool)
                hts_tpool_destroy(p.pool);
}


int main(int argc, char *argv[]){
        int c;
        char* inputPath;
        char* outputPath;
	char* filteredPath;
	char* vcfPath;
	int nthreads = 2;
	int minQual = 10;
	int clusterMargin = 1000;
	double ratioThreshold = 0.3;
        char *program;
        (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
        while (~(c=getopt(argc, argv, "i:o:f:v:r:m:t:q:h"))) {
                switch (c) {
                        case 'i':
                                inputPath = optarg;
                                break;
                        case 'o':
                                outputPath = optarg;
                                break;
			case 'f':
				filteredPath = optarg;
				break;
		 	case 'v':
                                vcfPath = optarg;
                                break;
			case 't':
				nthreads = atoi(optarg);
				break;
			case 'm':
                                clusterMargin = atoi(optarg);
                                break;
			case 'r':
				ratioThreshold = atof(optarg);
				break;
			case 'q':
				minQual = atoi(optarg);
				break;
                        default:
                                if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
                        help:
                                fprintf(stderr, "\nUsage: %s  -i <INPUT_BAM> -v <VCF> -o <OUTPUT_BAM> \n\t Filter the reads that contain the alternative alleles of the snps in the given VCF\n\n", program);
                                fprintf(stderr, "Options:\n");
                                fprintf(stderr, "         -i         input bam file (CS tag is required)\n");
                                fprintf(stderr, "         -o         output bam file\n");
				fprintf(stderr, "         -f         output bam file that contains the removed reads with alternative alleles\n");
				fprintf(stderr, "         -v         vcf file containing biallelic snps (output of 'bcftools view -Ov -f PASS -m2 -M2 -v snps')\n");
				fprintf(stderr, "         -t         number of threads (for bam I/O) [Default: 2]\n");
				fprintf(stderr, "         -m         margin creating snp clusters [Default: 1000]\n");
				fprintf(stderr, "         -r         threshold for cluster mismatch ratio [Default: 0.3]\n");
				fprintf(stderr, "         -q         minimum base quality for a mismatch [Default: 10]\n");
                                return 1;
                }
        }
	stHash* snpTable = getSnpTable(vcfPath);
	fprintf(stderr, "snp table is created\n");
	/*
	stHashIterator* it = stHash_getIterator(snpTable);
	char* key;
	stList* list;
	stIntTuple* locTuple;
	while((key = stHash_getNext(it)) != NULL){
		printf("%s\n",key);
		list = stHash_search(snpTable, key);
		for (int i=0 ; i < stList_length(list); i++){
			locTuple = stList_get(list, i);
                        int loc = stIntTuple_get(locTuple, 0);
			printf("$ %d\n",loc);
		}
	}*/
	//printf("#\n");
	filterReads(snpTable, inputPath, outputPath, filteredPath, nthreads, clusterMargin, ratioThreshold, minQual);
	stHash_destruct(snpTable);
}
