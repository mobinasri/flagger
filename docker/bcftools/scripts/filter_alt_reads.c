#include "sam.h"
#include <time.h>
#include "bgzf.h"
#include <regex.h>
#include "sonLib.h"
#include <assert.h>
#include "faidx.h"
#include <getopt.h>

#define CS_PATTERN "(:([0-9]+))|(([+-\\*])([a-z]+))"

//TODO: make a separate library for ptCigarIt
typedef struct {
	// the constant attributes
        uint32_t* cigar;
	int n; // number of operations
	bool is_rev;
	// the attributes below may change in each iteration
	uint8_t op; // current cigar operation
	int len; // length of the current operation
        int idx;
	// cs tag attributes
	bool use_cs;
	char* cs;
	int offset_cs;
	regex_t preg_cs;
	//location attributes
        int rds_f; // read start forward
	int rde_f; // read end forward
	int sqs; // seq start
	int sqe; // seq end
	int rfs; // ref start
	int rfe; // ref end
}ptCigarIt;


void ptCigarIt_destruct(ptCigarIt* cigar_it){
	regfree(&cigar_it->preg_cs);
	free(cigar_it);
}

ptCigarIt* ptCigarIt_construct(bam1_t* b, bool use_cs){
	ptCigarIt* cigar_it = (ptCigarIt*) malloc(sizeof(ptCigarIt));
	cigar_it->cigar = bam_get_cigar(b);
	cigar_it->len = 0;
	cigar_it->op = -1;
	cigar_it->is_rev = bam_is_rev(b);
	cigar_it->n = b->core.n_cigar;
	cigar_it->idx = -1;
	cigar_it->sqs = 0;
	cigar_it->sqe = -1;
	cigar_it->rfs = b->core.pos;
	cigar_it->rfe = b->core.pos - 1;
	int lclip = 0;
	int rclip = 0;
	if (bam_cigar_op(cigar_it->cigar[0]) == BAM_CHARD_CLIP) {
		lclip = bam_cigar_oplen(cigar_it->cigar[0]);
        }
        if (bam_cigar_op(cigar_it->cigar[cigar_it->n - 1]) == BAM_CHARD_CLIP) {
        	rclip = bam_cigar_oplen(cigar_it->cigar[cigar_it->n - 1]);
        }
	cigar_it->rds_f = cigar_it->is_rev ? rclip + lclip + b->core.l_qseq : 0;
	cigar_it->rde_f = cigar_it->is_rev ? rclip + lclip + b->core.l_qseq  - 1 : -1;
	cigar_it->use_cs = use_cs;
	int error;
	if(use_cs){
		cigar_it->offset_cs = 0;
		cigar_it->cs = bam_aux_get(b, "cs") + 1; // skip "Z" at the beginning
		if ((error = regcomp(&cigar_it->preg_cs, CS_PATTERN, REG_EXTENDED)) != 0) {
                	printf("regcomp() failed, returning nonzero (%d)\n", error);
                	exit(1);
        	}
		//printf("%s\n",cigar_it->cs);
	}
	return cigar_it;
}


int ptCigarIt_next_cs(ptCigarIt* cigar_it){
        regmatch_t pm;
	char* cs_shifted = cigar_it->cs + cigar_it->offset_cs;
	int error = regexec(&cigar_it->preg_cs, cs_shifted, 1, &pm, REG_NOTBOL);
	if (error != 0) return 0;
	int rd_step;
        int sq_step;
        int rf_step;
	char match_len[20];
	/*for(int i=pm.rm_so; i<pm.rm_eo; i++){
		printf("%c", cs_shifted[i]);
	}
	printf("\n");*/
        switch(cs_shifted[pm.rm_so]) {
		case ':':
			cigar_it->op = BAM_CEQUAL;
			memcpy(match_len, cs_shifted + 1, pm.rm_eo - pm.rm_so - 1);
			match_len[pm.rm_eo - pm.rm_so - 1] = '\0';
			cigar_it->len = atoi(match_len);
			rd_step = cigar_it->len;
                        sq_step = cigar_it->len;
                        rf_step = cigar_it->len;
                        break;
		case '*':
			cigar_it->op = BAM_CDIFF;
			cigar_it->len = 1;
                        rd_step = cigar_it->len;
                        sq_step = cigar_it->len;
                        rf_step = cigar_it->len;
                        break;
                case '+':
			cigar_it->op = BAM_CINS;
			cigar_it->len = pm.rm_eo - pm.rm_so - 1;
                        rd_step = cigar_it->len;
                        sq_step = cigar_it->len;
                        rf_step = 0;
                        break;
                case '-':
			cigar_it->op = BAM_CDEL;
			cigar_it->len = pm.rm_eo - pm.rm_so - 1;
                        rd_step = 0;
                        sq_step = 0;
                        rf_step = cigar_it->len;
                        break;
        }
        //shift read interval on forward strand
        if (cigar_it->is_rev) {
                cigar_it->rde_f = cigar_it->rds_f - 1;
                cigar_it->rds_f -= rd_step;
        }
        else {
                cigar_it->rds_f = cigar_it->rde_f + 1;
                cigar_it->rde_f += rd_step;
        }
        //shift seq interval
        cigar_it->sqs = cigar_it->sqe + 1;
        cigar_it->sqe += sq_step;
        //shift ref interval
        cigar_it->rfs = cigar_it->rfe + 1;
        cigar_it->rfe += rf_step;
	//update offset for cs tag
	cigar_it->offset_cs += pm.rm_eo;
	return cigar_it->len;
}

int ptCigarIt_next(ptCigarIt* cigar_it){
	if (cigar_it->idx == cigar_it->n - 1) return 0;
	uint32_t* cigar = cigar_it->cigar;
	if(cigar_it->use_cs){ // if it is set to use CS tag instead of cigar itself
		bool first_clip = cigar_it->idx == -1 && 
				  (bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP ||
                           	   bam_cigar_op(cigar[0]) == BAM_CHARD_CLIP);
		bool second_clip = cigar_it->idx == 0 &&
                                  (bam_cigar_op(cigar[1]) == BAM_CSOFT_CLIP ||
                                   bam_cigar_op(cigar[1]) == BAM_CHARD_CLIP);
		if(!(first_clip || second_clip)){//if no soft- or hard-clip at the beginning
			int cs_return = ptCigarIt_next_cs(cigar_it);
			// if CS is finished and there exist soft- or hard- clipping 
			// in one operation before the last one
			if(cs_return == 0 &&
			   (bam_cigar_op(cigar[cigar_it->n-2]) == BAM_CSOFT_CLIP ||
                           bam_cigar_op(cigar[cigar_it->n-2]) == BAM_CHARD_CLIP)){
				cigar_it->idx = cigar_it->idx = cigar_it->n-3;
				cigar_it->use_cs = false;
				return ptCigarIt_next(cigar_it);
                        }
			// if CS is finished and there exist soft- or hard- clipping
                        // in the last operation
			else if(cs_return == 0 &&
				(bam_cigar_op(cigar[cigar_it->n-1]) == BAM_CSOFT_CLIP ||
                           	bam_cigar_op(cigar[cigar_it->n-1]) == BAM_CHARD_CLIP )){
                                cigar_it->idx = cigar_it->idx = cigar_it->n-2;
                                cigar_it->use_cs = false;
                                return ptCigarIt_next(cigar_it);
                        }
			else return cs_return;
		}
	}
	cigar_it->idx += 1;
	int rd_step;
	int sq_step;
	int rf_step;
	cigar_it->op = bam_cigar_op(cigar[cigar_it->idx]);
	cigar_it->len = bam_cigar_oplen(cigar[cigar_it->idx]);
	switch(cigar_it->op) {
		case BAM_CMATCH:
		case BAM_CEQUAL:
		case BAM_CDIFF:
			rd_step = cigar_it->len;
			sq_step = cigar_it->len;
			rf_step = cigar_it->len;	
      			break;
		case BAM_CINS:
			rd_step = cigar_it->len;
                        sq_step = cigar_it->len;
                        rf_step = 0;
                        break;
   		case BAM_CDEL:
			rd_step = 0;
			sq_step = 0;
			rf_step = cigar_it->len;
			break;
		case BAM_CSOFT_CLIP:
			rd_step = cigar_it->len;
			sq_step = cigar_it->len;
			rf_step = 0;
			break;
		case BAM_CHARD_CLIP:
			rd_step = cigar_it->len;
                        sq_step = 0;
                        rf_step = 0;
      			break;
	}
	//shift read interval on forward strand
	if (cigar_it->is_rev) {
		cigar_it->rde_f = cigar_it->rds_f - 1;
		cigar_it->rds_f -= rd_step;
	}
	else {
		cigar_it->rds_f = cigar_it->rde_f + 1;
		cigar_it->rde_f += rd_step;
	}
	//shift seq interval
	cigar_it->sqs = cigar_it->sqe + 1;
	cigar_it->sqe += sq_step;
	//shift ref interval
	cigar_it->rfs = cigar_it->rfe + 1;
	cigar_it->rfe += rf_step;
	return cigar_it->len;
}


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


void filterReads(stHash* snpTable, char* inputPath, char* outputPath){
	samFile* fp = sam_open(inputPath, "r");
	samFile* fo = sam_open(outputPath, "wb");
        sam_hdr_t* sam_hdr = sam_hdr_read(fp);
	sam_hdr_write(fo, sam_hdr);
        bam1_t* b = bam_init1();
	ptCigarIt* cigarIt;
	int len;
	uint8_t* seq;
	char* readName;
	int tid;
	char* contigName;
	char contigNamePrev[50];
	contigNamePrev[0] = '\0';
	stList* snpList;
	int startSnpIndex=0; int snpIndex=0;
	int loc;
	stIntTuple* locTuple;
	bool writeFlag = true;
	while( sam_read1(fp, sam_hdr, b) > -1){
		if ((b->core.flag & BAM_FSECONDARY) > 0) continue;
		readName = bam_get_qname(b);
		seq = bam_get_seq(b);
		tid = b->core.tid;
                contigName = sam_hdr_tid2name(sam_hdr, tid);
		if (strcmp(contigName, contigNamePrev) != 0){
			snpList = stHash_search(snpTable, contigName);
			startSnpIndex = 0;
		}
		if (snpList == NULL) { // This contig has no snps
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
		//printf("Start Iterating!\n");
		cigarIt = ptCigarIt_construct(b, true);
		snpIndex = startSnpIndex;
		writeFlag = true;
		while(ptCigarIt_next(cigarIt)){
			if (cigarIt->op == BAM_CDIFF) {
				len = cigarIt->rfe - cigarIt->rfs + 1;
                                for(int i=0; i < len; i++){
					while(loc < (cigarIt->rfs + i) & snpIndex < stList_length(snpList) - 1){
						snpIndex++;
                        			locTuple = stList_get(snpList, snpIndex);
                        			loc = stIntTuple_get(locTuple, 0);
					}
					if (cigarIt->rfs + i == loc){
						writeFlag = false;
						break;
					}
					if (loc < (cigarIt->rfs + i)){
                       				break;
                			}
					//printf("##%d - %d\n", cigarIt->rfs, cigarIt->rfe);
					//printf("#%d\n", cigarIt->len);
				}
			}
		}
		if(writeFlag){
			if (sam_write1(fo, sam_hdr, b) == -1) {
                        	fprintf(stderr, "Couldn't write %s\n", readName);
			}
		}
	}
	sam_close(fo);
	sam_close(fp);
}


int main(int argc, char *argv[]){
        int c;
        char* inputPath;
        char* outputPath;
	char* vcfPath;
        char *program;
        (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
        while (~(c=getopt(argc, argv, "i:o:v:h"))) {
                switch (c) {
                        case 'i':
                                inputPath = optarg;
                                break;
                        case 'o':
                                outputPath = optarg;
                                break;
		 	case 'v':
                                vcfPath = optarg;
                                break;
                        default:
                                if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
                        help:
                                fprintf(stderr, "\nUsage: %s  -i <INPUT_BAM> -v <VCF> -o <OUTPUT_BAM> \n", program);
                                fprintf(stderr, "Options:\n");
                                fprintf(stderr, "         -i         input bam file (CS tag is required)\n");
                                fprintf(stderr, "         -o         output bam file\n");
				fprintf(stderr, "         -v         vcf file containing biallelic snps\n");
                                return 1;
                }
        }
	stHash* snpTable = getSnpTable(vcfPath);
	fprintf(stderr, "snp table is created\n");
	stHashIterator* it = stHash_getIterator(snpTable);
	char* key;
	/**while((key = stHash_getNext(it)) != NULL){
		printf("%s\n",key);
	}**/
	filterReads(snpTable, inputPath, outputPath);
	stHash_destruct(snpTable);
}
