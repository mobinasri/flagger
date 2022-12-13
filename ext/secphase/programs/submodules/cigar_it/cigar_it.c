#include "sam.h"
#include "cigar_it.h"
#include "common.h"
#include <regex.h>
#include <getopt.h>


void ptCigarIt_destruct(ptCigarIt* cigar_it){
	if(cigar_it->use_cs) regfree(&cigar_it->preg_cs);
	if(cigar_it->use_md) regfree(&cigar_it->preg_md);
	free(cigar_it);
}

ptCigarIt* ptCigarIt_construct(bam1_t* b, bool use_cs, bool use_md){
	if((use_cs || use_md) == false) {
		fprintf(stderr, "At least one of the CS and MD modes should be used for cigar iteration");
		exit(1);
	}
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
	cigar_it->match_remain=0;
	cigar_it->use_cs = false;
	cigar_it->use_md = false;
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
	int error;
	cigar_it->cs = bam_aux_get(b, "cs");
	cigar_it->md = bam_aux_get(b, "MD");
	if(cigar_it->cs && use_cs){
		cigar_it->cs += 1; // skip "Z" at the beginning
		cigar_it->use_cs = true;
		cigar_it->offset_cs = 0;
		if ((error = regcomp(&cigar_it->preg_cs, CS_PATTERN, REG_EXTENDED)) != 0) {
                        printf("regcomp() failed, returning nonzero (%d)\n", error);
                        exit(1);
                }
	}
	else if (cigar_it->md && use_md){
		cigar_it->md += 1; // skip "Z" at the beginning
		cigar_it->use_md = true;
                cigar_it->offset_md = 0;
                if ((error = regcomp(&cigar_it->preg_md, MD_PATTERN, REG_EXTENDED)) != 0) {
                        printf("regcomp() failed, returning nonzero (%d)\n", error);
                        exit(1);
                }
	}
	else{
		fprintf(stderr, "At least one of the MD or CS tags should be present!\n");
		exit(1);
	}
	return cigar_it;
}


int ptCigarIt_next_md(ptCigarIt* cigar_it){
        regmatch_t pm;
        char* md_shifted = cigar_it->md + cigar_it->offset_md;
        int error = regexec(&cigar_it->preg_md, md_shifted, 1, &pm, REG_NOTBOL);
        if (error != 0) return 0;
        /*int rd_step;
        int sq_step;
        int rf_step;*/
        char match_len[20];
        /*for(int i=pm.rm_so; i<pm.rm_eo; i++){
                printf("%c", md_shifted[i]);
        }
        printf("\n");*/
	if(md_shifted[pm.rm_so] == '0'){ // if it starts with 0 it means two consecutive mismacthes 
                        cigar_it->op = BAM_CDIFF;
                        //printf("equal\n");
                        cigar_it->len = 0;
                        /*rd_step = cigar_it->len;
                        sq_step = cigar_it->len;
                        rf_step = cigar_it->len;*/
        }
	else if(md_shifted[pm.rm_so] <= '9'){ // if it starts with a number
                        cigar_it->op = BAM_CEQUAL;
			//printf("equal\n");
                        memcpy(match_len, md_shifted, pm.rm_eo - pm.rm_so);
			match_len[pm.rm_eo - pm.rm_so] = '\0';
                        cigar_it->len = atoi(match_len);
                        /*rd_step = cigar_it->len;
                        sq_step = cigar_it->len;
                        rf_step = cigar_it->len;*/
        }
	else if(md_shifted[pm.rm_so] < 90) { // if it is an alphabet letter
                        cigar_it->op = BAM_CDIFF;
                        cigar_it->len = 1 + (pm.rm_eo - pm.rm_so - 1) / 2;
                        /*rd_step = cigar_it->len;
                        sq_step = cigar_it->len;
                        rf_step = cigar_it->len;*/
	}
        else if(md_shifted[pm.rm_so] == 94) { // 94 is the ascii code of '^'
                        cigar_it->op = BAM_CDEL;
                        cigar_it->len = pm.rm_eo - pm.rm_so - 1;
                        /*rd_step = 0;
                        sq_step = 0;
                        rf_step = cigar_it->len;*/
        }
	//printf("%d\t%d\n", cigar_it->op, cigar_it->len);
        /*
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
        */
        //update offset for cs tag
        cigar_it->offset_md += pm.rm_eo;
	// between consecutive mismatches in MD tag there placed an additional zero
	if (cigar_it->len == 0) ptCigarIt_next_md(cigar_it);
	//printf("## %d\n",cigar_it->len);
        return cigar_it->len;
}



int ptCigarIt_next_cs(ptCigarIt* cigar_it){
        regmatch_t pm;
	char* cs_shifted = cigar_it->cs + cigar_it->offset_cs;
	int error = regexec(&cigar_it->preg_cs, cs_shifted, 1, &pm, REG_NOTBOL);
	if (error != 0) return 0;
	/*int rd_step;
        int sq_step;
        int rf_step;*/
	char match_len[20];
	/*printf("\n\t#");
	for(int i=pm.rm_so; i<pm.rm_eo; i++){
		printf("%c", cs_shifted[i]);
	}
	printf("\n");*/
        switch(cs_shifted[pm.rm_so]) {
		case ':':
			cigar_it->op = BAM_CEQUAL;
			memcpy(match_len, cs_shifted + 1, pm.rm_eo - pm.rm_so - 1);
			match_len[pm.rm_eo - pm.rm_so - 1] = '\0';
			cigar_it->len = atoi(match_len);
			/*rd_step = cigar_it->len;
                        sq_step = cigar_it->len;
                        rf_step = cigar_it->len;*/
                        break;
		case '*':
			cigar_it->op = BAM_CDIFF;
			cigar_it->len = (pm.rm_eo - pm.rm_so + 1) / 3;
                        /*rd_step = cigar_it->len;
                        sq_step = cigar_it->len;
                        rf_step = cigar_it->len;*/
                        break;
                case '+':
			cigar_it->op = BAM_CINS;
			cigar_it->len = pm.rm_eo - pm.rm_so - 1;
                        /*rd_step = cigar_it->len;
                        sq_step = cigar_it->len;
                        rf_step = 0;*/
                        break;
		case '-':
			cigar_it->op = BAM_CDEL;
			cigar_it->len = pm.rm_eo - pm.rm_so - 1;
                        /*rd_step = 0;
                        sq_step = 0;
                        rf_step = cigar_it->len;*/
                        break;
        }
	/*
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
	*/
	//update offset for cs tag
	cigar_it->offset_cs += pm.rm_eo;
	return cigar_it->len;
}

int ptCigarIt_next(ptCigarIt* cigar_it){
	if (cigar_it->idx == cigar_it->n - 1) return 0;
	uint32_t* cigar = cigar_it->cigar;
	cigar_it->idx += 1;
	int rd_step;
	int sq_step;
	int rf_step;
	uint8_t op = bam_cigar_op(cigar[cigar_it->idx]);
	int len = bam_cigar_oplen(cigar[cigar_it->idx]);
	int md_len; // an intermediate variable for MD tag
	//printf("###%d\t%d\n", op, len);
	
	switch(op) {
		case BAM_CMATCH:
		case BAM_CEQUAL:
                case BAM_CDIFF:
			if(cigar_it->match_remain == 0) cigar_it->match_remain = len;
			if(cigar_it->use_cs) {
				ptCigarIt_next_cs(cigar_it);
				cigar_it->match_remain -= cigar_it->len;
				if (0 < cigar_it->match_remain){
                                        cigar_it->idx -= 1; // return back to "M" until it is fully iterated by MD or CS
                                }
			}
			else if(cigar_it->use_md){
				if(0 <= cigar_it->match_remain) ptCigarIt_next_md(cigar_it);
				// In the following lines
				// we handle the case when
				// Cigar ops and MD ops does not end 
				// in the same location
				if(cigar_it->match_remain < 0){ // Cigar ops should be proceeded
					cigar_it->op = BAM_CEQUAL;
					cigar_it->len = min(len, -1 * cigar_it->match_remain);
					cigar_it->match_remain += len;
				}
				else { // MD ops should be proceeded
					md_len = cigar_it->len;
                                	cigar_it->len = min(cigar_it->len, cigar_it->match_remain);
					cigar_it->match_remain -= md_len;
				}
				if (0 < cigar_it->match_remain){
					cigar_it->idx -= 1; // return back to "M" until it is fully iterated by MD or CS
				}
			}
			rd_step = cigar_it->len;
			sq_step = cigar_it->len;
			rf_step = cigar_it->len;
      			break;
		case BAM_CINS:
			if(cigar_it->use_cs) ptCigarIt_next_cs(cigar_it);
			//MD does not have "I"
			cigar_it->len = len;
			cigar_it->op = op;
			rd_step = len;
                        sq_step = len;
                        rf_step = 0;
                        break;
   		case BAM_CDEL:
			if(cigar_it->use_cs) ptCigarIt_next_cs(cigar_it);
			else if(cigar_it->use_md) ptCigarIt_next_md(cigar_it);
			rd_step = 0;
			sq_step = 0;
			rf_step = len;
			break;
		case BAM_CSOFT_CLIP:
			cigar_it->len = len;
                        cigar_it->op = op;
			rd_step = len;
			sq_step = len;
			rf_step = 0;
			break;
		case BAM_CHARD_CLIP:
			cigar_it->len = len;
                        cigar_it->op = op;
			rd_step = len;
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
