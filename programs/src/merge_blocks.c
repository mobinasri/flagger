#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include "sonLib.h"
#include "common.h"
#include "block_it.h"


void mergeBlocks(Format_t format1, char* block1Path, 
		 Format_t format2, char* block2Path,
		 char* outputPath){
    FILE* fp1; FILE* fp2; FILE* fo;
    
    fp1 = fopen(block1Path, "r");
    fp2 = fopen(block2Path, "r");
    fo = fopen(outputPath, "w");
    
    if (fp1 == NULL){
	printf("Couldn't open %s", block1Path);
        exit(EXIT_FAILURE);
    }
    if (fp2 == NULL){
        printf("Couldn't open %s", block2Path);
        exit(EXIT_FAILURE);
    }
    if (fo == NULL){
        printf("Couldn't open %s", outputPath);
        exit(EXIT_FAILURE);
    }

    // Block for the first coverage file
    Block_t* block1 = Block_construct(format1);
    char preContig1[200];
    // Block for the second coverage file
    Block_t* block2 = Block_construct(format2);
    char preContig2[200];
    assert(Block_next(fp2, block2) == 1);
    while (Block_next(fp1, block1) == 1) {
	    if (strcmp(preContig1, block1->ctg) != 0){
		    fprintf(fo, ">%s %d\n", block1->ctg, block1->ctgLen);
            }
	    strcpy(preContig1, block1->ctg);
	    while(strcmp(block1->ctg, block2->ctg) == 0 &&
	          block2->e <= block1->e){
		    assert(block1->s <= block2->e);
                    fprintf(fo, "%d\t%d", max(block1->s, block2->s), block2->e);
		    // write block1 attrbs
                    for(int i=0; i < block1->attrbsLen; i++){
                            fprintf(fo, "\t%s", block1->attrbs[i]);
                    }
		    // write block2 attrbs
                    for(int i=0; i < block2->attrbsLen; i++){
                            fprintf(fo, "\t%s", block2->attrbs[i]);
                    }
                    fprintf(fo, "\n");
                    strcpy(preContig2, block2->ctg);
		    // if file is finished or contig has changed
                    if(Block_next(fp2, block2) != 1){
                            break;
                    }
            }
            if(strcmp(block1->ctg, block2->ctg) == 0 &&
	       block2->s <= block1->e && 
	       block1->e < block2->e){
                    fprintf(fo, "%d\t%d", max(block1->s, block2->s), block1->e);
		    // write block1 attrbs
                    for(int i=0; i < block1->attrbsLen; i++){
                            fprintf(fo, "\t%s", block1->attrbs[i]);
                    }
		    // write block2 attrbs
                    for(int i=0; i < block2->attrbsLen; i++){
                            fprintf(fo, "\t%s", block2->attrbs[i]);
                    }
                    fprintf(fo, "\n");
            }
    }
    free(block1);
    free(block2);
    fclose(fp1);
    fclose(fp2);
    fflush(fo);
    fclose(fo);
}

Format_t getFormat(char* path){
	// 3 is hardcoded here since both BED and COV have three letters
	char* prefix = path + strlen(path) - 3;
	if (strcmp(prefix, "bed") == 0){
		return BED;
	}
	else if (strcmp(prefix, "cov") == 0){
		return COV;
	}
	else{
		printf("format should be either cov or bed\n");
		exit(EXIT_FAILURE);
	}
}

int main(int argc, char *argv[]) {
   int c;
   char* path1;
   char* path2;
   char* outputPath;
   char *program;
   (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
   while (~(c=getopt(argc, argv, "a:b:o:h"))) {
		switch (c) {
			case 'a':
                                path1 = optarg;
                                break;
			case 'b':
                                path2 = optarg;
                                break;
			case 'o':
                                outputPath = optarg;
                                break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: %s  -a <COV/BED_FILE> -b <COV/BED_FILE> -o <OUT_COV_FILE> \n", program);
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -a         first input file (either .cov or .bed)\n");
				fprintf(stderr, "         -b         first input file (either .cov or .bed)\n");
				fprintf(stderr, "         -o         output path for merged cov file\n");
				return 1;	
		}		
   }
   mergeBlocks(getFormat(path1), path1, 
	       getFormat(path2), path2, 
	       outputPath);
   return 0;
}
