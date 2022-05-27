#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>

int getBlockTypeIndex(float* probArray){
    float max = -1;
    int index = -1;
    for(int i = 2; i < 6; i++){
        if(max < probArray[i]){
	    index = i;
	    max = probArray[i];
	}
    }
    return index - 2;
}

int readNextBlock(FILE* fileReader, char** contig, int*contigLength, int* blockStart, int* blockEnd, int* coverage){
    size_t len = 0;
    char* line = NULL;
    char* token;
    int start, end;
    ssize_t read = getline(&line, &len, fileReader);
    if (read == -1){
	 printf("## FILE END\n");
         *contig = NULL;
         *blockStart = -1;
	 *blockEnd = -1;
	 *coverage = -1;
	 return 0;
    }
    if (line[0] == '>') {
	token = strtok(line, " ");
        strcpy(*contig, token + 1); // skip '>' and copy
	printf("%s\n",*contig);
	token = strtok(NULL, " ");
	*contigLength = atoi(token);
        readNextBlock(fileReader, contig, contigLength, blockStart, blockEnd, coverage);
    }
    else {
	token = strtok(line, "\t");
        *blockStart = atoi(token);
	token = strtok(NULL, "\t");
	*blockEnd = atoi(token);
	token = strtok(NULL, "\t");
        *coverage = atoi(token);
    }
    free(line);
    return 1;
}

void findBlocks(char* covPath, char* prefix){
    FILE* fp; FILE* fo = NULL;
    char outputPath[200];
    
    fp = fopen(covPath, "r");
    
    if (fp == NULL)
        exit(EXIT_FAILURE);

    char* contig = malloc(50);
    int contigLength = 0;
    char* preContig = malloc(50);
    int blockStart=0;
    int start=0, end=0, cov=0;

    while (readNextBlock(fp, &contig, &contigLength, &start, &end, &cov) == 1) {
	if ((strcmp(preContig, contig) != 0)){
		if (fo != NULL){
			fflush(fo);
			fclose(fo);
		}
		strcpy(outputPath, prefix);
    		strcat(outputPath, ".");
		strcat(outputPath, contig);
		strcat(outputPath, ".cov");
		fo = fopen(outputPath, "w");
		fprintf(fo, ">%s %d\n", contig, contigLength);
		strcpy(preContig, contig);
	}
	fprintf(fo, "%d\t%d\t%d\n", start, end, cov);
    }
    fflush(fo);
    fclose(fo);
    fclose(fp);
}
int main(int argc, char *argv[]) {
   int c;
   char* covPath;
   char* prefix;
   char *program;
   (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
   while (~(c=getopt(argc, argv, "t:c:p:h"))) {
		switch (c) {
			case 'c':
                                covPath = optarg;
                                break;
			case 'p':
                                prefix = optarg;
                                break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: %s  -t <TABLE> -c <COVERAGE> - <PREFIX> \n", program);
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -c         coverage file\n");
				fprintf(stderr, "         -p         prefix for the output cov files\n");
				return 1;	
		}		
   }
   findBlocks(covPath, prefix);
   return 0;
}
