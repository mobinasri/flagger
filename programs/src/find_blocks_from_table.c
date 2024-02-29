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

int readNextBlock(FILE* fileReader, char** contig, int* blockStart, int* blockEnd, int* coverage){
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
        readNextBlock(fileReader, contig, blockStart, blockEnd, coverage);
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

float** readProbTable(char* tablePath, int* maxCovObserved){
    FILE* fp;
    char* line = NULL;
    size_t len = 0;
    ssize_t read;
    char* token;
    int cov=0;
    float hapProb=0.0;
    int maxCov = 2000;
    float** table = malloc( (maxCov + 1) * sizeof(float*));
    for(int i = 0; i < (maxCov + 1); i++){
	    table[i] = malloc( 6 * sizeof(float));
    }
    fp = fopen(tablePath, "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);
    while ((read = getline(&line, &len, fp)) != -1) {
        if (line[0] == '#') continue;
	token = strtok(line, "\t");
        cov = atoi(token);
        if (cov > maxCov) {
            table = (float**) realloc(table, (2 * cov + 1) * sizeof(float*));
            // initialize just added elements
            for(int i = (maxCov + 1); i <= 2 * cov; i++){
                table[i] = malloc( 6 * sizeof(float));
		for(int j = 0; j < 6; j++){
			table[i][j] = 0;
		}
            }
            maxCov = 2 * cov;
        }
	if (cov > *maxCovObserved) {
		*maxCovObserved = cov;
	}
	token = strtok(NULL, "\t");
        table[cov][0] = strtof(token, NULL);
        token = strtok(NULL, "\t");
        table[cov][1] = strtof(token, NULL);
        token = strtok(NULL, "\t");
        table[cov][2] = strtof(token, NULL);
        token = strtok(NULL, "\t");
        table[cov][3] = strtof(token, NULL);
        token = strtok(NULL, "\t");
        table[cov][4] = strtof(token, NULL);
        token = strtok(NULL, "\t");
        table[cov][5] = strtof(token, NULL);
    }
    return table;
    //for(int i = 0; i <= maxCov; i++){
    //    printf("%d\t%f\n", i, table[i]);
    //}

}

void findBlocks(float** probTable, char* covPath, char* prefix, int maxCovObserved){
    FILE* fp; FILE* fo;
    char suffixes[4][20] = {".error.bed", ".duplicated.bed", ".haploid.bed", ".collapsed.bed"};
    char outputPath[1000];
    FILE** foArray = malloc( 4 * sizeof(FILE*));
    for(int i = 0; i < 4; i++){
	    strcpy(outputPath, prefix);
	    strcat(outputPath, suffixes[i]);
            foArray[i] = fopen(outputPath, "w");
	    if (foArray[i] == NULL)
		    exit(EXIT_FAILURE);
    }
    int* cov2index = malloc( (maxCovObserved + 1) * sizeof(int));
    printf("%d\n",maxCovObserved);
    for (int i = 0; i <= maxCovObserved; i++){
	    cov2index[i] = getBlockTypeIndex(probTable[i]);
	    printf("%d -> %d\n", i, cov2index[i]);
    }
    
    fp = fopen(covPath, "r");
    
    if (fp == NULL)
        exit(EXIT_FAILURE);

    char* contig = malloc(100);
    char* preContig = malloc(100);
    int blockStart=0, preEnd=0, preCov=0;
    int start=0, end=0, cov=0;

    if(readNextBlock(fp, &preContig, &blockStart, &preEnd, &preCov) != 1){
	    for(int i = 0; i < 4; i++){
                fclose(foArray[i]);
            }
	    fclose(fp);
	    return;
    }
    strcpy(contig, preContig);
    while (readNextBlock(fp, &contig, &start, &end, &cov) == 1) {
	if ((strcmp(preContig, contig) == 0) && ((preEnd + 1) == start) && (cov2index[preCov] == cov2index[cov])){
		preEnd = end;
		preCov = cov;
		continue;
	}
	else {
		fprintf(foArray[cov2index[preCov]], "%s\t%d\t%d\n", preContig, blockStart - 1, preEnd);
		strcpy(preContig, contig);
		blockStart = start;
		preEnd = end;
		preCov = cov;
	}
    }
    fprintf(foArray[cov2index[preCov]], "%s\t%d\t%d\n", preContig, blockStart - 1, preEnd);
    for(int i = 0; i < 4; i++){
	    fflush(foArray[i]);
	    fclose(foArray[i]);
    }
    fclose(fp);
}
int main(int argc, char *argv[]) {
   int c;
   char* tablePath;
   char* covPath;
   char* prefix;
   char *program;
   (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
   while (~(c=getopt(argc, argv, "t:c:p:h"))) {
		switch (c) {
			case 't':
				tablePath = optarg; 
				break;
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
				fprintf(stderr, "         -t         probability table\n");
				fprintf(stderr, "         -c         coverage file\n");
				fprintf(stderr, "         -p         prefix for the output bed file\n");
				return 1;	
		}		
   }
   printf("Reading the prob table\n");
   int maxCovObserved = 0;
   float** probTable = readProbTable(tablePath, &maxCovObserved);
   //for(int i = 0; i <= 100; i++){
   //    printf("%d\t%f\n", i, hapTable1[i]);
   //}
   findBlocks(probTable, covPath, prefix, maxCovObserved);
   free(probTable);
   return 0;
}
