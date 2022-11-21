#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include "sonLib.h"
#include "common.h"

stHash* getEffectiveContigLens(char* contigLenPath){
	stHash* contigLenTable = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, free);
	char* contigName; int contigSize;
	FILE* f = fopen(contigLenPath, "r");
        size_t len = 0;
        char* line = NULL;
        char* token;
        int n;
        while(getline(&line, &len, f) != -1) {
                token = strtok(line, "\t");
		contigName = malloc(50* sizeof(char));
                strcpy(contigName, token);
                token = strtok(NULL, "\t");
                int* contigSize = malloc(sizeof(int));
		*contigSize = atoi(token);
		fprintf(stderr,"%s,%d\n", contigName, *contigSize);
		stHash_insert(contigLenTable, contigName, contigSize);
	}
	return contigLenTable;
}


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
         *contig = NULL;
         *blockStart = -1;
	 *blockEnd = -1;
	 *coverage = -1;
	 return 0;
    }
    if (line[0] == '>') {
	token = strtok(line, " ");
        strcpy(*contig, token + 1); // skip '>' and copy
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

void splitCov(char* covPath, char* prefix, stHash* contigLenTable, int windowLen){
    FILE* fp; FILE* fo = NULL;
    char outputPath[200];
    
    fp = fopen(covPath, "r");
    
    if (fp == NULL)
        exit(EXIT_FAILURE);

    char* contig = malloc(50); contig[0] = '\0';
    int contigLength = 0;
    char* preContig = malloc(50); preContig[0] = '\0';
    int blockStart=0;
    int start=0, end=0, cov=0;
    int remainingWindowLen = 0;
    int remainingContigLen = 0;
    int windowIdx = 0;
    int nWindows;
    int fileStart;
    int* effContigLenPtr;
    int effContigLen;
    int sum=0;
    fprintf(stderr, "iterating over cov blocks\n");
    while (readNextBlock(fp, &contig, &contigLength, &start, &end, &cov) == 1) {
	    if ((strcmp(preContig, contig) != 0) || contig[0] == '\0'){
		    fprintf(stderr, "contig changed\n");
		    // open the first cov file for the current contig
		    // tmp.cov will be renamed to ${ctg}_${start}_${end}.cov once finished and closed
		    fo = fopen("tmp.cov", "w");
		    // write the header of the current contig
		    fprintf(fo, ">%s %d\n", contig, contigLength);
		    // get the effective total length of the current contig
		    effContigLenPtr = stHash_search(contigLenTable, contig);
		    fprintf(stderr, "contig len fetched\n");
		    if(effContigLenPtr == NULL) fprintf(stderr, "%s, NULL\n",contig);
                    effContigLen = *effContigLenPtr;
		    remainingContigLen = effContigLen;
		    nWindows = (int) ((double) effContigLen / windowLen);
		    nWindows = nWindows == 0 ? 1 : nWindows;
		    fprintf(stderr, "contig len =%d, nwindows = %d\n", effContigLen, nWindows);
		    windowIdx = 0;
		    remainingWindowLen = (nWindows == 1) ? remainingContigLen : windowLen;
		    fileStart = start;
		    sum =0;
	    }
	    int blockLen = end - start + 1;

	    if (blockLen < remainingWindowLen){
		    remainingWindowLen -= blockLen;
		    remainingContigLen -= blockLen;
		    sum += blockLen;
		    //fprintf(stderr, "%d\n", remainingContigLen);
		    fprintf(fo, "%d\t%d\t%d\n", start, end, cov);
	    }
	    // Note that the length of the last block of contig should be equal to remainingWindowLen
	    else if(remainingWindowLen <= blockLen){
		    // save last block in the currect cov file and close the currect cov file
		    fprintf(fo, "%d\t%d\t%d\n", start, start + remainingWindowLen - 1, cov);
		    fflush(fo);
                    fclose(fo);
		    // rename the saved cov file since we now know the end location
		    sprintf(outputPath, "%s.%s_%d_%d.cov", prefix, contig, fileStart, start + remainingWindowLen - 1);
		    printf("%s\t%d\t%d\n", contig, fileStart, start + remainingWindowLen - 1);
		    fprintf(stderr, "window finished, %s_%d_%d\n", contig, fileStart, start + remainingWindowLen - 1);
		    fprintf(stderr, "sum=%d,remainingContigLen = %d, remainingWindowLen= %d\n", sum, remainingContigLen, remainingWindowLen);
		    rename("tmp.cov", outputPath);
		    // if we have more blocks for the current contig, open a new tmp cov file
		    if (remainingWindowLen < remainingContigLen){
			    fo = fopen("tmp.cov", "w");
			    fileStart = start + remainingWindowLen;
			    fprintf(fo, ">%s %d\n", contig, contigLength);
			    fprintf(fo, "%d\t%d\t%d\n", start + remainingWindowLen, end, cov);
			    windowIdx += 1;
			    if (windowIdx < nWindows - 1){
				    remainingWindowLen = windowLen;
			    }
			    else{ // windowIdx == nWindows - 1 (last window of the current contig)
				    remainingWindowLen = remainingContigLen - blockLen;
			    }
			    remainingContigLen -= blockLen;
		    }
	    }
	    strcpy(preContig, contig);
    }
    // close the given coverage file
    fclose(fp);
}

int main(int argc, char *argv[]) {
   int c;
   int windowLen=5e6;
   char* contigLenPath;
   char* covPath;
   char* prefix;
   char *program;
   (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
   while (~(c=getopt(argc, argv, "c:p:f:s:h"))) {
		switch (c) {
			case 'c':
                                covPath = optarg;
                                break;
			case 'p':
                                prefix = optarg;
                                break;
			case 'f':
                                contigLenPath = optarg;
                                break;
			case 's':
                                windowLen = atoi(optarg);
                                break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: %s  -f <CTG_LEN_FILE> -c <COVERAGE> -p <PREFIX> -s <WINDOW_LENGTH> \n", program);
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -c         coverage file\n");
				fprintf(stderr, "         -f         a 2-column tab-delimited file (1st: contig name, 2st: effective contig length)\n");
				fprintf(stderr, "         -p         prefix for the output cov files\n");
				fprintf(stderr, "         -s         window size[default : 5Mb]\n");
				return 1;	
		}		
   }
   fprintf(stderr, "reading contig lens\n");
   stHash* contigLenTable = getEffectiveContigLens(contigLenPath);
   fprintf(stderr, "splitting cov\n");
   splitCov(covPath, prefix, contigLenTable, windowLen);
   fprintf(stderr, "desctructing contig table\n");
   stHash_destruct(contigLenTable);
   return 0;
}
