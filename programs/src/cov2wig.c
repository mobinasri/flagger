#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include "sonLib.h"
#include "common.h"

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

void splitCov(char* covPath, int segmentLen, char* wigPath, char* name, int threshold){
    FILE* fp; FILE* fo;
    
    fp = fopen(covPath, "r");
    fo = fopen(wigPath, "w");
    
    if (fp == NULL)
        exit(EXIT_FAILURE);

    char* contig = malloc(200);
    int contigLength = 0;
    char* preContig = malloc(200);
    int start=0, end=0, cov=0;
    int segmentStart = 0; int segmentEnd = -1;
    int sumCov=0;
    float avgCov=0.0;
    while (readNextBlock(fp, &contig, &contigLength, &start, &end, &cov) == 1) {
		if ((strcmp(preContig, contig) != 0)) {
			if (segmentStart <= segmentEnd) {
				avgCov = (double) sumCov / (segmentEnd - segmentStart + 1);
                                avgCov = avgCov < threshold ? avgCov : threshold;
                                fprintf(fo,"%.2f\n", avgCov);
			}
			fprintf(fo, "track type=\"wiggle_0\" name=\"%s\"\nfixedStep chrom=%s start=%d step=%d span=%d\n", name, contig, start, segmentLen, segmentLen);
			segmentStart = start;
			segmentEnd = segmentStart - 1;
			sumCov = 0;
		}
		//iterate over segments untill the end of the currect coverage block is covered
		while ( segmentLen <= end - segmentStart + 1 )  {
			//printf("%d\t%d\n",blockLen, (end - blockStart + 1));
			segmentEnd = segmentStart + segmentLen - 1;
			sumCov += (segmentEnd - max(start, segmentStart) + 1) * cov;
			avgCov = (double) sumCov / (segmentEnd - segmentStart + 1);
                        avgCov = avgCov < threshold ? avgCov : threshold;
                        fprintf(fo,"%.2f\n", avgCov);
			segmentStart += segmentLen;
			segmentEnd = segmentStart - 1;
			sumCov = 0;
		}
		sumCov += (end - segmentEnd) * cov;
		segmentEnd = end;
		strcpy(preContig, contig);
    }
    fclose(fp);
    fflush(fo);
    fclose(fo);
}

int main(int argc, char *argv[]) {
   int c;
   int segmentSize=1024;
   char* trackName;
   char* faiPath;
   char* covPath;
   char* wigPath;
   int threshold=250;
   char *program;
   (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
   while (~(c=getopt(argc, argv, "i:o:t:f:s:n:h"))) {
		switch (c) {
			case 'i':
                                covPath = optarg;
                                break;
			case 'o':
                                wigPath = optarg;
                                break;
			case 'f':
                                faiPath = optarg;
                                break;
			case 's':
                                segmentSize = atoi(optarg);
                                break;
			case 'n':
                                trackName = optarg;
                                break;
			case 't':
				threshold = atoi(optarg);
				break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: %s  -f <FAI_FILE> -i <COVERAGE> -o <WIG> \n", program);
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -i         input coverage file\n");
				fprintf(stderr, "         -f         fai file\n");
				fprintf(stderr, "         -o         output path for wig file\n");
				fprintf(stderr, "         -s         segment size[default : 1024]\n");
				fprintf(stderr, "         -t         does not allow any coverage larger than this threshold[default : 250]\n");
				fprintf(stderr, "         -n         track name\n");
				return 1;	
		}		
   }
   splitCov(covPath, segmentSize, wigPath, trackName, threshold);
   return 0;
}
