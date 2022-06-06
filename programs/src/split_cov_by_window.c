#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include "sonLib.h"
#include "common.h"


typedef struct{
	char* contigName;
	int contigLength;
	int start; //inclusive 1-based
	int end; //inclusive 1-based
}Window;


Window* Window_construct(char* contigName, int length, int start, int end){
	Window* wnd = malloc(sizeof(Window));
	wnd->contigName = malloc(50);
	strcpy(wnd->contigName, contigName);
	wnd->contigLength = length;
	wnd->start = start;
	wnd->end = end;
	return wnd;
}

void Window_destruct(void* window){
	Window* wnd = window;
        free(wnd->contigName);
	wnd->contigName = NULL;
	free(wnd);
}



stList* getWindows(char* faiPath, int windowSize){
	stList* windows = stList_construct3(0, Window_destruct);
	char contigName[50]; int s, e, contigSize;
	FILE* faif = fopen(faiPath, "r");
	size_t len = 0;
    	char* line = NULL;
    	char* token;
	int n;
	while(getline(&line, &len, faif) != -1) {
    		token = strtok(line, "\t");
        	strcpy(contigName, token);
        	token = strtok(NULL, "\t");
        	contigSize = atoi(token);
		n = contigSize / windowSize;
		for(int i=0; i < n; i++){
			s = i * windowSize + 1;
			e = ( s + 2 * windowSize - 1 <= contigSize) ? (i + 1) * windowSize : contigSize;
			stList_append(windows, Window_construct(contigName, contigSize, s, e)); //inclusive 1-based
		}
		if (n == 0){
			s = 1;
			e = contigSize;
			stList_append(windows, Window_construct(contigName, contigSize, s, e)); //inclusive 1-based
		}
	}
	return windows;
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

void splitCov(char* covPath, char* prefix, stList* windows){
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
    int wnd_idx = 0;
    Window* wnd = stList_get(windows, wnd_idx);
    printf("%s\t%d\t%d\n", wnd->contigName, wnd->start, wnd->end);
    sprintf(outputPath, "%s.%s_%d_%d.cov", prefix, wnd->contigName, wnd->start, wnd->end);
    fo = fopen(outputPath, "w");
    while (readNextBlock(fp, &contig, &contigLength, &start, &end, &cov) == 1) {
	    if ((strcmp(preContig, contig) != 0)) fprintf(fo, ">%s %d\n", contig, contigLength);
	    fprintf(fo, "%d\t%d\t%d\n", max(start, wnd->start), min(wnd->end, end), cov);
	    while (wnd && strcmp(wnd->contigName, contig) == 0 && wnd->end <= end) {
			    wnd_idx++;
			    if (wnd_idx < stList_length(windows)) {
			    	wnd = stList_get(windows, wnd_idx);
			    }
			    else {
				    wnd = NULL;
				    break;
			    }
			    printf("%s\t%d\t%d\n", wnd->contigName, wnd->start, wnd->end);
			    fflush(fo);
			    fclose(fo);
			    sprintf(outputPath, "%s.%s_%d_%d.cov", prefix, wnd->contigName, wnd->start, wnd->end);
			    fo = fopen(outputPath, "w");
			    if (strcmp(wnd->contigName, contig) == 0) fprintf(fo, ">%s %d\n", contig, contigLength);
			    if (strcmp(wnd->contigName, contig) == 0 && wnd->start <= end){
				    fprintf(fo, "%d\t%d\t%d\n", max(start, wnd->start), min(wnd->end, end), cov);
			    }
	    }
	    strcpy(preContig, contig);
    }
    fflush(fo);
    fclose(fo);
    fclose(fp);
}

int main(int argc, char *argv[]) {
   int c;
   int windowSize=5e6;
   char* faiPath;
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
                                faiPath = optarg;
                                break;
			case 's':
                                windowSize = atoi(optarg);
                                break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: %s  -f <FAI_FILE> -c <COVERAGE> -p <PREFIX> -s <SEGMENT_SIZE> \n", program);
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -c         coverage file\n");
				fprintf(stderr, "         -f         fai file\n");
				fprintf(stderr, "         -p         prefix for the output cov files\n");
				fprintf(stderr, "         -s         window size[default : 5Mb]\n");
				return 1;	
		}		
   }
   stList* windows = getWindows(faiPath, windowSize);
   splitCov(covPath, prefix, windows);
   stList_destruct(windows);
   return 0;
}
