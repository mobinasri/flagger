#include<stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <inttypes.h>

void cov2counts(char* inPath, char* outPath){
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    char * token;
    int start=0, end=0, cov=0;
    int maxCov = 1000;
    int64_t* counts = malloc( (maxCov + 1) * sizeof(int64_t));
    // initialize counts
    for(int i = 0; i <= maxCov; i++){
	    counts[i] = 0;
    }
    fp = fopen(inPath, "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);
    while ((read = getline(&line, &len, fp)) != -1) {
        if (line[0] == '>') continue;
	line[strlen(line)-1] = '\0';
	//printf("%s\n", line);
	// get start location
	token = strtok(line, "\t");
	start = atoi(token);
	// get end location
        token = strtok(NULL, "\t");
	end = atoi(token);
	// get coverage
	token = strtok(NULL, "\t");
        cov = atoi(token);
	if (cov > maxCov) {
		counts = (int64_t*) realloc(counts, (cov + 1) * sizeof(int64_t));// since we expect a small number of coverages highar than 1000 it's seems efficient to increment the size one by one
		// initialize just added elements
		for(int i = (maxCov + 1); i <= cov; i++){
			counts[i] = 0;
		}
		maxCov = cov;
	}	
	counts[cov] += end - start + 1;
    }
    fclose(fp);
    // write to output
    fp = fopen(outPath, "w+");
    if (fp == NULL)
        exit(EXIT_FAILURE);
    for(int i = 0; i <= maxCov; i++){
            fprintf(fp, "%d\t%" PRId64 "\n", i, counts[i]);
    }
    fclose(fp);
}
int main(int argc, char *argv[]) {
   int c;
   char* inputPath;
   char* outputPath;
   char *program;
   (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
   while (~(c=getopt(argc, argv, "i:o:h"))) {
		switch (c) {
			case 'i':
				inputPath = optarg; 
				break;
			case 'o':
				outputPath = optarg; 
				break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: %s  -i -i <INPUT> -o <OUTPUT> \n", program);
				fprintf(stderr, "Options:\n");	
				fprintf(stderr, "         -i         input .cov file\n");	
				fprintf(stderr, "         -o         output .counts file\n");
				return 1;	
		}		
   }
   cov2counts(inputPath, outputPath);
   return 0;
}
