#include<stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>


int readContigSizes(char* faiPath, char*** contigNamesP, int** contigSizesP){
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    char * token;
    char* contigName;
    int contigSize;
    int numberOfContigs = 0;
    int arraySize = 100;
    char** contigNames = (char**) malloc( arraySize * sizeof(char*) );
    int* contigSizes = (int*) malloc( arraySize * sizeof(int) );
    fp = fopen(faiPath, "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);
    while ((read = getline(&line, &len, fp)) != -1) {
	line[strlen(line) - 1] = '\0';
        // get contig name
        token = strtok(line, "\t");
        contigName = token;
        // get contig size
        token = strtok(NULL, "\t");
        contigSize = atoi(token);
	numberOfContigs += 1;
        if (arraySize <= numberOfContigs) {
		arraySize *= 2;
                contigSizes = (int*) realloc(contigSizes, arraySize * sizeof(int));
		contigNames = (char**) realloc(contigNames, arraySize * sizeof(char*));
        }
	contigNames[numberOfContigs - 1] = (char*) malloc((strlen(contigName) + 1) * sizeof(char));
	strcpy(contigNames[numberOfContigs - 1], contigName);
        contigSizes[numberOfContigs - 1] = contigSize;
    }
    *contigNamesP = contigNames;
    *contigSizesP = contigSizes;
    fclose(fp);
    return numberOfContigs;
}

// Needs improvement maybe using a hashtable
int getContigSize(char* contigName, char** contigNames, int* contigSizes, int numberOfContigs){
	int idx = 0;
	while(strcmp(contigName,contigNames[idx]) != 0) idx++;
	if (idx >= numberOfContigs) exit(EXIT_FAILURE);
	return contigSizes[idx];
}

void depth2cov(char* depthPath, char** contigNames, int* contigSizes, int numberOfContigs, char* covPath){
    FILE * fout;
    FILE * fin;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    char * token;
    int pos=0, startPos=0, cov=0;
    int posNew=0, covNew=0;
    char contigName[50];
    char contigNameNew[50];
    fin = fopen(depthPath, "r");
    if (fin == NULL)
        exit(EXIT_FAILURE);
    fout = fopen(covPath, "w");
    if (fout == NULL)
        exit(EXIT_FAILURE);
    read = getline(&line, &len, fin);
    // Read the first position and its coverage
    line[strlen(line) - 1] = '\0';
    // get contig name
    token = strtok(line, "\t");
    strcpy(contigName, token);
    // get position
    token = strtok(NULL, "\t");
    startPos = atoi(token);
    // get coverage
    token = strtok(NULL, "\t");
    cov = atoi(token);
    pos = startPos;
    printf(">%s\n", contigName);
    fprintf(fout, ">%s %d\n", contigName, getContigSize(contigName, contigNames, contigSizes, numberOfContigs));
    //int p =0;
    while((read = getline(&line, &len, fin)) != -1) {
	line[strlen(line) - 1] = '\0';
        // get contig name
        token = strtok(line, "\t");
	strcpy(contigNameNew, token);
        // get position
        token = strtok(NULL, "\t");
        posNew = atoi(token);
        // get coverage
        token = strtok(NULL, "\t");
        covNew = atoi(token);
	// If a new contig is starting here, write the previous coverage
        if (strcmp(contigNameNew, contigName) != 0) {
            fprintf(fout, "%d\t%d\t%d\n", startPos, pos, cov);
            // Write the name of the new contig and its size
            fprintf(fout, ">%s %d\n", contigNameNew, getContigSize(contigNameNew, contigNames, contigSizes, numberOfContigs));
	    printf(">%s\n", contigNameNew);
            pos = posNew;
            startPos = pos;
	    strcpy(contigName, contigNameNew);
            cov = covNew;
	}
        // If contig and coverage is not changed update the position and continue reading
	else if (cov == covNew) {
             pos = posNew;
	}
	// If contig is the same and coverage changed write the previous coverage and corresponding block
        else {
             fprintf(fout, "%d\t%d\t%d\n", startPos, pos, cov);
             pos = posNew;
             startPos = pos;
             cov = covNew;
        }
    }
    // Write the coverage of the last block
    fprintf(fout, "%d\t%d\t%d\n", startPos, pos, cov);
    fclose(fin);
    fclose(fout);
}
int main(int argc, char *argv[]) {
   int c;
   char* faiPath;
   char* depthPath;
   char* covPath;
   char *program;
   (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
   while (~(c=getopt(argc, argv, "f:d:o:"))) {
		switch (c) {
			case 'f':
				faiPath = optarg; 
				break;
			case 'd':
                                depthPath = optarg;
                                break;
			case 'o':
				covPath = optarg; 
				break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: %s  -f <FAI> -d <DEPTH> -o <OUTPUT> \n", program);
				fprintf(stderr, "Options:\n");	
				fprintf(stderr, "         -f         input .fai file\n");
			        fprintf(stderr, "         -d         input .depth file\n");	
				fprintf(stderr, "         -o         output .cov file\n");
				return 1;	
		}		
   }
   char** contigNames;
   int* contigSizes;
   int numberOfContigs;
   printf("Reading fai file ...\n");
   numberOfContigs = readContigSizes(faiPath, &contigNames, &contigSizes);
   for(int i = 0; i < numberOfContigs; i++){
	   printf("%s\t%d\n",contigNames[i], contigSizes[i]);
   }
   printf("Coverting depth to cov ... \n");
   depth2cov(depthPath, contigNames, contigSizes, numberOfContigs, covPath);
   //release memory
   for(int i = 0; i < numberOfContigs; i++){
           free(contigNames[i]);
   }
   free(contigNames);
   free(contigSizes);
   return 0;
}
