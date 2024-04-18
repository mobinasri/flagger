#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include "sonLib.h"
#include "common.h"
#include "chunk.h" 


int main(int argc, char *argv[]) {
   int c;
   char *covPath;
   char *faiPath = NULL;
   int chunkCanonicalLen = 20000000;
   char *program;
   (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
   while (~(c=getopt(argc, argv, "i:l:f:h"))) {
		switch (c) {
			case 'i':
                                covPath = optarg;
                                break;
			case 'f':
				faiPath = optarg;
				break;
			case 'l':
                                chunkCanonicalLen = atoi(optarg);
                                break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
			help:	
				fprintf(stderr, "\nUsage: %s  -i <COV_OR_BED_FILE> -l <CHUNK_LEN> [-f <FAI>] \n", program);
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -i         path to input coverage file (can be either cov/cov.gz/bed/bed.gz)\n");
				fprintf(stderr, "         -f         path to fai file (It can be skipped if the coverage file is cov/cov.gz)\n");
				fprintf(stderr, "         -l         chunk length [Default = 20000000 (20Mb)]\n");
				return 1;	
		}		
   }
   stList* templateChunks = ChunksCreator_createCovIndex(covPath, faiPath, chunkCanonicalLen);
   char covIndexPath[1000];
   sprintf(covIndexPath, "%s.index", covPath);
   ChunksCreator_writeCovIndex(templateChunks, covIndexPath);
   for(int i = 0; i < stList_length(templateChunks); i++){
	   Chunk* chunk = stList_get(templateChunks, i);
	   printf("[%s]Chunk %d:\t%s\t%d\t%d\t%ld\n", get_timestamp(), i, chunk->ctg, chunk->s, chunk->e, chunk->fileOffset);
   }
   stList_destruct(templateChunks);
}
