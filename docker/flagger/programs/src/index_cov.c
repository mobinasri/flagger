#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include "sonLib.h"
#include "common.h"
#include "hmm.h" 


int main(int argc, char *argv[]) {
   int c;
   char covPath[200];
   int chunkLen;
   char *program;
   (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
   while (~(c=getopt(argc, argv, "c:l:h"))) {
		switch (c) {
			case 'c':
                                strcpy(covPath, optarg);
                                break;
			case 'l':
                                chunkLen = atoi(optarg);
                                break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
			help:	
				fprintf(stderr, "\nUsage: %s  -c <COV_FILE> -l <CHUNK_LEN> \n", program);
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -c         path to .cov file\n");
				fprintf(stderr, "         -l         chunk length\n");
				return 1;	
		}		
   }
   stList* chunks = createCovIndex(covPath, chunkLen);
   char covIndexPath[200];
   sprintf(covIndexPath, "%s.index", covPath);
   writeCovIndex(chunks, covIndexPath);
   stList* parsedChunks = parseCovIndex(covIndexPath);
   for(int i = 0; i < stList_length(parsedChunks); i++){
	   Chunk* chunk = stList_get(parsedChunks, i);
	   printf("%d:\t%s\t%d\t%d\t%ld\n", i, chunk->ctg, chunk->s, chunk->e, chunk->fileOffset);
   }
   stList_destruct(chunks);
   stList_destruct(parsedChunks);
}
