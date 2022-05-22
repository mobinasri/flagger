#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "sonLib.h"
#include "common.h"
#include "block_it.h"


Block_t* Block_construct(Format_t format){
	Block_t* b = malloc(sizeof(Block_t));
	b->format = format;
	b->ctgLen = -1;
	b->s=-1; b->e=-1;
	b->attrbs = NULL;
	b->attrbsLen = 0;
}

void Block_destruct(Block_t* block){
	// free block attrbs
    	for(int i = 0; i < block->attrbsLen; i++){
                 free(block->attrbs[i]);
	}
    	free(block->attrbs);
    	block->attrbs = NULL;
	free(block);
}

int Block_next(FILE* fileReader, Block_t* block){
	if (block->format == COV){
		return readNextBlockCov(fileReader, block);
	}
	else if(block->format == BED){
		return readNextBlockBed(fileReader, block);
	}
	else{
		printf("ERROR: FORMAT should be either BED or COV!\n");
		exit(EXIT_FAILURE);
	}
}

int readNextBlockBed(FILE* fileReader, Block_t* block){
    size_t len = 0;
    char* line = NULL;
    char* token;
    ssize_t read = getline(&line, &len, fileReader);
    line[read-1] = '\0';
    // free block attrbs to fill it with new ones
    for(int i = 0; i < block->attrbsLen; i++){
                 free(block->attrbs[i]);
    }
    free(block->attrbs);
    block->attrbs = NULL;
    block->attrbsLen = 0;
    block->ctgLen = -1;
    // if this is the end of the file
    if (read == -1){
         block->ctg[0] = '\0';
         block->s = -1;
         block->e = -1;
         return 0;
    }
    Splitter* splitter = Splitter_construct(line, '\t');
    token = Splitter_getToken(splitter);
    strcpy(block->ctg, token);
    token = Splitter_getToken(splitter);
    block->s = atoi(token) + 1; // start make it 1-based 
    token = Splitter_getToken(splitter);
    block->e = atoi(token); // end 1-based already
    while ((token = Splitter_getToken(splitter)) != NULL){
	    block->attrbsLen += 1;
	    if(block->attrbsLen == 1){
		    block->attrbs = malloc(1 * sizeof(char*));
	    }
	    else{ // increase the size of the attrbs if there is more attrbs
		    realloc(block->attrbs, block->attrbsLen * sizeof(char*));
	    }
	    // save the currect attrb
	    block->attrbs[block->attrbsLen - 1] = malloc(strlen(token) + 1);
	    strcpy(block->attrbs[block->attrbsLen - 1], token);
    }
    Splitter_destruct(splitter);
    free(line);
    return 1;
}

int readNextBlockCov(FILE* fileReader, Block_t* block){
    size_t len = 0;
    char* line = NULL;
    char* token;
    int start, end;
    // free block attrbs to fill it with new ones
    for(int i = 0; i < block->attrbsLen; i++){
                 free(block->attrbs[i]);
    }
    free(block->attrbs);
    block->attrbs = NULL;
    block->attrbsLen = 0;
    ssize_t read = getline(&line, &len, fileReader);
    line[read-1]='\0';
    // if this is the end of the file
    if (read == -1){
         block->ctg[0] = '\0';
	 block->ctgLen = 0;
         block->s = -1;
         block->e = -1;
         return 0;
    }
    if (line[0] == '>') { // contig name and size is after '>'
	Splitter* splitter = Splitter_construct(line, ' ');
	token = Splitter_getToken(splitter);
        strcpy(block->ctg, token + 1); // skip '>' and copy
	token = Splitter_getToken(splitter);
	block->ctgLen = atoi(token);
	Splitter_destruct(splitter);
        readNextBlockCov(fileReader, block);
    }
    else {
	Splitter* splitter = Splitter_construct(line, '\t');
	token = Splitter_getToken(splitter);
        block->s = atoi(token); // start 1-based already
	token = Splitter_getToken(splitter);
	block->e = atoi(token); // end 1-based already
	while ((token = Splitter_getToken(splitter)) != NULL){
            block->attrbsLen += 1;
	    if(block->attrbsLen == 1){
                    block->attrbs = malloc(1 * sizeof(char*));
            }
            else{// increase the size of the attrbs if there is more attrbs
                    realloc(block->attrbs, block->attrbsLen * sizeof(char*));
            }
	    // save the current attrb
            block->attrbs[block->attrbsLen - 1] = malloc(strlen(token) + 1);
            strcpy(block->attrbs[block->attrbsLen - 1], token);
	}
	Splitter_destruct(splitter);
    }
    free(line);
    return 1;
}
