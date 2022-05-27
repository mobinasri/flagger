#include <stdio.h>

#ifndef BLOCK_IT_H
#define BLOCK_IT_H

typedef enum Format_t {COV, BED} Format_t;

typedef struct Block_t{
	Format_t format;
        char ctg[50];
        int ctgLen;
        int s; // 1-based
        int e; // 1-based
        char** attrbs;
        int attrbsLen;
} Block_t;

// Given the format make an empty block
Block_t* Block_construct(Format_t format);

void Block_destruct(Block_t* block);

// Read next block in BED or COV
int Block_next(FILE* fileReader, Block_t* block);

#endif
