#ifndef COMMON_H
#define COMMON_H
#include <stdint.h>
#include <string.h>
#include <sys/stat.h>   // stat
#include <stdbool.h>    // bool type
#include <assert.h>

char* get_timestamp();

bool file_exists (char *filename);

typedef struct Splitter{
	char* str;
	int offset;
	char* token;
	char delimiter;
} Splitter;

Splitter* Splitter_construct(char* str, char delimiter);

void Splitter_destruct(Splitter* splitter);

char* Splitter_getToken(Splitter* splitter);

int min(int a, int b);

int max(int a, int b);

uint8_t maxCharArray(uint8_t* a, int len);

uint8_t minCharArray(uint8_t* a, int len);

int maxIntArray(int* a, int len);

int minIntArray(int* a, int len);

#endif


