#ifndef COMMON_H
#define COMMON_H
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdbool.h>

#define ARRAY_SIZE(arr) (sizeof((arr)) / sizeof((arr)[0]))
#define TIMESTAMP_SIZE 40


//#define DEBUG

#ifdef DEBUG
//void DEBUG_PRINT(const char *, ...);
void DEBUG_PRINT(const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
}

#else
static inline void DEBUG_PRINT(const char *fmt, ...) {};
#endif

typedef struct Splitter{
	char* str;
	int offset;
	char* token;
	char delimiter;
} Splitter;


bool file_exists (char *filename);

char* get_timestamp();

char* copyString(char* str);

Splitter* Splitter_construct(char* str, char delimiter);

void Splitter_destruct(Splitter* splitter);

char* Splitter_getToken(Splitter* splitter);

int min(int a, int b);

int max(int a, int b);

uint8_t maxCharArray(uint8_t* a, int len);

uint8_t minCharArray(uint8_t* a, int len);

#endif


