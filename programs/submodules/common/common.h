#include <stdint.h>
#include <string.h>
#include <sys/stat.h>   // stat
#include <stdbool.h>    // bool type
#include <assert.h>
#include <stdio.h>
#include <time.h>

#ifndef COMMON_H
#define COMMON_H

#define ARRAY_SIZE(arr) (sizeof((arr)) / sizeof((arr)[0]))
#define TIMESTAMP_SIZE 40
#define PI 3.14159

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

int getFirstIndexWithNonZeroBitFromRight(int32_t a);

char* get_timestamp();

bool file_exists (char *filename);

char* read_whole_file(char* file_path, long* length_ptr, char* mode);

char* copyString(char* str);

void removeSpacesInPlace(char* s);

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


double *Double_construct1DArray(int length);
double **Double_construct2DArray(int length1, int length2);
void Double_destruct1DArray(double *array);
void Double_destruct2DArray(double **array, int length1);
void Double_fill1DArray(double *array, int length, double value);
void Double_fill2DArray(double **array, int length1, int length2, double value);
double Double_sum1DArray(double *array, int length);
double Double_sum2DArray(double **array, int length1, int length2);
double *Double_copy1DArray(double *src, int length);
double **Double_copy2DArray(double **src, int length1, int length2);
void Double_multiply1DArray(double *array, int length, double factor);
void Double_multiply2DArray(double **array, int length1, int length2, double factor);
double Double_getMaxValue1DArray(double *array, int length);
double Double_getMaxValue2DArray(double **array, int length1, int length2);


int *Int_construct1DArray(int length);
int **Int_construct2DArray(int length1, int length2);
void Int_destruct1DArray(int *array);
void Int_destruct2DArray(int **array, int length1);
void Int_fill1DArray(int *array, int length, int value);
void Int_fill2DArray(int **array, int length1, int length2, int value);
int Int_sum1DArray(int *array, int length);
int Int_sum2DArray(int **array, int length1, int length2);
int *Int_copy1DArray(int *src, int length);
int **Int_copy2DArray(int **src, int length1, int length2);
void Int_multiply1DArray(int *array, int length, int factor);
void Int_multiply2DArray(int **array, int length1, int length2, int factor);
int Int_getMaxValue1DArray(int *array, int length);
int Int_getMaxValue2DArray(int **array, int length1, int length2);

uint8_t maxCharArray(uint8_t* a, int len);

uint8_t minCharArray(uint8_t* a, int len);

int maxIntArray(int* a, int len);

int minIntArray(int* a, int len);

// Takes a string and a delimiter that separates the integer numbers in the string
// Returns the array of numbers and sets the array size
int *Splitter_getIntArray(char *str, char delimiter, int *arraySize);

// Takes a string and a delimiter that separates the double numbers in the string
// Returns the array of numbers and sets the array size
double *Splitter_getDoubleArray(char *str, char delimiter, double *arraySize);

char* String_copy(const char* src);
char *String_joinDoubleArray(double *array, int length, char delimiter);
char *String_joinStringArray(const char** array, int elementMaxSize, int length, char delimiter);

#endif


