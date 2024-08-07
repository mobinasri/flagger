#include "common.h"
#include <stdio.h>
#include <stdlib.h>


int getFirstIndexWithNonZeroBitFromRight(int32_t a){
    int n=0;
    for(;n < 32; n++){
        if ((a >> n) & 1){
            return n;
        }
    }
    return -1;
}

char* copyString(char* str){
    char* copy  = (char*) malloc(strlen(str) + 1);
    strcpy(copy, str);
    return copy;
}

void removeSpacesInPlace(char* s) {
    //https://stackoverflow.com/questions/1726302/remove-spaces-from-a-string-in-c
    char* d = s;
    do {
        while (*d == ' ') {
            ++d;
        }
    } while (*s++ = *d++);
}

char* read_whole_file(char* file_path, long* length_ptr, char* mode){
    char * buffer = 0;
    long length;
    FILE * f = fopen (file_path, mode);

    if (f)
    {
      fseek (f, 0, SEEK_END);
      length = ftell (f);
      fseek (f, 0, SEEK_SET);
      buffer = malloc (length);
      if (buffer)
      {
        fread (buffer, 1, length, f);
      }
      fclose (f);
    }
    else{
        return NULL;
        *length_ptr = 0;
    }
    *length_ptr = length;
    return buffer;
}

char* get_timestamp() {
    static char timestamp[TIMESTAMP_SIZE + 1];  // static variable to hold the timestamp string
    time_t t = time(NULL);
    struct tm *tm = localtime(&t);
    snprintf(timestamp, TIMESTAMP_SIZE + 1, "%04d-%02d-%02d %02d:%02d:%02d",
             tm->tm_year + 1900, tm->tm_mon + 1, tm->tm_mday,
             tm->tm_hour, tm->tm_min, tm->tm_sec);
    return timestamp;
}

bool file_exists(char *filename) {
    struct stat buffer;
    return (stat(filename, &buffer) == 0);
}

int min(int a, int b) {
    return a < b ? a : b;
}

int max(int a, int b) {
    return b < a ? a : b;
}

double *Double_construct1DArray(int length){
    double *array = (double*) malloc(length * sizeof(double));
    Double_fill1DArray(array, length, 0.0);
    return array;
}

double **Double_construct2DArray(int length1, int length2){
    double **array = malloc(length1 * sizeof(double*));
    for(int i = 0; i < length1; i++) {
        array[i] = Double_construct1DArray(length2);
    }
    return array;
}

void Double_destruct1DArray(double *array){
    free(array);
}

void Double_destruct2DArray(double **array, int length1){
    for(int i = 0; i < length1; i++) {
        Double_destruct1DArray(array[i]);
    }
    free(array);
}


void Double_fill1DArray(double *array, int length, double value){
    for (int i = 0; i < length; i++) {
        array[i] = value;
    }
}

void Double_fill2DArray(double **array, int length1, int length2, double value){
    for (int i = 0; i < length1; i++) {
        Double_fill1DArray(array[i], length2, value);
    }
}

double Double_sum1DArray(double *array, int length){
    double sum = 0.0;
    for (int i = 0; i < length; i++) {
        sum += array[i];
    }
}

double Double_sum2DArray(double **array, int length1, int length2){
    double sum = 0.0;
    for (int i = 0; i < length1; i++) {
        sum += Double_sum1DArray(array[i], length2);
    }
}

double *Double_copy1DArray(double *src, int length){
    double *dest = Double_construct1DArray(length);
    memcpy(src, dest, length * sizeof(double));
    return dest;
}

double **Double_copy2DArray(double **src, int length1, int length2){
    double **dest = malloc(length1 * sizeof(double*));
    for (int i = 0; i < length1; i++) {
        dest[i] = Double_copy1DArray(src[i], length2);
    }
    return dest;
}

int *Int_construct1DArray(int length){
    int *array = malloc(length * sizeof(int));
    Int_fill1DArray(array, length, 0);
    return array;
}

int **Int_construct2DArray(int length1, int length2){
    int **array = malloc(length1 * sizeof(int*));
    for(int i=0; i < length1; i++) {
        array[i] = Int_construct1DArray(length2);
    }
    return array;
}

void Int_fill1DArray(int *array, int length, int value){
    for (int i = 0; i < length; i++) {
        array[i] = value;
    }
}

void Int_fill2DArray(int **array, int length1, int length2, int value){
    for (int i = 0; i < length1; i++) {
        Int_fill1DArray(array[i], length2, value);
    }
}

int Int_sum1DArray(int *array, int length){
    int sum = 0;
    for (int i = 0; i < length; i++) {
        sum += array[i];
    }
    return sum;
}

int Int_sum2DArray(int **array, int length1, int length2){
    int sum = 0;
    for (int i = 0; i < length1; i++) {
        sum += Int_sum1DArray(array[i], length2);
    }
    return sum;
}

void Int_destruct1DArray(int *array){
    free(array);
}

void Int_destruct2DArray(int **array, int length1){
    for (int i = 0; i < length1; i++) {
        Int_destruct1DArray(array[i]);
    }
}

int *Int_copy1DArray(int *src, int length){
    int *dest = Int_construct1DArray(length);
    memcpy(src, dest, length * sizeof(int));
    return dest;
}

int **Int_copy2DArray(int **src, int length1, int length2){
    int **dest = malloc(length1 * sizeof(int*));
    for (int i = 0; i < length1; i++) {
        dest[i] = Int_copy1DArray(src[i], length2);
    }
    return dest;
}


uint8_t maxCharArray(uint8_t *a, int len) {
    assert(len > 0);
    uint8_t m = a[0];
    for (int i = 0; i < len; i++) {
        m = a[i] < m ? m : a[i];
    }
    return m;
}

uint8_t minCharArray(uint8_t *a, int len) {
    assert(len > 0);
    uint8_t m = a[0];
    for (int i = 0; i < len; i++) {
        m = m < a[i] ? m : a[i];
    }
    return m;
}

int maxIntArray(int *a, int len) {
    assert(len > 0);
    int m = a[0];
    for (int i = 0; i < len; i++) {
        m = a[i] < m ? m : a[i];
    }
    return m;
}

int minIntArray(int *a, int len) {
    assert(len > 0);
    int m = a[0];
    for (int i = 0; i < len; i++) {
        m = m < a[i] ? m : a[i];
    }
    return m;
}

Splitter *Splitter_construct(char *str, char delimiter) {
    Splitter *splitter = malloc(sizeof(Splitter));
    splitter->str = malloc((strlen(str) + 1) * sizeof(char));
    strcpy(splitter->str, str);
    splitter->token = malloc((strlen(str) + 1) * sizeof(char));
    splitter->delimiter = delimiter;
    splitter->offset = 0;
}

void Splitter_destruct(Splitter *splitter) {
    free(splitter->str);
    free(splitter->token);
    free(splitter);
}

char *Splitter_getToken(Splitter *splitter) {
    int i = splitter->offset;
    int j = 0;
    while (splitter->str[i] != '\0' && splitter->str[i] != splitter->delimiter) {
        splitter->token[j] = splitter->str[i];
        j++;
        i++;
    }
    splitter->offset = splitter->str[i] == splitter->delimiter ? i + 1 : i;
    splitter->token[j] = '\0';
    if (j == 0) { // end of the string
        free(splitter->token);
        splitter->token = NULL;
    }
    return splitter->token;
}

int *Splitter_getIntArray(char *str, char delimiter, int *arraySize){
    Splitter *splitter = Splitter_construct(str, delimiter);
    char *token;
    int i = 0;
    int *intArray = NULL;
    while ((token = Splitter_getToken(splitter)) != NULL) {
        intArray = (int*) realloc(intArray, (i + 1) * sizeof(int));
        intArray[i] = atoi(token);
        i++;
    }
    Splitter_destruct(splitter);
    *arraySize = i;
    return intArray;
}

double *Splitter_getDoubleArray(char *str, char delimiter, double *arraySize){
    Splitter *splitter = Splitter_construct(str, delimiter);
    char *token;
    int i = 0;
    double *doubleArray = NULL;
    while ((token = Splitter_getToken(splitter)) != NULL) {
        doubleArray = (double*) realloc(doubleArray, (i + 1) * sizeof(double));
        doubleArray[i] = atof(token);
        i++;
    }
    Splitter_destruct(splitter);
    *arraySize = i;
    return doubleArray;
}