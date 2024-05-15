#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include "sonLib.h"


#include <sys/resource.h>
#include <sys/time.h>

// System functions are copied from https://github.com/chhylp123/hifiasm/blob/master/sys.cpp

double System_getCpuTime(void) {
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);
    return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

double System_getRealTimePoint(void) {
    struct timeval tp;
    struct timezone tzp;
    gettimeofday(&tp, &tzp);
    return tp.tv_sec + tp.tv_usec * 1e-6;
}


long System_getPeakRSS(void) {
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);
#ifdef __linux__
    return r.ru_maxrss * 1024;
#else
    return r.ru_maxrss;
#endif
}

double System_getPeakRSSInGB(void) {
    return System_getPeakRSS() / 1073741824.0;
}

double System_getCpuUsage(double cputime, double realtime) {
    return (cputime + 1e-9) / (realtime + 1e-9);
}


char *extractFileExtension(char *filePath) {
    int len = strlen(filePath);
    int i = len - 1;
    for (; 0 <= i; i--) {
        if (filePath[i] == '.') {
            if (strcmp(filePath + i, ".gz") != 0 &&
                strcmp(filePath + i, ".tar") != 0 &&
                strcmp(filePath + i, ".tar.gz") != 0 &&
                strcmp(filePath + i, ".zip") != 0)
                break;
        }
    }
    char *extension = malloc(len - i);
    strcpy(extension, filePath + i + 1);
    return extension;
}

int getFirstIndexWithNonZeroBitFromRight(int32_t a) {
    int n = 0;
    for (; n < 32; n++) {
        if ((a >> n) & 1) {
            return n;
        }
    }
    return -1;
}

char *copyString(char *str) {
    char *copy = (char *) malloc(strlen(str) + 1);
    strcpy(copy, str);
    return copy;
}

void removeSpacesInPlace(char *s) {
    //https://stackoverflow.com/questions/1726302/remove-spaces-from-a-string-in-c
    char *d = s;
    do {
        while (*d == ' ') {
            ++d;
        }
    } while (*s++ = *d++);
}

char *read_whole_file(char *file_path, long *length_ptr, char *mode) {
    char *buffer = 0;
    long length;
    FILE *f = fopen(file_path, mode);

    if (f) {
        fseek(f, 0, SEEK_END);
        length = ftell(f);
        fseek(f, 0, SEEK_SET);
        buffer = malloc(length);
        if (buffer) {
            fread(buffer, 1, length, f);
        }
        fclose(f);
    } else {
        return NULL;
        *length_ptr = 0;
    }
    *length_ptr = length;
    return buffer;
}

char *get_timestamp() {
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

bool folder_exists(char *folderpath) {
    DIR *dir = opendir(folderpath);
    if (dir) {
        closedir(dir);
        return true;
    } else {
        return false;
    }
}


int min(int a, int b) {
    return a < b ? a : b;
}

int max(int a, int b) {
    return b < a ? a : b;
}

double *Double_construct1DArray(int length) {
    double *array = (double *) malloc(length * sizeof(double));
    Double_fill1DArray(array, length, 0.0);
    return array;
}

double **Double_construct2DArray(int length1, int length2) {
    double **array = malloc(length1 * sizeof(double *));
    for (int i = 0; i < length1; i++) {
        array[i] = Double_construct1DArray(length2);
    }
    return array;
}

void Double_destruct1DArray(double *array) {
    free(array);
}

void Double_destruct2DArray(double **array, int length1) {
    for (int i = 0; i < length1; i++) {
        Double_destruct1DArray(array[i]);
    }
    free(array);
}


void Double_fill1DArray(double *array, int length, double value) {
    for (int i = 0; i < length; i++) {
        array[i] = value;
    }
}

void Double_fill2DArray(double **array, int length1, int length2, double value) {
    for (int i = 0; i < length1; i++) {
        Double_fill1DArray(array[i], length2, value);
    }
}

double Double_sum1DArray(double *array, int length) {
    double sum = 0.0;
    for (int i = 0; i < length; i++) {
        sum += array[i];
    }
    return sum;
}

double Double_sum2DArray(double **array, int length1, int length2) {
    double sum = 0.0;
    for (int i = 0; i < length1; i++) {
        sum += Double_sum1DArray(array[i], length2);
    }
    return sum;
}

double *Double_copy1DArray(double *src, int length) {
    double *dest = Double_construct1DArray(length);
    memcpy(dest, src, length * sizeof(double));
    return dest;
}

double **Double_copy2DArray(double **src, int length1, int length2) {
    double **dest = malloc(length1 * sizeof(double *));
    for (int i = 0; i < length1; i++) {
        dest[i] = Double_copy1DArray(src[i], length2);
    }
    return dest;
}


void Double_multiply1DArray(double *array, int length, double factor) {
    for (int i = 0; i < length; i++) {
        array[i] *= factor;
    }
}

void Double_multiply2DArray(double **array, int length1, int length2, double factor) {
    for (int i = 0; i < length1; i++) {
        Double_multiply1DArray(array[i], length2, factor);
    }
}

double Double_getMaxValue1DArray(double *array, int length) {
    assert(length > 0);
    double maxValue = array[0];
    for (int i = 0; i < length; i++) {
        if (maxValue < array[i]) {
            maxValue = array[i];
        }
    }
    return maxValue;
}

double Double_getMaxValue2DArray(double **array, int length1, int length2) {
    assert(length1 > 0 && length2 > 0);
    double maxValue = array[0][0];
    for (int i = 0; i < length1; i++) {
        double maxInRow = Double_getMaxValue1DArray(array[i], length2);
        if (maxValue < maxInRow) {
            maxValue = maxInRow;
        }
    }
    return maxValue;
}

int Double_getArgMaxIndex1DArray(double *array, int length) {
    assert(length > 0);
    double maxValue = array[0];
    int index = 0;
    for (int i = 0; i < length; i++) {
        if (maxValue < array[i]) {
            maxValue = array[i];
            index = i;
        }
    }
    return index;
}


int *Int_construct1DArray(int length) {
    int *array = malloc(length * sizeof(int));
    Int_fill1DArray(array, length, 0);
    return array;
}

int **Int_construct2DArray(int length1, int length2) {
    int **array = malloc(length1 * sizeof(int *));
    for (int i = 0; i < length1; i++) {
        array[i] = Int_construct1DArray(length2);
    }
    return array;
}

void Int_fill1DArray(int *array, int length, int value) {
    for (int i = 0; i < length; i++) {
        array[i] = value;
    }
}

void Int_fill2DArray(int **array, int length1, int length2, int value) {
    for (int i = 0; i < length1; i++) {
        Int_fill1DArray(array[i], length2, value);
    }
}

int Int_sum1DArray(int *array, int length) {
    int sum = 0;
    for (int i = 0; i < length; i++) {
        sum += array[i];
    }
    return sum;
}

int Int_sum2DArray(int **array, int length1, int length2) {
    int sum = 0;
    for (int i = 0; i < length1; i++) {
        sum += Int_sum1DArray(array[i], length2);
    }
    return sum;
}

void Int_destruct1DArray(int *array) {
    free(array);
}

void Int_destruct2DArray(int **array, int length1) {
    for (int i = 0; i < length1; i++) {
        Int_destruct1DArray(array[i]);
    }
}

int *Int_copy1DArray(int *src, int length) {
    int *dest = Int_construct1DArray(length);
    memcpy(dest, src, length * sizeof(int));
    return dest;
}

int **Int_copy2DArray(int **src, int length1, int length2) {
    int **dest = malloc(length1 * sizeof(int *));
    for (int i = 0; i < length1; i++) {
        dest[i] = Int_copy1DArray(src[i], length2);
    }
    return dest;
}

void Int_multiply1DArray(int *array, int length, int factor) {
    for (int i = 0; i < length; i++) {
        array[i] *= factor;
    }
}

void Int_multiply2DArray(int **array, int length1, int length2, int factor) {
    for (int i = 0; i < length1; i++) {
        Int_multiply1DArray(array[i], length2, factor);
    }
}

int Int_getMaxValue1DArray(int *array, int length) {
    assert(length > 0);
    int maxValue = array[0];
    for (int i = 0; i < length; i++) {
        if (maxValue < array[i]) {
            maxValue = array[i];
        }
    }
    return maxValue;

}

int Int_getModeValue1DArray(int *array, int length, int minValue, int maxValue) {
    int numberOfPossibleValues = maxValue - minValue + 1;
    int *counts = Int_construct1DArray(numberOfPossibleValues);

    for (int i = 0; i < length; i++) {
        int index = array[i] < minValue ? 0 : array[i] - minValue;
        index = numberOfPossibleValues <= index ? numberOfPossibleValues - 1 : index;
        counts[index] += 1;
    }
    int modeValue = minValue;
    int maxCount = counts[0]; // 0 is for minValue
    for (int i = 1; i < numberOfPossibleValues; i++) {
        if (maxCount < counts[i]) {
            modeValue = minValue + i;
            maxCount = counts[i];
        }
    }
    Int_destruct1DArray(counts);
    return modeValue;
}


int Int_getMaxValue2DArray(int **array, int length1, int length2) {
    assert(length1 > 0 && length2 > 0);
    int maxValue = array[0][0];
    for (int i = 0; i < length1; i++) {
        int maxInRow = Int_getMaxValue1DArray(array[i], length2);
        if (maxValue < maxInRow) {
            maxValue = maxInRow;
        }
    }
    return maxValue;
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
    return splitter;
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

int *Splitter_getIntArray(char *str, char delimiter, int *arraySize) {
    Splitter *splitter = Splitter_construct(str, delimiter);
    char *token;
    int i = 0;
    int *intArray = NULL;
    while ((token = Splitter_getToken(splitter)) != NULL) {
        intArray = (int *) realloc(intArray, (i + 1) * sizeof(int));
        intArray[i] = atoi(token);
        i++;
    }
    Splitter_destruct(splitter);
    *arraySize = i;
    return intArray;
}

double *Splitter_getDoubleArray(char *str, char delimiter, int *arraySize) {
    Splitter *splitter = Splitter_construct(str, delimiter);
    char *token;
    int i = 0;
    double *doubleArray = NULL;
    while ((token = Splitter_getToken(splitter)) != NULL) {
        doubleArray = (double *) realloc(doubleArray, (i + 1) * sizeof(double));
        doubleArray[i] = atof(token);
        i++;
    }
    Splitter_destruct(splitter);
    *arraySize = i;
    return doubleArray;
}


char *String_copy(const char *src) {
    char *dest = malloc(strlen(src) + 1);
    strcpy(dest, src);
    return dest;
}

char *String_joinDoubleArray(double *array, int length, char delimiter) {
    char *str = malloc(length * 10 * sizeof(char));
    str[0] = '\0';
    for (int i = 0; i < length - 1; i++) {
        sprintf(str + strlen(str), "%.2e%c", array[i], delimiter);
    }
    sprintf(str + strlen(str), "%.2e", array[length - 1]);
    str[strlen(str)] = '\0';
    return str;
}

char *String_joinIntArray(int *array, int length, char delimiter) {
    char *str = malloc(length * 10 * sizeof(char));
    str[0] = '\0';
    for (int i = 0; i < length - 1; i++) {
        sprintf(str + strlen(str), "%d%c", array[i], delimiter);
    }
    sprintf(str + strlen(str), "%d", array[length - 1]);
    str[strlen(str)] = '\0';
    return str;
}


char *String_joinStringArray(const char **array, int elementMaxSize, int length, char delimiter) {
    char *str = malloc(length * elementMaxSize * sizeof(char));
    str[0] = '\0';
    for (int i = 0; i < length - 1; i++) {
        sprintf(str + strlen(str), "%s%c", array[i], delimiter);
    }
    sprintf(str + strlen(str), "%s", array[length - 1]);
    return str;
}


stList *Splitter_parseLinesIntoList(const char *filepath) {
    stList *lines = stList_construct3(0, free);
    FILE *fp = fopen(filepath, "r");
    if (fp == NULL) {
        fprintf(stderr, "Error: Unable to open %s\n", filepath);
    }
    int len;
    int read;
    char *token;
    char *line = malloc(1000);
    while ((read = getline(&line, &len, fp)) > 0) {
        //fprintf(stderr, "%s %d\n",line, read);
        if (0 < read && line[read - 1] == '\n') line[read - 1] = '\0';
        Splitter *splitter = Splitter_construct(line, ' ');
        // get contig name
        token = Splitter_getToken(splitter);
        stList_append(lines, copyString(token));
        Splitter_destruct(splitter);
    }
    fclose(fp);
    free(line);
    return lines;
}


IntBinArray *IntBinArray_constructSingleBin(int start, int end, char *name) {
    IntBinArray *binArray = malloc(sizeof(IntBinArray));
    binArray->starts = Int_construct1DArray(1);
    binArray->starts[0] = start;
    binArray->ends = Int_construct1DArray(1);
    binArray->ends[0] = end;
    binArray->names = stList_construct3(0, free);
    stList_append(binArray->names, copyString(name));
    binArray->numberOfBins = 1;
    return binArray;
}

// it takes a path to a tsv file (tab-delimited)
// like this
//#start    end name
//0 1000    0-1Kb
//1000  2000    1-2Kb
//2000  10000   2Kb<
IntBinArray *IntBinArray_constructFromFile(const char *filePath){
    IntBinArray *binArray = malloc(sizeof(IntBinArray));
    stList *lines = Splitter_parseLinesIntoList(filePath);
    int numberOfBins = 0;
    for(int i=0; i < stList_length(lines); i++) {
        // get line
        char *line = stList_get(lines, i);
        if (line[0] != '#') numberOfBins += 1;
    }
    // set attributes
    binArray->names = stList_construct3(0, free);
    binArray->starts = Int_construct1DArray(numberOfBins);
    binArray->ends = Int_construct1DArray(numberOfBins);
    binArray->numberOfBins = numberOfBins;
    int binIndex = 0;
    for(int i=0; i < stList_length(lines); i++){
        // get line
        char *line = stList_get(lines, i);
        if(line[0] == '#') continue;
        Splitter *splitter = Splitter_construct(line,'\t');
        // fetch columns using atof to handle scientific notation here
        int start = (int) atof(Splitter_getToken(splitter)); //first column
        int end = (int) atof(Splitter_getToken(splitter)); // second column
        char *name = copyString(Splitter_getToken(splitter)); // third column
        // set attributes
        stList_append(binArray->names, name);
        binArray->starts[binIndex] = start;
        binArray->ends[binIndex] = end;
        binIndex += 1; // bin index might be different from i due to header
        Splitter_destruct(splitter);
    }
    stList_destruct(lines);
    return binArray;
}

void IntBinArray_checkBins(IntBinArray *binArray) {
    for(int i=1; i < binArray->numberOfBins; i++){
        if(binArray->starts[i] < binArray->ends[i-1]) {
            fprintf(stderr, "[%s] Error: Bins cannot have overlap with each other.\n");
            exit(EXIT_FAILURE);
        }
    }
}

int IntBinArray_getBinIndex(IntBinArray *binArray, int value) {
    for(int i=0; i < binArray->numberOfBins; i++){
        if(binArray->starts[i] <= value && value < binArray->ends[i]) return i;
    }
    fprintf(stderr, "[%s] Error: Value %d does not fit into the given bins.\n", get_timestamp());
    exit(EXIT_FAILURE);
}

char *IntBinArray_getBinNameByIndex(IntBinArray *binArray, int binIndex) {
    return (char *)stList_get(binArray->names, binIndex);
}

char *IntBinArray_getBinNameByValue(IntBinArray *binArray, int value) {
    int binIndex = IntBinArray_getBinIndex(binArray, value);
    return IntBinArray_getBinNameByIndex(binArray, binIndex);
}

void IntBinArray_destruct(IntBinArray *binArray) {
    free(binArray->starts);
    free(binArray->ends);
    if (binArray->names != NULL) stList_destruct(binArray->names);
    free(binArray);
}

