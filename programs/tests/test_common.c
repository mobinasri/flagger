#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "sonLib.h"

/*
char *String_joinDoubleArray(double *array, int length, char delimiter);
char *String_joinIntArray(int *array, int length, char delimiter);
char *String_joinStringArray(const char** array, int elementMaxSize, int length, char delimiter);
stList *Splitter_parseLinesIntoList(const char *filepath);
*/

bool test_Splitter_parseLinesIntoList(const char *filepath) {
    bool correct = true;
    stList *list = Splitter_parseLinesIntoList(filepath);
    if (stList_length(list) != 2) return false;
    char *string1 = stList_get(list, 0);
    char *string2 = stList_get(list, 1);
    correct &= (strcmp(string1, "line1") == 0);
    correct &= (strcmp(string2, "line2") == 0);
    stList_destruct(list);
    return correct;
}


bool test_String_joinDoubleArray() {
    double x[3];
    x[0] = 1e-3;
    x[1] = 0.1;
    x[2] = 2.33;
    char *str = String_joinDoubleArray(x, 3, ',');
    bool correct = strcmp(str, "1.00e-03,1.00e-01,2.33e+00") == 0;
    free(str);
    return correct;
}

bool test_String_joinIntArray() {
    int x[3];
    x[0] = 12;
    x[1] = 0;
    x[2] = 123;
    char *str = String_joinIntArray(x, 3, ',');
    bool correct = strcmp(str, "12,0,123") == 0;
    free(str);
    return correct;
}

bool test_String_joinStringArray() {
    char **x = malloc(3 * sizeof(char *));
    x[0] = copyString("s1");
    x[1] = copyString("s2");
    x[2] = copyString("s3");
    char *str = String_joinStringArray(x, 20, 3, ',');
    bool correct = strcmp(str, "s1,s2,s3") == 0;
    free(str);
    free(x[0]);
    free(x[1]);
    free(x[2]);
    free(x);
    return correct;
}


bool test_Splitter_getDoubleArray() {
    int len = -1;
    char *str = copyString("1e-3,0.1,2.33");
    double *array = Splitter_getDoubleArray(str, ',', &len);
    if (len != 3) return false;
    bool correct = (array[0] == 1e-3) && (array[1] == 0.1) && (array[2] == 2.33);
    free(str);
    free(array);
    return correct;

}

bool test_Splitter_getIntArray() {
    int len = -1;
    char *str = copyString("12,0,123");
    int *array = Splitter_getIntArray(str, ',', &len);
    if (len != 3) return false;
    bool correct = (array[0] == 12) && (array[1] == 0) && (array[2] == 123);
    free(str);
    free(array);
    return correct;

}


int main(int argc, char *argv[]) {
    char test1Path[1000] = "tests/test_files/test_common_1.txt";
    if (test_Splitter_parseLinesIntoList(test1Path) == true) {
        fprintf(stderr, "test_Splitter_parseLinesIntoList for %s passed!\n", test1Path);
    } else {
        fprintf(stderr, "test_Splitter_parseLinesIntoList for %s failed!\n", test1Path);
        return 1;
    }

    if (test_String_joinDoubleArray()) {
        fprintf(stderr, "test_String_joinDoubleArray passed!\n");
    } else {
        fprintf(stderr, "test_String_joinDoubleArray failed!\n");
        return 1;
    }

    if (test_String_joinIntArray()) {
        fprintf(stderr, "test_String_joinIntArray passed!\n");
    } else {
        fprintf(stderr, "test_String_joinIntArray failed!\n");
        return 1;
    }

    if (test_String_joinStringArray()) {
        fprintf(stderr, "test_String_joinStringArray passed!\n");
    } else {
        fprintf(stderr, "test_String_joinStringArray failed!\n");
        return 1;
    }

    if (test_Splitter_getIntArray()) {
        fprintf(stderr, "test_Splitter_getIntArray passed!\n");
    } else {
        fprintf(stderr, "test_Splitter_getIntArray failed!\n");
        return 1;
    }

    if (test_Splitter_getDoubleArray()) {
        fprintf(stderr, "test_Splitter_getDoubleArray passed!\n");
    } else {
        fprintf(stderr, "test_Splitter_getDoubleArray failed!\n");
        return 1;
    }
    return 0;

}
