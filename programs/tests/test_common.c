#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "sonLib.h"


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
    bool correct = strcmp(str, "1.00000e-03,1.00000e-01,2.33000e+00") == 0;
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

bool test_IntBinArray_constructFromFile(char *filePath) {
    IntBinArray *binArray = IntBinArray_constructFromFile(filePath);
    bool correct = true;
    // check some bin indices
    correct &= (IntBinArray_getBinIndex(binArray, 0) == 0);
    correct &= (IntBinArray_getBinIndex(binArray, 50) == 0);
    correct &= (IntBinArray_getBinIndex(binArray, 100) == 1);
    correct &= (IntBinArray_getBinIndex(binArray, 200) == 2);
    correct &= (IntBinArray_getBinIndex(binArray, 1000) == 2);
    // check some names
    correct &= (strcmp(IntBinArray_getBinNameByIndex(binArray, 0), "0-100") == 0);
    correct &= (strcmp(IntBinArray_getBinNameByIndex(binArray, 1), "100-200") == 0);
    correct &= (strcmp(IntBinArray_getBinNameByIndex(binArray, 2), "200<") == 0);
    // by value
    correct &= (strcmp(IntBinArray_getBinNameByValue(binArray, 0), "0-100") == 0);
    correct &= (strcmp(IntBinArray_getBinNameByValue(binArray, 1000), "200<") == 0);
    IntBinArray_destruct(binArray);
    return correct;

}

bool test_IntBinArray_getBinIndices(char *filePath) {
    IntBinArray *binArray = IntBinArray_constructFromFile(filePath);
    bool correct = true;

    int len1 = 0;
    int *indices1 = IntBinArray_getBinIndices(binArray, 50, &len1);
    correct &= len1 == 2;
    correct &= indices1[0] == 0;
    correct &= indices1[1] == 1;
    free(indices1);

    int len2 = 0;
    int *indices2 = IntBinArray_getBinIndices(binArray, 80, &len2);
    correct &= len2 == 1;
    correct &= indices2[0] == 1;
    free(indices2);

    IntBinArray_destruct(binArray);
    return correct;

}


bool test_IntBinArray_constructSingleBin() {
    IntBinArray *binArray = IntBinArray_constructSingleBin(0, 1e9, "all");
    bool correct = true;
    // check some bin indices
    correct &= (IntBinArray_getBinIndex(binArray, 0) == 0);
    correct &= (IntBinArray_getBinIndex(binArray, 500) == 0);
    // check some names
    correct &= (strcmp(IntBinArray_getBinNameByIndex(binArray, 0), "all") == 0);
    // by value
    correct &= (strcmp(IntBinArray_getBinNameByValue(binArray, 500), "all") == 0);
    IntBinArray_destruct(binArray);
    return correct;

}

bool test_Double_equality2DArray() {
    double **array1 = Double_construct2DArray(2, 3);
    double **array2 = Double_construct2DArray(2, 3);
    double **array3 = Double_construct2DArray(2, 3);

    array1[0][0] = 1.0;
    array1[1][2] = 2.0;

    array2[0][0] = 1.0;
    array2[1][2] = 2.0;

    array3[0][0] = 1.0;
    array3[1][2] = 3.0;

    bool correct = true;

    correct &= Double_equality2DArray(array1, array2, 2, 3) == true;
    correct &= Double_equality2DArray(array1, array3, 2, 3) == false;

    return correct;

}


bool test_Double_getString2DArray() {
    double **array = Double_construct2DArray(2, 3);
    array[0][0] = 1.0;
    array[1][2] = 2.0;

    char *str = Double_getString2DArray(array, 2, 3, "%.2f", 10, '\t');

    bool correct = strcmp(str, "1.00\t0.00\t0.00\n0.00\t0.00\t2.00") == 0;

    free(str);
    return correct;
}


int main(int argc, char *argv[]) {
    char test1Path[1000] = "tests/test_files/common/test_common_1.txt";
    char test2Path[1000] = "tests/test_files/common/test_common_bin_array.txt";
    char test3Path[1000] = "tests/test_files/common/test_common_bin_array_2.txt";

    bool allTestsPassed = true;

    // test 1
    bool test1Passed = test_Splitter_parseLinesIntoList(test1Path);
    printf("[common] Test parsing file into stList:");
    printf(test1Passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");
    allTestsPassed &= test1Passed;

    // test 2
    bool test2Passed = test_String_joinDoubleArray();
    printf("[common] Test String_joinDoubleArray:");
    printf(test2Passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");
    allTestsPassed &= test2Passed;

    // test 3
    bool test3Passed = test_String_joinIntArray();
    printf("[common] Test String_joinIntArray:");
    printf(test3Passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");
    allTestsPassed &= test3Passed;


    // test 4
    bool test4Passed = test_String_joinStringArray();
    printf("[common] Test String_joinStringArray:");
    printf(test4Passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");
    allTestsPassed &= test4Passed;

    // test 5
    bool test5Passed = test_Splitter_getIntArray();
    printf("[common] Test Splitter_getIntArray:");
    printf(test5Passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");
    allTestsPassed &= test5Passed;

    // test 6
    bool test6Passed = test_Splitter_getDoubleArray();
    printf("[common] Test Splitter_getDoubleArray:");
    printf(test6Passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");
    allTestsPassed &= test6Passed;

    // test 7
    bool test7Passed = test_IntBinArray_constructFromFile(test2Path);
    printf("[common] Test IntBinArray_constructFromFile:");
    printf(test7Passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");
    allTestsPassed &= test7Passed;

    // test 8
    bool test8Passed = test_IntBinArray_constructSingleBin();
    printf("[common] Test IntBinArray_constructSingleBin:");
    printf(test8Passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");
    allTestsPassed &= test8Passed;

    // test 9
    bool test9Passed = test_Double_equality2DArray();
    printf("[common] Test Double_equality2DArray:");
    printf(test9Passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");
    allTestsPassed &= test9Passed;

    // test 10
    bool test10Passed = test_Double_getString2DArray();
    printf("[common] Test Double_getString2DArray:");
    printf(test10Passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");
    allTestsPassed &= test10Passed;

    // test 11
    bool test11Passed = test_IntBinArray_getBinIndices(test3Path);
    printf("[common] Test IntBinArray_getBinIndices:");
    printf(test11Passed ? "\x1B[32m OK \x1B[0m\n" : "\x1B[31m FAIL \x1B[0m\n");
    allTestsPassed &= test11Passed;


    if (allTestsPassed)
        return 0;
    else
        return 1;
}
