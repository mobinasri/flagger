#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#ifndef DATA_TYPES_H
#define DATA_TYPES_H

typedef struct VectorChar{
        int dim; // Dimension
        uint8_t* data;
} VectorChar;

typedef struct VectorDouble{
        int dim; // Dimension
        double* data;
} VectorDouble;

typedef struct MatrixChar{
        int dim1; // Number of rows
        int dim2; // Number of columns
        uint8_t** data;
} MatrixChar;

typedef struct MatrixDouble{
        int dim1; // Number of rows
        int dim2; // Number of columns
        double** data;
} MatrixDouble;

// VectorChar functions

VectorChar* VectorChar_construct0(int dim);

VectorChar* VectorChar_construct1(int dim, uint8_t* data);

VectorChar** VectorChar_constructArray1D(int len, int dim);

VectorChar*** VectorChar_constructArray2D(int len1, int len2, int dim);

void VectorChar_setValue(VectorChar* vec, uint8_t value);

void VectorChar_divideByValue(VectorChar* vec, uint8_t value);

VectorChar* VectorChar_copy(VectorChar* src);

void VectorChar_copyInPlace(VectorChar* dest, VectorChar* src);

char* VectorChar_toString(VectorChar* vec);

void VectorChar_destruct(VectorChar* vec);

void VectorChar_destructArray1D(VectorChar** vecArray, int len);

void VectorChar_destructArray2D(VectorChar*** vecArray, int len1, int len2);

// VectorDouble functions

VectorDouble* VectorDouble_construct0(int dim);

VectorDouble* VectorDouble_construct1(int dim, double* data);

VectorDouble** VectorDouble_constructArray1D(int len, int dim);

VectorDouble*** VectorDouble_constructArray2D(int len1, int len2, int dim);

VectorDouble**** VectorDouble_constructArray3D(int len1, int len2, int len3, int dim);

void VectorDouble_setValue(VectorDouble* vec, double value);

void VectorDouble_divideByValue(VectorDouble* vec, double value);

VectorDouble* VectorDouble_copy(VectorDouble* src);

void VectorDouble_copyInPlace(VectorDouble* dest, VectorDouble* src);

char* VectorDouble_toString(VectorDouble* vec);

void VectorDouble_destruct(VectorDouble* vec);

void VectorDouble_destructArray1D(VectorDouble** vecArray, int len);

void VectorDouble_destructArray2D(VectorDouble*** vecArray, int len1, int len2);

void VectorDouble_destructArray3D(VectorDouble**** vecArray, int len1, int len2, int len3);


//MatrixChar functions

MatrixChar* MatrixChar_construct0(int dim1, int dim2);

MatrixChar* MatrixChar_construct1(int dim1, int dim2, uint8_t** data);

MatrixChar** MatrixChar_constructArray1D(int len, int dim1, int dim2);

MatrixChar*** MatrixChar_constructArray2D(int len1, int len2, int dim1, int dim2);

void MatrixChar_setValue(MatrixChar* mat, uint8_t value);

void MatrixChar_divideByValue(MatrixChar* mat, uint8_t value);

MatrixChar* MatrixChar_copy(MatrixChar* src);

void MatrixChar_copyInPlace(MatrixChar* dest, MatrixChar* src);

char* MatrixChar_toString(MatrixChar* mat);

void MatrixChar_destruct(MatrixChar* mat);

void MatrixChar_destructArray1D(MatrixChar** matArray, int len);

void MatrixChar_destructArray2D(MatrixChar*** matArray, int len1, int len2);

//MatrixDouble functions

MatrixDouble* MatrixDouble_construct0(int dim1, int dim2);

MatrixDouble* MatrixDouble_construct1(int dim1, int dim2, double** data);

MatrixDouble** MatrixDouble_constructArray1D(int len, int dim1, int dim2);

MatrixDouble*** MatrixDouble_constructArray2D(int len1, int len2, int dim1, int dim2);

void MatrixDouble_setValue(MatrixDouble* mat, double value);

void MatrixDouble_divideByValue(MatrixDouble* mat, double value);

MatrixDouble* MatrixDouble_copy(MatrixDouble* src);

void MatrixDouble_copyInPlace(MatrixDouble* dest, MatrixDouble* src);

char* MatrixDouble_toString(MatrixDouble* mat);

void MatrixDouble_destruct(MatrixDouble* mat);

void MatrixDouble_destructArray1D(MatrixDouble** matArray, int len);

void MatrixDouble_destructArray2D(MatrixDouble*** matArray, int len1, int len2);

MatrixDouble* MatrixDouble_parseFromFile(char* filePath, int dim1, int dim2);


#endif
