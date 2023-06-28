#include "data_types.h"
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"


// Functions for VectorChar

VectorChar* VectorChar_construct0(int dim){
	VectorChar* vec = (VectorChar*) malloc(sizeof(VectorChar));
	vec->dim = dim;
	vec->data = (uint8_t*) malloc(dim * sizeof(uint8_t));
	memset(vec->data, 0, dim);
	return vec;
}

VectorChar* VectorChar_construct1(int dim, uint8_t* data){
        VectorChar* vec = (VectorChar*) malloc(sizeof(VectorChar));
	vec->dim = dim;
        vec->data = (uint8_t*) malloc(dim * sizeof(uint8_t));
        memcpy(vec->data, data, dim);
        return vec;
}

VectorChar** VectorChar_constructArray1D(int len, int dim){
	VectorChar** vecArray = (VectorChar**) malloc(len * sizeof(VectorChar*));
	for(int i = 0; i < len; i++){
        	vecArray[i] = VectorChar_construct0(dim);
	}
        return vecArray;
}


VectorChar*** VectorChar_constructArray2D(int len1, int len2, int dim){
        VectorChar*** vecArray = (VectorChar***) malloc(len1 * sizeof(VectorChar**));
        for(int i = 0; i < len1; i++){
                vecArray[i] = VectorChar_constructArray1D(len2, dim);
        }
        return vecArray;
}


void VectorChar_destructArray1D(VectorChar** vecArray, int len){
        for(int i = 0; i < len; i++){
                VectorChar_destruct(vecArray[i]);
        }
	free(vecArray);
}

void VectorChar_destructArray2D(VectorChar*** vecArray, int len1, int len2){
        for(int i = 0; i < len1; i++){
                VectorChar_destructArray1D(vecArray[i], len2);
        }
        free(vecArray);
}

void VectorChar_setValue(VectorChar* vec, uint8_t value){
	for(int i = 0; i < vec->dim; i++){
		vec->data[i] = value;
	}
}

void VectorChar_divideByValue(VectorChar* vec, uint8_t value){
        assert(value != 0);
	float tmp;
	for(int i = 0; i < vec->dim; i++){
		tmp = (float) vec->data[i] / value;
                vec->data[i] = round(tmp);
        }
}

VectorChar* VectorChar_copy(VectorChar* src){
        VectorChar* dest = VectorChar_construct1(src->dim, src->data);
        return dest;
}

void VectorChar_copyInPlace(VectorChar* dest, VectorChar* src){
	assert(dest->dim == src->dim);
	for(int i = 0; i < dest->dim; i++){
                dest->data[i] = src->data[i];
        }
}

char* VectorChar_toString(VectorChar* vec){
	char tmp[50];
	char* str = malloc(vec->dim * 10);
	str[0] = '[';
	int offset = 1;
	for(int i=0; i< vec->dim; i++){
		sprintf(tmp, "%d", vec->data[i]);
		memcpy(str + offset, tmp, strlen(tmp));
		offset += strlen(tmp);
		if (i < vec->dim - 1){
			str[offset++] = ',';
			str[offset++] = '\n';
			str[offset++] = ' ';
		}
	}
	str[offset++] = ']';
	str[offset] = '\0';
	return str;
}

void VectorChar_destruct(VectorChar* vec){
        free(vec->data);
	vec->data = NULL;
	free(vec);
}


// Functions for VectorDouble

VectorDouble* VectorDouble_construct0(int dim){
        VectorDouble* vec = (VectorDouble*) malloc(sizeof(VectorDouble));
	vec->dim = dim;
        vec->data = (double*) malloc(dim * sizeof(double));
        memset(vec->data, 0, dim * sizeof(double));
        return vec;
}

VectorDouble* VectorDouble_construct1(int dim, double* data){
        VectorDouble* vec = (VectorDouble*) malloc(sizeof(VectorDouble));
	vec->dim = dim;
        vec->data = (double*) malloc(dim * sizeof(double));
        memcpy(vec->data, data, dim * sizeof(double));
        return vec;
}

VectorDouble** VectorDouble_constructArray1D(int len, int dim){
        VectorDouble** vecArray = (VectorDouble**) malloc(len * sizeof(VectorDouble*));
        for(int i = 0; i < len; i++){
                vecArray[i] = VectorDouble_construct0(dim);
        }
        return vecArray;
}


VectorDouble*** VectorDouble_constructArray2D(int len1, int len2, int dim){
        VectorDouble*** vecArray = (VectorDouble***) malloc(len1 * sizeof(VectorDouble**));
        for(int i = 0; i < len1; i++){
                vecArray[i] = VectorDouble_constructArray1D(len2, dim);
        }
        return vecArray;
}

VectorDouble**** VectorDouble_constructArray3D(int len1, int len2, int len3, int dim){
    VectorDouble**** vecArray = (VectorDouble****) malloc(len1 * sizeof(VectorDouble***));
    for(int i = 0; i < len1; i++){
        vecArray[i] = VectorDouble_constructArray2D(len2, len3, dim);
    }
    return vecArray;
}

void VectorDouble_destructArray1D(VectorDouble** vecArray, int len){
        for(int i = 0; i < len; i++){
                VectorDouble_destruct(vecArray[i]);
        }
        free(vecArray);
}

void VectorDouble_destructArray2D(VectorDouble*** vecArray, int len1, int len2){
        for(int i = 0; i < len1; i++){
                VectorDouble_destructArray1D(vecArray[i], len2);
        }
        free(vecArray);
}

void VectorDouble_destructArray3D(VectorDouble**** vecArray, int len1, int len2, int len3){
    for(int i = 0; i < len1; i++){
        VectorDouble_destructArray2D(vecArray[i], len2, len3);
    }
    free(vecArray);
}

void VectorDouble_setValue(VectorDouble* vec, double value){
        for(int i = 0; i < vec->dim; i++){
                vec->data[i] = value;
        }
}

void VectorDouble_divideByValue(VectorDouble* vec, double value){
        assert(value != 0.0);
        for(int i = 0; i < vec->dim; i++){
                vec->data[i] /= value;
        }
}

VectorDouble* VectorDouble_copy(VectorDouble* src){
        VectorDouble* dest = VectorDouble_construct1(src->dim, src->data);
        return dest;
}

void VectorDouble_copyInPlace(VectorDouble* dest, VectorDouble* src){
        assert(dest->dim == src->dim);
        for(int i = 0; i < dest->dim; i++){
                dest->data[i] = src->data[i];
        }
}

char* VectorDouble_toString(VectorDouble* vec){
        char tmp[50];
        char* str = malloc(vec->dim * 10);
        str[0] = '[';
        int offset = 1;
        for(int i=0; i< vec->dim; i++){
                sprintf(tmp, "%.2f", vec->data[i]);
                memcpy(str + offset, tmp, strlen(tmp));
                offset += strlen(tmp);
                if (i < vec->dim - 1){
                        str[offset++] = ',';
                        str[offset++] = '\n';
                        str[offset++] = ' ';
                }
        }
        str[offset++] = ']';
        str[offset] = '\0';
        return str;
}

void VectorDouble_destruct(VectorDouble* vec){
        free(vec->data);
        vec->data = NULL;
        free(vec);
}


// Functions for MatrixChar

MatrixChar* MatrixChar_construct0(int dim1, int dim2){
	MatrixChar* mat = (MatrixChar*) malloc(sizeof(MatrixChar));
	mat->dim1 = dim1;
	mat->dim2 = dim2;
	mat->data = (uint8_t**) malloc(dim1 * sizeof(uint8_t*));
	for(int i = 0; i < dim1; i++){
		mat->data[i] = (uint8_t*) malloc(dim2 * sizeof(uint8_t*));
		memset(mat->data[i], 0, dim2 * sizeof(uint8_t));
	}
	return mat;
}

MatrixChar* MatrixChar_construct1(int dim1, int dim2, uint8_t** data){
        MatrixChar* mat = (MatrixChar*) malloc(sizeof(MatrixChar));
	mat->dim1 = dim1;
	mat->dim2 = dim2;
        mat->data = (uint8_t**) malloc(dim1 * sizeof(uint8_t*));
        for(int i = 0; i < dim1; i++){
                mat->data[i] = (uint8_t*) malloc(dim2 * sizeof(uint8_t));
                memcpy(mat->data[i], data[i], dim2 * sizeof(uint8_t));
        }
        return mat;
}

MatrixChar** MatrixChar_constructArray1D(int len, int dim1, int dim2){
        MatrixChar** matArray = (MatrixChar**) malloc(len * sizeof(MatrixChar*));
        for(int i = 0; i < len; i++){
                matArray[i] = MatrixChar_construct0(dim1, dim2);
        }
        return matArray;
}


MatrixChar*** MatrixChar_constructArray2D(int len1, int len2, int dim1, int dim2){
        MatrixChar*** matArray = (MatrixChar***) malloc(len1 * sizeof(MatrixChar**));
        for(int i = 0; i < len1; i++){
                matArray[i] = MatrixChar_constructArray1D(len2, dim1, dim2);
        }
        return matArray;
}

void MatrixChar_destructArray1D(MatrixChar** matArray, int len){
        for(int i = 0; i < len; i++){
                MatrixChar_destruct(matArray[i]);
        }
        free(matArray);
}

void MatrixChar_destructArray2D(MatrixChar*** matArray, int len1, int len2){
        for(int i = 0; i < len1; i++){
                MatrixChar_destructArray1D(matArray[i], len2);
        }
        free(matArray);
}

void MatrixChar_setValue(MatrixChar* mat, uint8_t value){
        for(int i1 = 0; i1 < mat->dim1; i1++){
		for(int i2 = 0; i2 < mat->dim2; i2++){
                	mat->data[i1][i2] = value;
		}
        }
}

void MatrixChar_divideByValue(MatrixChar* mat, uint8_t value){
	assert(value != 0);
        for(int i1 = 0; i1 < mat->dim1; i1++){
                for(int i2 = 0; i2 < mat->dim2; i2++){
                        mat->data[i1][i2] /= value;
                }
        }
}

MatrixChar* MatrixChar_copy(MatrixChar* src){
        MatrixChar* dest = MatrixChar_construct1(src->dim1, src->dim2, src->data);
        return dest;
}

void MatrixChar_copyInPlace(MatrixChar* dest, MatrixChar* src){
	assert(dest->dim1 == src->dim1);
	assert(dest->dim2 == src->dim2);

        for(int i1 = 0; i1 < dest->dim1; i1++){
                for(int i2 = 0; i2 < dest->dim2; i2++){
                        dest->data[i1][i2] = src->data[i1][i2];
                }
        }
}

void MatrixChar_destruct(MatrixChar* mat){
	for(int i = 0; i < mat->dim1; i++){
		free(mat->data[i]);
		mat->data[i] = NULL;
	}
        free(mat->data);
        mat->data = NULL;
        free(mat);
}

char* MatrixChar_toString(MatrixChar* mat){
        char tmp[50];
        char* str = malloc(mat->dim1 * mat->dim2 * 10);
        str[0] = '[';
        int offset = 1;
        for(int i=0; i< mat->dim1; i++){
		for(int j=0; j< mat->dim2; j++){
                	sprintf(tmp, "%d", mat->data[i][j]);
                	memcpy(str + offset, tmp, strlen(tmp));
                	offset += strlen(tmp);
			if(j < mat->dim2 - 1){
				str[offset++] = ',';
				str[offset++] = ' ';
			}
		}
                if (i < mat->dim1 - 1){
                        str[offset++] = ',';
                        str[offset++] = '\n';
                        str[offset++] = ' ';
                }
        }
        str[offset++] = ']';
        str[offset] = '\0';
        return str;
}

// Functions for MatrixDouble
MatrixDouble* MatrixDouble_construct0(int dim1, int dim2){
        MatrixDouble* mat = (MatrixDouble*) malloc(sizeof(MatrixDouble));
	mat->dim1 = dim1;
	mat->dim2 = dim2;
        mat->data = (double**) malloc(dim1 * sizeof(double*));
        for(int i = 0; i < dim1; i++){
                mat->data[i] = (double*) malloc(dim2 * sizeof(double*));
                memset(mat->data[i], 0, dim2 * sizeof(double));
        }
        return mat;
}

MatrixDouble* MatrixDouble_construct1(int dim1, int dim2, double** data){
        MatrixDouble* mat = (MatrixDouble*) malloc(sizeof(MatrixDouble));
	mat->dim1 = dim1;
	mat->dim2 = dim2;
        mat->data = (double**) malloc(dim1 * sizeof(double*));
        for(int i = 0; i < dim1; i++){
                mat->data[i] = (double*) malloc(dim2 * sizeof(double*));
                memcpy(mat->data[i], data[i], dim2 * sizeof(double));
        }
        return mat;
}

MatrixDouble** MatrixDouble_constructArray1D(int len, int dim1, int dim2){
        MatrixDouble** matArray = (MatrixDouble**) malloc(len * sizeof(MatrixDouble*));
        for(int i = 0; i < len; i++){
                matArray[i] = MatrixDouble_construct0(dim1, dim2);
        }
        return matArray;
}


MatrixDouble*** MatrixDouble_constructArray2D(int len1, int len2, int dim1, int dim2){
        MatrixDouble*** matArray = (MatrixDouble***) malloc(len1 * sizeof(MatrixDouble**));
        for(int i = 0; i < len1; i++){
                matArray[i] = MatrixDouble_constructArray1D(len2, dim1, dim2);
        }
        return matArray;
}

void MatrixDouble_destructArray1D(MatrixDouble** matArray, int len){
        for(int i = 0; i < len; i++){
                MatrixDouble_destruct(matArray[i]);
        }
        free(matArray);
}

void MatrixDouble_destructArray2D(MatrixDouble*** matArray, int len1, int len2){
        for(int i = 0; i < len1; i++){
                MatrixDouble_destructArray1D(matArray[i], len2);
        }
        free(matArray);
}

void MatrixDouble_setValue(MatrixDouble* mat, double value){
        for(int i1 = 0; i1 < mat->dim1; i1++){
                for(int i2 = 0; i2 < mat->dim2; i2++){
                        mat->data[i1][i2] = value;
                }
        }
}

void MatrixDouble_divideByValue(MatrixDouble* mat, double value){
        assert(value != 0.0);
	for(int i1 = 0; i1 < mat->dim1; i1++){
                for(int i2 = 0; i2 < mat->dim2; i2++){
                        mat->data[i1][i2] /= value;
                }
        }
}

MatrixDouble* MatrixDouble_copy(MatrixDouble* src){
        MatrixDouble* dest = MatrixDouble_construct1(src->dim1, src->dim2, src->data);
        return dest;
}

void MatrixDouble_copyInPlace(MatrixDouble* dest, MatrixDouble* src){
        assert(dest->dim1 == src->dim1);
        assert(dest->dim2 == src->dim2);

        for(int i1 = 0; i1 < dest->dim1; i1++){
                for(int i2 = 0; i2 < dest->dim2; i2++){
                        dest->data[i1][i2] = src->data[i1][i2];
                }
        }
}

char* MatrixDouble_toString(MatrixDouble* mat){
        char tmp[50];
        char* str = malloc(mat->dim1 * mat->dim2 * 30);
        str[0] = '[';
        int offset = 1;
        for(int i=0; i< mat->dim1; i++){
                for(int j=0; j< mat->dim2; j++){
                        sprintf(tmp, "%.3e", mat->data[i][j]);
                        memcpy(str + offset, tmp, strlen(tmp));
                        offset += strlen(tmp);
                        if(j < mat->dim2 - 1){
                                str[offset++] = ',';
                                str[offset++] = ' ';
                        }
                }
                if (i < mat->dim1 - 1){
                        str[offset++] = ',';
                        str[offset++] = '\n';
                        str[offset++] = ' ';
                }
        }
        str[offset++] = ']';
        str[offset] = '\0';
        return str;
}

void MatrixDouble_destruct(MatrixDouble* mat){
        for(int i = 0; i < mat->dim1; i++){
                free(mat->data[i]);
                mat->data[i] = NULL;
        }
        free(mat->data);
        mat->data = NULL;
        free(mat);
}


MatrixDouble* MatrixDouble_parseFromFile(char* filePath, int dim1, int dim2){
	MatrixDouble* mat = MatrixDouble_construct0(dim1, dim2);
	FILE* f = fopen(filePath, "r");
	size_t len = 0;
    	char* line = NULL;
    	char* token;
    	ssize_t read; 
	int i = 0;
	int j = 0;
	Splitter* splitter;
	while ((read = getline(&line, &len, f)) != -1){
		line[strlen(line)-1] = '\0';
		splitter = Splitter_construct(line, '\t');
		j = 0;
    		while ((token = Splitter_getToken(splitter)) != NULL){
			mat->data[i][j] = atof(token);
			fprintf(stderr, "(%d,%d)%s\n",i,j,token);
			j += 1;
		}
		Splitter_destruct(splitter);
		i += 1;
	}
	fprintf(stderr, "done\n");
	fclose(f);
	return mat;
}
