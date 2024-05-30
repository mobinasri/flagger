#ifndef COUNT_DATA_H
#define COUNT_DATA_H

#include <string.h>
#include <stdint.h>
#include <assert.h>
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "stdbool.h"
#include "common.h"


/*! @typedef
 * @abstract Structure for storing counts of integer values ranging from 0 up to a specified countsLength.
 * @field counts        An array containing the counts for each integer value. Each index represents a value, and
 *                      each count can be a floating-point number.
 * @field countsLength  The length of the counts array.
 * @field totalCount    The total number of values added to this structure.
 * @field sum           The sum of all values calculated by multiplying each index by its associated count.
 *                      In other words, it represents the sum of all values added to this structure.
 */
typedef struct CountData{
    double *counts;
    int countsLength;
    int totalCount;
    double sum;
}CountData;

/*
 * Construct a CountData structure
 */
CountData *CountData_construct(int length);

CountData *CountData_copy(CountData *src);

/*
 * Construct a 1D array of CountData structures
 */
CountData **CountData_construct1DArray(int length, int arrayLength);

CountData **CountData_copy1DArray(CountData **srcArray, int arrayLength);

/*
 * Increment the count of a specific value
 */
void CountData_increment(CountData * countData, int value, double count);

/*
 * Set all counts and related attributes to zero
 */
void CountData_reset(CountData* countData);

/*
 * Get the most frequent value (index) from the given inclusive interval
 */
int CountData_getMostFrequentValue(CountData * countData, int minValue, int maxValue);

/*
 * Get the count of the most frequent value (index) from the given inclusive interval
 */
double CountData_getMaxCount(CountData *countData, int minValue, int maxValue);

/*
 * Get total count from the given inclusive interval
 */
double CountData_getTotalCount(CountData *countData, int minValue, int maxValue);

/*
 * Destruct a CountData structure
 */
void CountData_destruct(CountData* countData);

/*
 * Destruct a 1D array of CountData structures
 */
void CountData_destruct1DArray(CountData **countData1DArray, int arrayLength);


#endif // COUNT_DATA_H
