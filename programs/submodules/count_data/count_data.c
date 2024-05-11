
#include "common.h"
#include "count_data.h"




/////////////////////////////
//// CountData Functions ////
/////////////////////////////




CountData *CountData_construct(int length){
    CountData *countData = malloc(sizeof(CountData));
    countData->counts = Double_construct1DArray(length);
    countData->countsLength = length;
    countData->totalCount = 0;
    countData->sum = 0.0;
    return countData;
}

CountData **CountData_construct1DArray(int length, int arrayLength){
    CountData **countData1DArray = malloc(arrayLength * sizeof(CountData*));
    for(int i=0; i < arrayLength; i++){
        countData1DArray[i] = CountData_construct(length);
    }
    return countData1DArray;
}

void CountData_destruct1DArray(CountData **countData1DArray, int arrayLength){
    for(int i=0; i < arrayLength; i++){
        CountData_destruct(countData1DArray[i]);
    }
    free(countData1DArray);
}

void CountData_increment(CountData *countData, int value, double count) {
    uint8_t effValue = value;
    if (countData->countsLength <= value) {
        effValue = countData->countsLength - 1;
    }
    countData->counts[effValue] += count;
    countData->sum += effValue * count;
    countData->totalCount += count;
}

void CountData_reset(CountData* countData){
    memset(countData->counts, 0, countData->countsLength * sizeof(double));
    countData->sum = 0.0;
    countData->totalCount = 0;
}

int CountData_getMostFrequentValue(CountData *countData, int minValue, int maxValue){
    int value = -1;
    double maxCount = 0.0;
    for(int i=minValue; i < min(maxValue, countData->countsLength); i++){
        if(maxCount < countData->counts[i]){
            value = i;
            maxCount = countData->counts[i];
        }
    }
    return value;
}

double CountData_getMaxCount(CountData *countData, int minValue, int maxValue){
    int value = -1;
    double maxCount = 0.0;
    for(int i=minValue; i < min(maxValue, countData->countsLength); i++){
        if(maxCount < countData->counts[i]){
            value = i;
            maxCount = countData->counts[i];
        }
    }
    return maxCount;
}

double CountData_getTotalCount(CountData *countData, int minValue, int maxValue){
    double totalCount = 0.0;
    for(int i=minValue; i < min(maxValue, countData->countsLength); i++){
            totalCount += countData->counts[i];
    }
    return totalCount;
}

void CountData_destruct(CountData* countData){
    free(countData->counts);
    free(countData);
}


