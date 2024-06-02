#include <stdint.h>
#include "data_types.h"
#include "ptBlock.h"
#include "common.h"
#include "track_reader.h"
#include "chunk.h"

#define MAX_COVERAGE 2048

static const char *const LABEL_COLORS[] = {"162,0,37",
                                    "250,104,0",
                                    "0,138,0",
                                    "170,0,255",
                                    "99, 99, 96",
                                    "250,200,0"};
static const char *const LABEL_NAMES[] = {"Err",
                                   "Dup",
                                   "Hap",
                                   "Col",
                                   "Unk",
                                   "Msj"};




Chunk *Chunk_construct(int chunkCanonicalLen) {
    Chunk *chunk = malloc(sizeof(Chunk));
    chunk->chunkCanonicalLen = chunkCanonicalLen;
    chunk->coverageInfoSeq = NULL;
    chunk->coverageInfoSeqLen = 0;
    chunk->ctg[0] = '\0';
    chunk->ctgLen = 0;
    chunk->s = -1;
    chunk->e = -1;
    chunk->windowLen = -1;
    chunk->windowItr = -1;
    chunk->windowSumCoverage = 0.0;
    chunk->windowSumCoverageHighMapq = 0.0;
    chunk->windowSumCoverageHighClip = 0.0;
    chunk->windowRegionArray = NULL;
    chunk->windowTruthArray = NULL;
    chunk->fileOffset = 0;
    return chunk;
}

Chunk *Chunk_constructWithAllocatedSeq(int chunkCanonicalLen, int windowLen, int maxSeqSize) {
    Chunk *chunk = malloc(sizeof(Chunk));
    chunk->chunkCanonicalLen = chunkCanonicalLen;
    chunk->coverageInfoMaxSeqSize = maxSeqSize; //chunkCanonicalLen * 2 / windowLen + 1;
    chunk->coverageInfoSeq = (CoverageInfo **) CoverageInfo_construct1DArray(chunk->coverageInfoMaxSeqSize);
    chunk->coverageInfoSeqLen = 0;
    chunk->ctg[0] = '\0';
    chunk->ctgLen = 0;
    chunk->s = -1;
    chunk->e = -1;
    chunk->windowLen = windowLen;
    chunk->windowItr = -1;
    chunk->windowSumCoverage = 0.0;
    chunk->windowSumCoverageHighMapq = 0.0;
    chunk->windowSumCoverageHighClip = 0.0;
    chunk->windowAnnotationFlag = 0ULL;
    chunk->windowRegionArray = (int *) malloc(windowLen * sizeof(int));
    chunk->windowTruthArray = (int *) malloc(windowLen * sizeof(int));
    chunk->windowPredictionArray = (int *) malloc(windowLen * sizeof(int));
    chunk->fileOffset = 0;
    return chunk;
}

stList *Chunk_constructListWithAllocatedSeq(stList *templateChunks, int windowLen) {
    stList *chunks = stList_construct3(0, (void (*)(void *)) Chunk_destruct);
    for (int c = 0; c < stList_length(templateChunks); c++) {
        Chunk *templateChunk = stList_get(templateChunks, c);
        int maxSeqSize = (templateChunk->e - templateChunk->s + 1) / windowLen + 1;
        Chunk *chunk = Chunk_constructWithAllocatedSeq(templateChunk->chunkCanonicalLen, windowLen, maxSeqSize);
        stList_append(chunks, chunk);
    }
    return chunks;
}

int Chunk_getMaximumCoverageValue(Chunk *chunk) {
    int maxCoverage = 0.0;
    for (int i = 0; i < chunk->coverageInfoSeqLen; i++) {
        CoverageInfo *coverageInfo = chunk->coverageInfoSeq[i];
        if(maxCoverage < coverageInfo->coverage){
            maxCoverage = coverageInfo->coverage;
        }
    }
    return maxCoverage;
}

int Chunk_cmp(const void *chunk_1_, const void *chunk_2_) {
    const Chunk *chunk_1 = (Chunk *) chunk_1_;
    const Chunk *chunk_2 = (Chunk *) chunk_2_;
    if (strcmp(chunk_1->ctg, chunk_2->ctg) == 0) {
        return chunk_1->s - chunk_2->s;
    } else {
        return strcmp(chunk_1->ctg, chunk_2->ctg);
    }
}

void Chunk_destruct(Chunk *chunk) {
    if (chunk->coverageInfoSeq) {
        CoverageInfo_destruct1DArray(chunk->coverageInfoSeq, chunk->coverageInfoMaxSeqSize);
    }
    /*if (chunk->windowSumCoverageInfo) {
        CoverageInfo_destruct(chunk->windowSumCoverageInfo);
    }*/
    free(chunk->windowRegionArray);
    free(chunk->windowTruthArray);
    free(chunk);
}


ChunksCreator *ChunksCreator_constructEmpty() {
    ChunksCreator *chunksCreator = malloc(sizeof(ChunksCreator));
    chunksCreator->covPath = NULL;
    chunksCreator->header = CoverageHeader_construct(NULL);

    chunksCreator->nextChunkIndexToRead = 0;
    chunksCreator->nThreads = 1;
    chunksCreator->chunkCanonicalLen = 0;

    chunksCreator->chunks = NULL;
    chunksCreator->templateChunks = NULL;
    chunksCreator->windowLen = 0;
    chunksCreator->mutex = malloc(sizeof(pthread_mutex_t));
    pthread_mutex_init(chunksCreator->mutex, NULL);
    return chunksCreator;
}


ChunksCreator *
ChunksCreator_constructFromCov(char *covPath, char *faiPath, int chunkCanonicalLen, int nThreads, int windowLen) {
    char *extension = extractFileExtension(covPath);
    if (strcmp(extension, "cov") != 0 &&
        strcmp(extension, "cov.gz") != 0 &&
        strcmp(extension, "bed") != 0 &&
        strcmp(extension, "bed.gz") != 0) {
        fprintf(stderr,
                "[%s][Error] Coverage file's extension can be either '.cov', '.cov.gz', 'bed' or 'bed.gz': %s\n",
                get_timestamp(), covPath);
        exit(EXIT_FAILURE);
    }
    ChunksCreator *chunksCreator = malloc(sizeof(ChunksCreator));
    char covIndexPath[1000];
    sprintf(covIndexPath, "%s.index", covPath);
    if (!file_exists(covIndexPath)) {
        fprintf(stderr, "[%s] Index file does not exist: expected path %s\n", get_timestamp(), covIndexPath);
        fprintf(stderr, "[%s] Constructing index in memory ... \n", get_timestamp());
        chunksCreator->templateChunks = ChunksCreator_createCovIndex(covPath, faiPath, chunkCanonicalLen);
        fprintf(stderr, "[%s] Index is constructed.\n", get_timestamp());
        ChunksCreator_writeCovIndex(chunksCreator->templateChunks, covIndexPath);
        fprintf(stderr, "[%s] Index is saved into %s (It will skip index construction for next runs).\n", get_timestamp(), covIndexPath);

    } else {
        fprintf(stderr, "[%s] Index file exists: %s\n", get_timestamp(), covIndexPath);
        fprintf(stderr, "[%s] Parsing index ... \n", get_timestamp());
        chunksCreator->templateChunks = ChunksCreator_parseCovIndex(covIndexPath);
        fprintf(stderr, "[%s] Index is parsed from disk.\n", get_timestamp());

        // check canonical chunk length and reindex if necessary
        Chunk *firstChunk = stList_get(chunksCreator->templateChunks, 0);
        if (firstChunk->chunkCanonicalLen != chunkCanonicalLen) {
            fprintf(stderr,
                    "[%s] Warning: Chunk length of the parsed index is %d which does not match the expected value %d\n",
                    get_timestamp(), firstChunk->chunkCanonicalLen, chunkCanonicalLen);
            stList_destruct(chunksCreator->templateChunks);
            fprintf(stderr, "[%s] Constructing index in memory for chunk length of %d  ...\n", get_timestamp(),
                    chunkCanonicalLen);
            chunksCreator->templateChunks = ChunksCreator_createCovIndex(covPath, faiPath, chunkCanonicalLen);
            fprintf(stderr, "[%s] New index is constructed.\n", get_timestamp());
        }
    }
    chunksCreator->covPath = copyString(covPath);
    // parse attributes from header lines
    fprintf(stderr, "[%s] Parsing header info for ChunksCreator.\n", get_timestamp());
    chunksCreator->header = CoverageHeader_construct(covPath);

    chunksCreator->nextChunkIndexToRead = 0;
    chunksCreator->nThreads = nThreads;
    chunksCreator->chunkCanonicalLen = chunkCanonicalLen;
    fprintf(stderr, "[%s] Creating empty chunks.\n", get_timestamp());
    // create empty chunks
    chunksCreator->chunks = Chunk_constructListWithAllocatedSeq(chunksCreator->templateChunks, windowLen);
    chunksCreator->windowLen = windowLen;
    chunksCreator->mutex = malloc(sizeof(pthread_mutex_t));
    pthread_mutex_init(chunksCreator->mutex, NULL);
    free(extension);
    return chunksCreator;
}

int ChunksCreator_getMaximumCoverageValue(ChunksCreator *chunksCreator){
    int maxCoverage = 0.0;
    if (chunksCreator->chunks != NULL) {
        for (int chunkIndex = 0; chunkIndex < stList_length(chunksCreator->chunks); chunkIndex++) {
            Chunk *chunk = stList_get(chunksCreator->chunks, chunkIndex);
            int maxInChunk = Chunk_getMaximumCoverageValue(chunk);
            if (maxCoverage < maxInChunk){
                maxCoverage = maxInChunk;
            }
        }
    }
    return maxCoverage;
}

void ChunksCreator_subsetChunksToContigs(ChunksCreator *chunksCreator, stList* contigList){
    stList *newChunksList = stList_construct3(0, Chunk_destruct);
    if (chunksCreator->chunks != NULL) {
        for (int chunkIndex = 0; chunkIndex < stList_length(chunksCreator->chunks); chunkIndex++) {
            Chunk *chunk = stList_get(chunksCreator->chunks, chunkIndex);
            if (stList_existInStringList(contigList, chunk->ctg)) {
                // delete access from old list
                // to avoid issues in freeing memory
                stList_set(chunksCreator->chunks, chunkIndex, NULL);
                // add chunk to new list
                stList_append(newChunksList, chunk);
            }
        }
    }
    // delete old list
    stList_destruct(chunksCreator->chunks);
    // update chunks
    chunksCreator->chunks = newChunksList;
    chunksCreator->nextChunkIndexToRead = 0;
}

// it will create a stList of Chunks with no coverage data
stList *ChunksCreator_createCovIndex(char *filePath, char *faiPath, int chunkCanonicalLen) {
    bool zeroBasedCoors = true;
    TrackReader *trackReader = TrackReader_construct(filePath, faiPath, zeroBasedCoors);
    if (trackReader->trackFileFormat == TRACK_FILE_FORMAT_BED ||
        trackReader->trackFileFormat == TRACK_FILE_FORMAT_BED_GZ) {
        if (trackReader->contigLengthTable == NULL) {
            TrackReader_destruct(trackReader);
            fprintf(stderr,
                    "[%s] Error: for creating a index for .bed or .bed.gz file (%s) you need to pass a fai file.\n",
                    get_timestamp(), filePath);
            exit(EXIT_FAILURE);
        }
    }
    int64_t preFileOffset = TrackReader_getFilePosition(trackReader);
    Chunk *chunk = NULL;
    Chunk *preChunk = NULL;
    stList *chunks = stList_construct3(0, (void (*)(void *)) Chunk_destruct);
    while (0 < TrackReader_next(trackReader)) {
        if (chunk == NULL || strcmp(chunk->ctg, trackReader->ctg) != 0) { // contig has changed
            chunk = Chunk_construct(chunkCanonicalLen);
            chunk->s = 0;
            chunk->e = trackReader->ctgLen < 2 * chunkCanonicalLen ? trackReader->ctgLen - 1 : chunkCanonicalLen - 1;
            strcpy(chunk->ctg, trackReader->ctg);
            chunk->ctgLen = trackReader->ctgLen;
            chunk->fileOffset = preFileOffset;
            stList_append(chunks, chunk);
            //fprintf(stderr, "%d\t%s\t%d\t%d\t%d\t%ld\n", stList_length(chunks), chunk->ctg, chunk->ctgLen, chunk->s, chunk->e, chunk->fileOffset);

        }
            // has reached the end of the previous chunk and there is at least one more chunk to add
        else if (chunk->e <= trackReader->e && chunk->e < trackReader->ctgLen - 1) {
            preChunk = chunk;
            chunk = Chunk_construct(chunkCanonicalLen);
            chunk->s = preChunk->e + 1;
            chunk->e =
                    trackReader->ctgLen < preChunk->e + 2 * chunkCanonicalLen ? trackReader->ctgLen - 1 : preChunk->e +
                                                                                                          chunkCanonicalLen;
            strcpy(chunk->ctg, trackReader->ctg);
            chunk->ctgLen = trackReader->ctgLen;
            if (chunk->s <= trackReader->e)
                chunk->fileOffset = preFileOffset;
            else
                chunk->fileOffset = TrackReader_getFilePosition(trackReader);
            stList_append(chunks, chunk);
            //fprintf(stderr, "%d\t%s\t%d\t%d\t%d\t%ld\n", stList_length(chunks), chunk->ctg, chunk->ctgLen, chunk->s, chunk->e, chunk->fileOffset);
        }
        preFileOffset = TrackReader_getFilePosition(trackReader);
    }
    TrackReader_destruct(trackReader);
    return chunks; // The code should reach here when there is no more trackReader left in the cov file
}


void ChunksCreator_writeCovIndex(stList *templateChunks, char *indexPath) {
    FILE *filePtr = fopen(indexPath, "w+");
    if (filePtr == NULL) {
        printf("Couldn't open %s\n", indexPath);
        exit(EXIT_FAILURE);
    }
    assert(0 < stList_length(templateChunks));
    Chunk *templateChunk = stList_get(templateChunks, 0);
    fprintf(filePtr, "%d\n", templateChunk->chunkCanonicalLen);
    for (int i = 0; i < stList_length(templateChunks); i++) {
        templateChunk = stList_get(templateChunks, i);
        fprintf(filePtr, "%s\t%d\t%d\t%d\t%ld\n",
                templateChunk->ctg,
                templateChunk->ctgLen,
                templateChunk->s,
                templateChunk->e,
                templateChunk->fileOffset);
    }
    fflush(filePtr);
    fclose(filePtr);
    fprintf(stderr, "[%s] Index file (%s) is written \n", get_timestamp(), indexPath);
}

stList *ChunksCreator_parseCovIndex(char *covIndexPath) {
    FILE *filePtr = fopen(covIndexPath, "r");
    if (filePtr == NULL) {
        printf("[%s] Error: Couldn't open %s\n", get_timestamp(), covIndexPath);
        exit(EXIT_FAILURE);
    }
    size_t len = 0;
    char *line = NULL;
    char *token;
    ssize_t read;
    stList *templateChunks = stList_construct3(0, (void (*)(void *)) Chunk_destruct);
    read = getline(&line, &len, filePtr);
    int chunkCanonicalLen = atoi(line);
    while ((read = getline(&line, &len, filePtr)) != -1) {
        Chunk *templateChunk = Chunk_construct(chunkCanonicalLen);
        // contig name
        token = strtok(line, "\t");
        strcpy(templateChunk->ctg, token);
        // contig length
        token = strtok(NULL, "\t");
        templateChunk->ctgLen = atoi(token);
        // start
        token = strtok(NULL, "\t");
        templateChunk->s = atoi(token);
        // end
        token = strtok(NULL, "\t");
        templateChunk->e = atoi(token);
        // file offset
        token = strtok(NULL, "\t");
        templateChunk->fileOffset = atol(token);
        // add chunk to the list
        stList_append(templateChunks, templateChunk);
    }
    return templateChunks;
}

void ChunksCreator_destruct(ChunksCreator *chunksCreator) {
    if (chunksCreator->header != NULL) {
        CoverageHeader_destruct(chunksCreator->header);
    }
    if (chunksCreator->chunks != NULL) {
        stList_destruct(chunksCreator->chunks);
    }
    if (chunksCreator->templateChunks != NULL) {
        stList_destruct(chunksCreator->templateChunks);
    }
    free(chunksCreator->mutex);
    free(chunksCreator->covPath);
    free(chunksCreator);
}


void ChunksCreator_sortChunks(ChunksCreator *chunksCreator) {
    stList_sort(chunksCreator->chunks, Chunk_cmp);
    stList_sort(chunksCreator->templateChunks, Chunk_cmp);
}

int Chunk_getWindowTruth(Chunk *chunk) {
    // we will not have more than 10 states/labels for sure
    return Int_getModeValue1DArray(chunk->windowTruthArray, (chunk->windowItr + 1), -1, 10);
}

int Chunk_getWindowPrediction(Chunk *chunk) {
    // we will not have more than 10 states/labels for sure
    return Int_getModeValue1DArray(chunk->windowPredictionArray, (chunk->windowItr + 1), -1, 10);
}


int Chunk_getWindowRegion(Chunk *chunk) {
    // we will not have more than 100 regions for sure
    return Int_getModeValue1DArray(chunk->windowRegionArray, (chunk->windowItr + 1), 0, 100);
}

int Chunk_addWindow(Chunk *chunk) {
    if (chunk->windowItr == -1) return 1;
    double coverage_avg = (double) chunk->windowSumCoverage / (chunk->windowItr + 1);
    double coverage_high_mapq_avg = (double) chunk->windowSumCoverageHighMapq / (chunk->windowItr + 1);
    double coverage_high_clip_avg = (double) chunk->windowSumCoverageHighClip / (chunk->windowItr + 1);

    // set coverage values
    chunk->coverageInfoSeq[chunk->coverageInfoSeqLen]->coverage =
            MAX_COVERAGE < round(coverage_avg) ? MAX_COVERAGE : round(coverage_avg);
    chunk->coverageInfoSeq[chunk->coverageInfoSeqLen]->coverage_high_mapq =
            MAX_COVERAGE < round(coverage_high_mapq_avg) ? MAX_COVERAGE : round(coverage_high_mapq_avg);
    chunk->coverageInfoSeq[chunk->coverageInfoSeqLen]->coverage_high_clip =
            MAX_COVERAGE < round(coverage_high_clip_avg) ? MAX_COVERAGE : round(coverage_high_clip_avg);

    // set annotation bits
    chunk->coverageInfoSeq[chunk->coverageInfoSeqLen]->annotation_flag = chunk->windowAnnotationFlag;

    // get consensus values for truth/prediction labels and region index
    // it will take the values with the highest frequencies
    int8_t truth = Chunk_getWindowTruth(chunk);
    int8_t prediction = Chunk_getWindowPrediction(chunk);
    int regionIndex = Chunk_getWindowRegion(chunk);

    // set region bits
    CoverageInfo_setRegionIndex(chunk->coverageInfoSeq[chunk->coverageInfoSeqLen], regionIndex);

    // set inference data
    // truth and prediction will be -1 if not provided in the input cov file
    CoverageInfo_addInferenceData(chunk->coverageInfoSeq[chunk->coverageInfoSeqLen], truth, prediction);

    chunk->coverageInfoSeqLen += 1;
    chunk->windowItr = -1;
    chunk->windowAnnotationFlag = 0ULL;

    chunk->windowSumCoverage = 0.0;
    chunk->windowSumCoverageHighMapq = 0.0;
    chunk->windowSumCoverageHighClip = 0.0;
    return 0;
}


int Chunk_addTrack(Chunk *chunk, TrackReader *trackReader) {
    // check if the contig name matches
    assert(strcmp(chunk->ctg, trackReader->ctg) == 0);
    int canonicalPosInChunk = chunk->windowLen * chunk->coverageInfoSeqLen + chunk->windowItr + 1;
    // check if the overlap starts exactly one base after where the already added sequence ends
    //fprintf(stderr, "%d + %d * %d + %d + 1 == %d\n", chunk->s , chunk->windowLen, chunk->coverageInfoSeqLen, chunk->windowItr, max(trackReader->s, chunk->s));
    int canonicalStart = chunk->s + canonicalPosInChunk;
    assert(canonicalStart == max(trackReader->s, chunk->s));
    int canonicalBasesToAdd = min(trackReader->e, chunk->e) - max(trackReader->s, chunk->s) + 1;
    if (canonicalBasesToAdd <= 0) return 0;
    for (int i = 0; i < canonicalBasesToAdd; i++) {
        // windowItr initial value is -1
        chunk->windowItr += 1;
        chunk->windowItr %= chunk->windowLen;
        // the attrbs in the trackReader include coverage values and region index
        chunk->windowSumCoverage += atof(trackReader->attrbs[0]);
        chunk->windowSumCoverageHighMapq += atof(trackReader->attrbs[1]);
        chunk->windowSumCoverageHighClip += atof(trackReader->attrbs[2]);
        // get annotation indices and update window flag for annotation
        int len = 0;
        int *annotationIndices = Splitter_getIntArray(trackReader->attrbs[3], ',', &len);
        // use OR operation to keep all overlapping annotations
        chunk->windowAnnotationFlag |= CoverageInfo_getAnnotationFlagFromArray(annotationIndices, len);
        free(annotationIndices);
        // the first attribute right after annotation indices is the region index
        chunk->windowRegionArray[chunk->windowItr] = atoi(trackReader->attrbs[4]);
        // parse truth label if it exists (optional attribute)
        chunk->windowTruthArray[chunk->windowItr] = 6 <= trackReader->attrbsLen ? atoi(trackReader->attrbs[5]) : -1;
        // parse prediction label if it exists (optional attribute)
        chunk->windowPredictionArray[chunk->windowItr] =
                7 <= trackReader->attrbsLen ? atoi(trackReader->attrbs[6]) : -1;
        if (chunk->windowItr == chunk->windowLen - 1) { // the window is fully iterated
            int ret = Chunk_addWindow(chunk);
            if (ret == 1) {
                fprintf(stderr, "[Warning] Empty window cannot be added to chunk!\n");
            }
        }
    }
    return canonicalBasesToAdd;
}


int ChunksCreator_parseChunks(ChunksCreator *chunksCreator) {
    // create a thread pool
    tpool_t *tm = tpool_create(chunksCreator->nThreads);
    fprintf(stderr, "[%s] Created a thread pool with %d threads for parsing chunks \n", get_timestamp(),
            chunksCreator->nThreads);
    for (int i = 0; i < stList_length(chunksCreator->templateChunks); i++) {
        // Add a new job to the thread pool
        tpool_add_work(tm,
                       ChunksCreator_parseOneChunk,
                       (void *) chunksCreator);
        //ChunksCreator_parseOneChunk((void*) chunksCreator);
    }
    fprintf(stderr, "[%s] Queued %d jobs for the thread pool (no more than %d jobs will be processed at a time)\n",
            get_timestamp(), stList_length(chunksCreator->templateChunks), chunksCreator->nThreads);
    tpool_wait(tm);
    tpool_destroy(tm);

    return 0;
}

void ChunksCreator_parseOneChunk(void *chunksCreator_) {
    ChunksCreator *chunksCreator = (ChunksCreator *) chunksCreator_;
    // Lock the mutex to stop other threads from updating nextChunkIndexToRead and fetching a new chunk
    pthread_mutex_lock(chunksCreator->mutex);
    if (stList_length(chunksCreator->templateChunks) <= chunksCreator->nextChunkIndexToRead) {
        return;
    }
    Chunk *templateChunk = stList_get(chunksCreator->templateChunks, chunksCreator->nextChunkIndexToRead);
    // Construct a trackReader for iteration
    bool zeroBasedCoors = true;
    TrackReader *trackReader = TrackReader_construct(chunksCreator->covPath, NULL, zeroBasedCoors);
    strcpy(trackReader->ctg, templateChunk->ctg);
    trackReader->ctgLen = templateChunk->ctgLen;
    // Open cov file and jump to the first trackReader of the chunk
    TrackReader_setFilePosition(trackReader, templateChunk->fileOffset);
    // Get the chunk that has to be filled with the emitted sequence
    // and set start and end coordinates and also the contig name from the given template chunk
    Chunk *chunk = stList_get(chunksCreator->chunks, chunksCreator->nextChunkIndexToRead);
    chunksCreator->nextChunkIndexToRead += 1;
    // Unlock the mutex to let other threads fetch a chunk to read
    pthread_mutex_unlock(chunksCreator->mutex);
    chunk->coverageInfoSeqLen = 0;
    chunk->windowItr = -1;
    chunk->s = templateChunk->s;
    chunk->e = templateChunk->e;
    strcpy(chunk->ctg, templateChunk->ctg);
    chunk->ctgLen = templateChunk->ctgLen;
    // iterate over the tracks in the cov file
    while (0 < TrackReader_next(trackReader)) {
        // if the trackReader overlaps the chunk
        if (chunk->s <= trackReader->e && trackReader->s <= chunk->e) {
            assert(0 < Chunk_addTrack(chunk, trackReader));
        }
        if (chunk->e <= trackReader->e) { // chunk is read completely
            if (chunk->windowItr != -1) { // handle a partially iterated window
                Chunk_addWindow(chunk);
            }
            break;
        }
    }
    TrackReader_destruct(trackReader);
}


void ChunksCreator_writeChunksIntoBinaryFile(ChunksCreator *chunksCreator, char *binPath) {
    stList *chunks = chunksCreator->chunks;
    int chunkCanonicalLen = chunksCreator->chunkCanonicalLen;
    int windowLen = chunksCreator->windowLen;

    if (strcmp(extractFileExtension(binPath), "bin") != 0) {
        fprintf(stderr,
                "[%s] Warning: %s should end with '.bin' as the proper extension for saving chunks in binary format.\n",
                get_timestamp(), binPath);
    }

    FILE *fp = fopen(binPath, "wb+");

    int numberOfChunks = stList_length(chunks);
    if (chunksCreator->header == NULL) {
        fprintf(stderr, "[%s] Error: header cannot be undefined for writing chunks in binary.\n", get_timestamp());
        exit(EXIT_FAILURE);
    }
    CoverageHeader *header = chunksCreator->header;
    int numberOfAnnotations = header->numberOfAnnotations;
    stList *annotationNames = header->annotationNames;
    int *regionCoverages = header->regionCoverages;
    int numberOfRegions = header->numberOfRegions;
    int numberOfLabels = header->numberOfLabels;
    bool isTruthAvailable = header->isTruthAvailable;
    bool isPredictionAvailable = header->isPredictionAvailable;

    // write number of annotations
    fwrite(&numberOfAnnotations, sizeof(int32_t), 1, fp);
    // write annotation names
    for (int i = 0; i < numberOfAnnotations; i++) {
        char * annotationName = (char *) stList_get(annotationNames, i);
        int annotationNameSize = strlen(annotationName) + 1;
        // write number of chars in the string (+ null character)
        fwrite(&annotationNameSize, sizeof(int32_t), 1, fp);
        // write all chars ( + null character)
        fwrite(annotationName, sizeof(char), annotationNameSize, fp);
    }
    // write number of regions
    fwrite(&numberOfRegions, sizeof(int32_t), 1, fp);
    // write region coverages
    fwrite(regionCoverages, sizeof(int32_t), numberOfRegions, fp);
    // write number of labels (will be non-zero if there are truth/prediction labels in the input coverage file)
    fwrite(&numberOfLabels, sizeof(int32_t), 1, fp);
    // write truth availability
    fwrite(&isTruthAvailable, sizeof(bool), 1, fp);
    // write truth availability
    fwrite(&isPredictionAvailable, sizeof(bool), 1, fp);
    // write chunk length attributes
    fwrite(&chunkCanonicalLen, sizeof(int32_t), 1, fp);
    fwrite(&windowLen, sizeof(int32_t), 1, fp);

    // write chunks
    for (int c = 0; c < numberOfChunks; c++) {
        Chunk *chunk = stList_get(chunks, c);
        fprintf(stderr, "[%s] Writing the chunk %s:%d-%d\n", get_timestamp(), chunk->ctg, chunk->s, chunk->e);
        // write the length of contig name + null character
        int32_t ctgNameLen = strlen(chunk->ctg) + 1;
        fwrite(&ctgNameLen, sizeof(int32_t), 1, fp);
        // write the contig name
        fwrite(chunk->ctg, sizeof(char), ctgNameLen, fp);
        // write start and end
        fwrite(&(chunk->s), sizeof(int32_t), 1, fp);
        fwrite(&(chunk->e), sizeof(int32_t), 1, fp);
        //write actual size of chunk
        fwrite(&(chunk->coverageInfoSeqLen), sizeof(int32_t), 1, fp);
        uint16_t *covArray1 = malloc(chunk->coverageInfoSeqLen * sizeof(uint16_t));
        uint16_t *covArray2 = malloc(chunk->coverageInfoSeqLen * sizeof(uint16_t));
        uint16_t *covArray3 = malloc(chunk->coverageInfoSeqLen * sizeof(uint16_t));
        uint64_t *annotationArray = malloc(chunk->coverageInfoSeqLen * sizeof(uint64_t));
        int8_t *truthArray = malloc(chunk->coverageInfoSeqLen * sizeof(int8_t));
        int8_t *predictionArray = malloc(chunk->coverageInfoSeqLen * sizeof(int8_t));
        for (int i = 0; i < chunk->coverageInfoSeqLen; i++) {
            // put coverage and annotation data in related buffers
            covArray1[i] = (uint16_t) chunk->coverageInfoSeq[i]->coverage;
            covArray2[i] = (uint16_t) chunk->coverageInfoSeq[i]->coverage_high_mapq;
            covArray3[i] = (uint16_t) chunk->coverageInfoSeq[i]->coverage_high_clip;
            annotationArray[i] = (uint64_t) chunk->coverageInfoSeq[i]->annotation_flag;
            // put optional inference data in buffers
            Inference *inference = chunk->coverageInfoSeq[i]->data;
            if (inference != NULL) {
                truthArray[i] = (int8_t) inference->truth;
                predictionArray[i] = (int8_t) inference->prediction;
            } else {
                truthArray[i] = -1;
                predictionArray[i] = -1;
            }
        }
        // write arrays
        fwrite(covArray1, sizeof(uint16_t), chunk->coverageInfoSeqLen, fp);
        fwrite(covArray2, sizeof(uint16_t), chunk->coverageInfoSeqLen, fp);
        fwrite(covArray3, sizeof(uint16_t), chunk->coverageInfoSeqLen, fp);
        fwrite(annotationArray, sizeof(uint64_t), chunk->coverageInfoSeqLen, fp);
        fwrite(truthArray, sizeof(int8_t), chunk->coverageInfoSeqLen, fp);
        fwrite(predictionArray, sizeof(int8_t), chunk->coverageInfoSeqLen, fp);
        // free buffers
        free(covArray1);
        free(covArray2);
        free(covArray3);
        free(annotationArray);
        free(truthArray);
        free(predictionArray);
    }
    fflush(fp);
    fclose(fp);
    fprintf(stderr, "[%s] Done writing chunks in the binary file, %s\n", get_timestamp(), binPath);
}


// pass the output of ChunksCreator_constructEmpty() to this function
void ChunksCreator_parseChunksFromBinaryFile(ChunksCreator *chunksCreator, char *binPath) {
    if (!file_exists(binPath)) {
        fprintf(stderr,
                "[%s ] Error: The bin file %s does not exist. Please use create_bin_chunks for creating the bin file.\n",
                get_timestamp(), binPath);
        exit(EXIT_FAILURE);
    }
    FILE *fp = fopen(binPath, "rb");

    CoverageHeader *header = chunksCreator->header;

    // read number of annotations
    fread(&header->numberOfAnnotations, sizeof(int32_t), 1, fp);
    // read annotation names
    stList_destruct(header->annotationNames);
    header->annotationNames = stList_construct3(header->numberOfAnnotations, free);
    for (int annotationIndex = 0; annotationIndex < header->numberOfAnnotations; annotationIndex++) {
        // read size of annotation name
        int annotationNameSize = 0;
        fread(&annotationNameSize, sizeof(int32_t), 1, fp);
        // read annotation name
        char *annotationName = malloc(annotationNameSize * sizeof(char));
        fread(annotationName, sizeof(char), annotationNameSize, fp);
        // add annotation name to list
        stList_set(header->annotationNames, annotationIndex, annotationName);
    }
    // read number of regions
    fread(&header->numberOfRegions, sizeof(int32_t), 1, fp);
    // read region coverages
    free(header->regionCoverages);
    header->regionCoverages = Int_construct1DArray(header->numberOfRegions);
    fread(header->regionCoverages, sizeof(int32_t), header->numberOfRegions, fp);
    // read number of labels (will be non-zero if there are truth/prediction labels in the input coverage file)
    fread(&header->numberOfLabels, sizeof(int32_t), 1, fp);
    // read truth availability one bit
    fread(&header->isTruthAvailable, sizeof(bool), 1, fp);
    // read prediction availability one bit
    fread(&header->isPredictionAvailable, sizeof(bool), 1, fp);
    // update region names
    CoverageHeader_updateRegionNames(header);

    // read chunk length attributes
    fread(&chunksCreator->chunkCanonicalLen, sizeof(int32_t), 1, fp);
    fread(&chunksCreator->windowLen, sizeof(int32_t), 1, fp);


    chunksCreator->chunks = stList_construct3(0, Chunk_destruct);
    int32_t ctgNameLen;
    int maxSeqSize = chunksCreator->chunkCanonicalLen * 2 / chunksCreator->windowLen + 1;
    // create buffers for coverage values
    uint16_t *covArray1 = malloc(sizeof(uint16_t) * maxSeqSize);
    uint16_t *covArray2 = malloc(sizeof(uint16_t) * maxSeqSize);
    uint16_t *covArray3 = malloc(sizeof(uint16_t) * maxSeqSize);
    // create a buffer for annotation/region binary 64 bit flags
    uint64_t *annotationArray = malloc(sizeof(uint64_t) * maxSeqSize);
    // create buffers for truth/prediction labels
    int8_t *truthArray = malloc(sizeof(int8_t) * maxSeqSize);
    int8_t *predictionArray = malloc(sizeof(int8_t) * maxSeqSize);
    // create a buffer for contig name
    char ctg[1000];
    int32_t start;
    int32_t end;
    int32_t seqLen;
    // start reading binary chunks
    while (fread(&ctgNameLen, sizeof(int32_t), 1, fp) > 0) {
        // read the contig name
        fread(ctg, sizeof(char), ctgNameLen, fp);
        // read start and end
        fread(&start, sizeof(int32_t), 1, fp);
        fread(&end, sizeof(int32_t), 1, fp);
        // read actual size of chunk
        fread(&seqLen, sizeof(int32_t), 1, fp);
        //create chunk with enough allocated coverage sequence
        // here the size of the allocated array in chunk is equal to the actual size of the coverage sequence
        Chunk *chunk = Chunk_constructWithAllocatedSeq(chunksCreator->chunkCanonicalLen, chunksCreator->windowLen,
                                                       seqLen);
        // set contig name
        strcpy(chunk->ctg, ctg);
        // set coordinates
        chunk->s = start;
        chunk->e = end;
        chunk->coverageInfoSeqLen = seqLen;
        // read arrays
        fread(covArray1, sizeof(uint16_t), chunk->coverageInfoSeqLen, fp);
        fread(covArray2, sizeof(uint16_t), chunk->coverageInfoSeqLen, fp);
        fread(covArray3, sizeof(uint16_t), chunk->coverageInfoSeqLen, fp);
        fread(annotationArray, sizeof(uint64_t), chunk->coverageInfoSeqLen, fp);
        fread(truthArray, sizeof(int8_t), chunk->coverageInfoSeqLen, fp);
        fread(predictionArray, sizeof(int8_t), chunk->coverageInfoSeqLen, fp);
        for (int i = 0; i < chunk->coverageInfoSeqLen; i++) {
            chunk->coverageInfoSeq[i]->coverage = covArray1[i];
            chunk->coverageInfoSeq[i]->coverage_high_mapq = covArray2[i];
            chunk->coverageInfoSeq[i]->coverage_high_clip = covArray3[i];
            chunk->coverageInfoSeq[i]->annotation_flag = annotationArray[i];
            // add inference data
            CoverageInfo_addInferenceData(chunk->coverageInfoSeq[i], truthArray[i], predictionArray[i]);
        }
        stList_append(chunksCreator->chunks, chunk);
    }
    // free buffers
    free(covArray1);
    free(covArray2);
    free(covArray3);
    free(annotationArray);
    free(truthArray);
    free(predictionArray);

    fclose(fp);
}

void ChunksCreator_writeChunksIntoBedGraph(ChunksCreator *chunksCreator,
                                           const char *outputPath,
                                           const char *trackName,
                                           u_int16_t (*getCoverageInfoAttribute)(CoverageInfo *),
                                           const char *color) {
    FILE *fp = fopen(outputPath, "w");
    if (fp == NULL) {
        fprintf(stderr, "[%s] Error: Failed to open file %s.\n", get_timestamp(), outputPath);
    }
    // since chunkCanonicalLen is set to maximum contig size
    // then each contig is parsed in a single chunk
    for (int chunkIndex = 0; chunkIndex < stList_length(chunksCreator->chunks); chunkIndex++) {
        Chunk *chunk = stList_get(chunksCreator->chunks, chunkIndex);
        fprintf(fp, "track type=bedGraph name=\"%s\" visibility=full color=%s\n", trackName, color);
        for (int i = 0; i < chunk->coverageInfoSeqLen; i++) {
            CoverageInfo *coverageInfo = chunk->coverageInfoSeq[i];
            int start = chunk->s + i * chunk->windowLen;
            int end = min(chunk->s + (i + 1) * chunk->windowLen, chunk->e + 1);
            fprintf(fp, "%s %d %d %d\n", chunk->ctg, start, end, getCoverageInfoAttribute(coverageInfo));
        }
        fprintf(fp, "\n");

    }
    fclose(fp);
}


ChunkIterator *ChunkIterator_construct(ChunksCreator *chunksCreator) {
    ChunkIterator *chunkIterator = malloc(sizeof(ChunkIterator));
    chunkIterator->numberOfChunks = stList_length(chunksCreator->chunks);
    chunkIterator->chunksCreator = chunksCreator;
    chunkIterator->nextChunkIndex = 0;
    chunkIterator->nextWindowIndex = 0;
    chunkIterator->block = ptBlock_construct(-1, -1, -1, -1, -1, -1);
    // set coverage info data initially to NULL
    // set the related functions
    ptBlock_set_data(chunkIterator->block,
                     NULL,
                     destruct_cov_info_data,
                     copy_cov_info_data,
                     extend_cov_info_data);
    return chunkIterator;
}

ChunkIterator *ChunkIterator_copy(ChunkIterator *src) {
    ChunkIterator *dest = malloc(sizeof(ChunkIterator));
    dest->numberOfChunks = src->numberOfChunks;
    dest->chunksCreator = src->chunksCreator;
    dest->nextChunkIndex = src->nextChunkIndex;
    dest->nextWindowIndex = src->nextWindowIndex;
    dest->block = ptBlock_construct(src->block->rfs, src->block->rfe, -1, -1, -1, -1);
    // set coverage info data initially to NULL
    // set the related functions
    ptBlock_set_data(dest->block,
                     src->block->data,
                     destruct_cov_info_data,
                     copy_cov_info_data,
                     extend_cov_info_data);
    return dest;
}

void ChunkIterator_reset(ChunkIterator *chunkIterator) {
    chunkIterator->nextChunkIndex = 0;
    chunkIterator->nextWindowIndex = 0;
    if(chunkIterator->block != NULL){
	    chunkIterator->block->data = NULL;
	    ptBlock_destruct(chunkIterator->block);
    }
    chunkIterator->block = ptBlock_construct(-1, -1, -1, -1, -1, -1);
    // set coverage info data initially to NULL
    // set the related functions
    ptBlock_set_data(chunkIterator->block,
                     NULL,
                     destruct_cov_info_data,
                     copy_cov_info_data,
                     extend_cov_info_data);
}


void ChunkIterator_destruct(ChunkIterator *chunkIterator) {
    if (chunkIterator->block != NULL) {
        chunkIterator->block->data = NULL; // this data is coming from another source (not copied)
        ptBlock_destruct(chunkIterator->block);
    }
    free(chunkIterator);
}

ptBlock *ChunkIterator_getNextPtBlock(ChunkIterator *chunkIterator, char *ctg_name) {
    if (chunkIterator->nextChunkIndex == chunkIterator->numberOfChunks) {
        ctg_name[0] = '\0';
	if (chunkIterator->block != NULL){
		chunkIterator->block->data = NULL;
		ptBlock_destruct(chunkIterator->block);
		chunkIterator->block = NULL;
	}
        return chunkIterator->block;
    }
    Chunk *chunk = stList_get(chunkIterator->chunksCreator->chunks, chunkIterator->nextChunkIndex);
    // this chunk is finished
    if (chunkIterator->nextWindowIndex == chunk->coverageInfoSeqLen) {
        chunkIterator->nextChunkIndex += 1;
        chunkIterator->nextWindowIndex = 0;
        return ChunkIterator_getNextPtBlock(chunkIterator, ctg_name);
    }

    int i = chunkIterator->nextWindowIndex;
    int start = chunk->s + i * chunk->windowLen; //0-based inclusive
    int end = min(chunk->s + (i + 1) * chunk->windowLen - 1, chunk->e); //0-based inclusive
    CoverageInfo *coverageInfo = chunk->coverageInfoSeq[chunkIterator->nextWindowIndex];

    // update block attributes based on the current window
    chunkIterator->block->data = coverageInfo;
    chunkIterator->block->rfs = start;
    chunkIterator->block->rfe = end;
    // update ctg name
    strcpy(ctg_name, chunk->ctg);

    // go to the next window, if it is out of bound
    // it will be handled in the next call
    chunkIterator->nextWindowIndex += 1;

    // return update block
    return chunkIterator->block;
}


void ChunksCreator_writePredictionIntoFinalBED(ChunksCreator *chunksCreator, char *outputPath, char *trackName) {

    // open file for writing bed
    FILE *fout = fopen(outputPath, "w");
    if (fout == NULL) {
        fprintf(stderr, "[%s] Error: %s cannot be opened.\n", get_timestamp(), outputPath);
        exit(EXIT_FAILURE);
    }
    // write first line of BED file
    fprintf(fout, "track name=%s visibility=1 itemRgb=\"On\"\n", trackName);

    ChunkIterator *iterator = ChunkIterator_construct(chunksCreator);

    ptBlock *block = NULL;
    int bedTrackStart = 0;
    int preEnd = 0;
    int preLabel = -1;
    char ctg[200];
    char preCtg[200];
    preCtg[0] = '\0';

    while ((block = ChunkIterator_getNextPtBlock(iterator, ctg)) != NULL) {

        // get coverage info for this block
        CoverageInfo *coverageInfo = (CoverageInfo *) block->data;
        int start = block->rfs;
        int end = block->rfe;

        // initialize start
        if(preLabel == -1 || preCtg[0] == '\0'){
            bedTrackStart = start;
        }

        // index 4 is for "Unk"
        int predictionLabel = 4;
        // get labels
        if (coverageInfo->data == NULL) {
            fprintf(stderr, "[%s] Warning: inference data %s:%d-%d does not exist for writing final bed.\n",
                    get_timestamp(),
                    start,
                    end,
                    ctg);
        }else {
            Inference *inference = coverageInfo->data;
            if (inference->prediction != -1) {
                predictionLabel = inference->prediction;
            }
        }

        bool labelChanged = preLabel != -1 && predictionLabel != preLabel;
        bool contigChanged = preCtg[0] != '\0' && strcmp(preCtg, ctg) != 0;
        if(labelChanged || contigChanged){
            fprintf(fout, "%s\t%d\t%d\t%s\t0\t+\t%d\t%d\t%s\n",
                    preCtg,
                    bedTrackStart,
                    preEnd + 1,
                    LABEL_NAMES[preLabel],
                    bedTrackStart,
                    preEnd + 1,
                    LABEL_COLORS[preLabel]);
            // update start
            bedTrackStart = start;
        }

        preEnd = end;
        preLabel = predictionLabel;
        strcpy(preCtg, ctg);
    }
    if (preLabel != -1) {
        fprintf(fout, "%s\t%d\t%d\t%s\t0\t+\t%d\t%d\t%s\n",
                preCtg,
                bedTrackStart,
                preEnd + 1,
                LABEL_NAMES[preLabel],
                bedTrackStart,
                preEnd + 1,
                LABEL_COLORS[preLabel]);
    }

    // close file
    fclose(fout);

    ChunkIterator_destruct(iterator);
}

int ChunksCreator_getTotalNumberOfChunks(ChunksCreator *chunksCreator) {
    return (int) stList_length(chunksCreator->chunks);
}

int64_t ChunksCreator_getTotalLength(ChunksCreator *chunksCreator) {
    ChunkIterator *iterator = ChunkIterator_construct(chunksCreator);

    ptBlock *block = NULL;

    char ctg[200];
    int64_t totalLength = 0;
    while ((block = ChunkIterator_getNextPtBlock(iterator, ctg)) != NULL) {
        int start = block->rfs;
        int end = block->rfe;
        totalLength += end - start + 1;
    }

    ChunkIterator_destruct(iterator);
    return totalLength;
}

