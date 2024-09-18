#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "sonLib.h"
#include "common.h"
#include "track_reader.h"
#include <zlib.h>

# define LINE_MAX_SIZE 8192   /* line length maximum */

TrackFileFormat TrackReader_getTrackFileFormat(char *filePath) {
    char *extension = extractFileExtension(filePath);
    TrackFileFormat trackFileFormat;
    if (strcmp(extension, "cov") == 0)
        trackFileFormat = TRACK_FILE_FORMAT_COV;
    else if (strcmp(extension, "cov.gz") == 0)
        trackFileFormat = TRACK_FILE_FORMAT_COV_GZ;
    else if (strcmp(extension, "bed") == 0)
        trackFileFormat = TRACK_FILE_FORMAT_BED;
    else if (strcmp(extension, "bed.gz") == 0)
        trackFileFormat = TRACK_FILE_FORMAT_BED_GZ;
    else
        trackFileFormat = TRACK_FILE_FORMAT_UNDEFINED;

    if (trackFileFormat == TRACK_FILE_FORMAT_UNDEFINED) {
        fprintf(stderr, "[Error] %s should be either cov, cov.gz, bed or bed.gz! \n", filePath);
        exit(EXIT_FAILURE);
    }
    free(extension);
    return trackFileFormat;
}

stList *TrackReader_parseHeaderLines(TrackReader *trackReader) {
    stList *headerLines = stList_construct3(0, free);
    char *line = malloc(LINE_MAX_SIZE);
    ssize_t read;
    // set pointer to the start of the file
    TrackReader_setFilePosition(trackReader, 0);
    while (0 < (read = TrackReader_readLine(trackReader, &line, LINE_MAX_SIZE))) {
        if (line[0] == '#') {
            stList_append(headerLines, copyString(line));
        }
    }
    return headerLines;
}

CoverageHeader *CoverageHeader_construct(char *filePath) {
    CoverageHeader *header = malloc(sizeof(CoverageHeader));

    header->numberOfAnnotations = 0;
    header->annotationNames = stList_construct3(0, free);
    header->regionNames = stList_construct3(0,free);
    header->numberOfRegions = 0;
    header->regionCoverages = NULL;
    header->numberOfLabels = 0;
    header->isTruthAvailable = false;
    header->isPredictionAvailable = false;
    header->startOnlyMode = false;
    header->averageAlignmentLength = 0;

    if (filePath != NULL) {
        TrackReader *trackReader = TrackReader_construct(filePath, NULL, false);
        header->headerLines = TrackReader_parseHeaderLines(trackReader);

        CoverageHeader_updateAnnotationNames(header);
        CoverageHeader_updateRegionCoverages(header);
        CoverageHeader_updateRegionNames(header);
        CoverageHeader_updateTruthAvailability(header);
        CoverageHeader_updatePredictionAvailability(header);
        CoverageHeader_updateNumberOfLabels(header);
        CoverageHeader_updateStartOnlyMode(header); // first update mode then average length
        CoverageHeader_updateAverageAlignmentLength(header);

        TrackReader_destruct(trackReader);
    } else {
        header->headerLines = stList_construct3(0, free);
    }

    return header;
}

void CoverageHeader_updateRegionNames(CoverageHeader *header){
    if(header->regionNames != NULL){
        stList_destruct(header->regionNames);
    }
    header->regionNames = stList_construct3(0, free);
    char name[100];
    for(int i=0; i < header->numberOfRegions; i++){
        sprintf(name, "region_%d", i);
        stList_append(header->regionNames, copyString(name));
    }
}

void CoverageHeader_writeIntoFile(CoverageHeader *header,
                                  void *filePtr,
                                  bool isCompressed) {
    for (int i = 0; i < stList_length(header->headerLines); i++) {
        char *line = stList_get(header->headerLines, i);
        if (isCompressed) {
            gzFile *gzFilePtr = filePtr;
            gzprintf(*gzFilePtr, "%s\n", line);
        } else {
            FILE *fp = filePtr;
            fprintf(fp, "%s\n", line);
        }
    }
}

CoverageHeader *CoverageHeader_constructByAttributes(stList *annotationNames,
                                                     int *regionCoverages,
                                                     int numberOfRegions,
                                                     int numberOfLabels,
                                                     bool isTruthAvailable,
                                                     bool isPredictionAvailable,
                                                     bool startOnlyMode,
                                                     int averageAlignmentLength) {

    CoverageHeader *header = CoverageHeader_construct(NULL);

    char line[1000];

    header->numberOfAnnotations = stList_length(annotationNames);
    header->numberOfRegions = numberOfRegions;
    header->numberOfLabels = numberOfLabels;
    header->isTruthAvailable = isTruthAvailable;
    header->isPredictionAvailable = isPredictionAvailable;
    header->averageAlignmentLength = averageAlignmentLength;
    header->startOnlyMode = startOnlyMode;
    CoverageHeader_updateRegionNames(header);

    // add header line for annotation length
    sprintf(line, "#annotation:len:%d", stList_length(annotationNames));
    stList_append(header->headerLines, copyString(line));

    // update annotation attributes
    for (int i = 0; i < stList_length(annotationNames); i++) {
        // update annotation name
        char *annotationName = (char *) stList_get(annotationNames, i);
        stList_append(header->annotationNames, copyString(annotationName));
        // add header line for annotation name
        sprintf(line, "#annotation:name:%d:%s", i, annotationName);
        stList_append(header->headerLines, copyString(line));
    }

    header->regionCoverages = Int_construct1DArray(numberOfRegions);
    // add header line for number of regions
    sprintf(line, "#region:len:%d", numberOfRegions);
    stList_append(header->headerLines, copyString(line));

    // update region attributes
    for (int i = 0; i < numberOfRegions; i++) {
        // update region coverage
        header->regionCoverages[i] = regionCoverages[i];
        // add header line for region coverage
        sprintf(line, "#region:coverage:%d:%d", i, regionCoverages[i]);
        stList_append(header->headerLines, copyString(line));
    }

    // add number of labels for truth/prediction
    if (0 < numberOfLabels) {
        // add header line for number of labels
        sprintf(line, "#label:len:%d", numberOfLabels);
        stList_append(header->headerLines, copyString(line));
    }

    // are truth labels available
    if (isTruthAvailable) {
        sprintf(line, "#truth:true");
        stList_append(header->headerLines, copyString(line));
    } else {
        sprintf(line, "#truth:false");
        stList_append(header->headerLines, copyString(line));
    }

    // are prediction labels available
    if (isPredictionAvailable) {
        sprintf(line, "#prediction:true");
        stList_append(header->headerLines, copyString(line));
    } else {
        sprintf(line, "#prediction:false");
        stList_append(header->headerLines, copyString(line));
    }

    // add number of labels for truth/prediction
    if (0 < averageAlignmentLength) {
        // add header line for number of labels
        sprintf(line, "#avg_alignment_len:%d", averageAlignmentLength);
        stList_append(header->headerLines, copyString(line));
    }

    // start-only mode
    if (startOnlyMode) {
        sprintf(line, "#start-only:true");
        stList_append(header->headerLines, copyString(line));
    } else {
        sprintf(line, "#start-only:false");
        stList_append(header->headerLines, copyString(line));
    }

    return header;
}


void CoverageHeader_destruct(CoverageHeader *header) {
    if (header->headerLines != NULL) {
        stList_destruct(header->headerLines);
    }
    if (header->annotationNames != NULL) {
        stList_destruct(header->annotationNames);
    }
    if (header->regionNames != NULL) {
        stList_destruct(header->regionNames);
    }
    free(header->regionCoverages);
    free(header);
}


void CoverageHeader_updateNumberOfAnnotations(CoverageHeader *header) {
    header->numberOfAnnotations = 0;
    stList *headerLines = header->headerLines;
    char *token;
    bool lineFound = false;
    for (int i = 0; i < stList_length(headerLines); i++) {
        char *headerLine = stList_get(headerLines, i);
        if (strncmp("#annotation:len", headerLine, strlen("#annotation:len")) == 0) {
            Splitter *splitter = Splitter_construct(headerLine, ':');
            token = Splitter_getToken(splitter); //#annotation
            token = Splitter_getToken(splitter); //len
            token = Splitter_getToken(splitter); //number
            header->numberOfAnnotations = atoi(token);
            Splitter_destruct(splitter);
            lineFound = true;
            break;
        }
    }
    if (lineFound == false) {
        fprintf(stderr, "[%s] Error: No '#annotation:len:' found in the header. annotation len should be at least 1.\n",
                get_timestamp());
        exit(EXIT_FAILURE);
    }
    if (header->numberOfAnnotations <= 0 && lineFound) {
        fprintf(stderr, "[%s] Error: The value of '#annotation:len:' in the header should be at least 1.\n",
                get_timestamp());
        exit(EXIT_FAILURE);
    }
}

void CoverageHeader_updateAnnotationNames(CoverageHeader *header) {
    CoverageHeader_updateNumberOfAnnotations(header);
    if (header->annotationNames != NULL) {
        stList_destruct(header->annotationNames);
    }
    header->annotationNames = stList_construct3(header->numberOfAnnotations, free);
    for (int annotationIndex = 0; annotationIndex < header->numberOfAnnotations; annotationIndex++) {
        stList_set(header->annotationNames, annotationIndex, copyString("NA"));
    }

    stList *headerLines = header->headerLines;
    char *token;
    int numberOfParsedNames = 0;
    for (int i = 0; i < stList_length(headerLines); i++) {
        char *headerLine = stList_get(headerLines, i);
        if (strncmp("#annotation:name:", headerLine, strlen("#annotation:name:")) == 0) {
            Splitter *splitter = Splitter_construct(headerLine, ':');
            token = Splitter_getToken(splitter); //'#annotation'
            token = Splitter_getToken(splitter); //'name'
            token = Splitter_getToken(splitter); // index
            int annotationIndex = atoi(token);
            token = Splitter_getToken(splitter); //name
            stList_set(header->annotationNames, annotationIndex, copyString(token));
            numberOfParsedNames += 1;
            Splitter_destruct(splitter);
        }
    }
    if (numberOfParsedNames != header->numberOfAnnotations) {
        fprintf(stderr,
                "[%s] Warning: Number of parsed annotation names (%d) does not match '#annotation:len:%d' header line. Unnamed annotations will have the name 'NA'.\n",
                get_timestamp(), numberOfParsedNames, header->numberOfAnnotations);
    }
}

void CoverageHeader_updateAverageAlignmentLength(CoverageHeader *header) {
    header->averageAlignmentLength = 0;
    stList *headerLines = header->headerLines;
    char *token;
    bool lineFound = false;
    for (int i = 0; i < stList_length(headerLines); i++) {
        char *headerLine = stList_get(headerLines, i);
        if (strncmp("#avg_alignment_len:", headerLine, strlen("#avg_alignment_len:")) == 0) {
            Splitter *splitter = Splitter_construct(headerLine, ':');
            token = Splitter_getToken(splitter); //#avg_alignment_len:
            token = Splitter_getToken(splitter); //number
            header->averageAlignmentLength = atoi(token);
            Splitter_destruct(splitter);
            lineFound = true;
            break;
        }
    }
    if (lineFound == false && header->startOnlyMode) {
        fprintf(stderr, "[%s] Error: No '#avg_alignment_len:' found in the header (required for start-only mode)!\n",
                get_timestamp());
        exit(EXIT_FAILURE);
    }
    if (header->averageAlignmentLength <= 0 && lineFound && header->startOnlyMode) {
        fprintf(stderr, "[%s] Error: The value of '#avg_alignment_len:' in the header should be > 0 (required for start-only mode).\n",
                get_timestamp());
        exit(EXIT_FAILURE);
    }
}


void CoverageHeader_updateNumberOfRegions(CoverageHeader *header) {
    header->numberOfRegions = 0;
    stList *headerLines = header->headerLines;
    char *token;
    bool lineFound = false;
    for (int i = 0; i < stList_length(headerLines); i++) {
        char *headerLine = stList_get(headerLines, i);
        if (strncmp("#region:len", headerLine, strlen("#region:len")) == 0) {
            Splitter *splitter = Splitter_construct(headerLine, ':');
            token = Splitter_getToken(splitter); //#region
            token = Splitter_getToken(splitter); //len
            token = Splitter_getToken(splitter); //number
            header->numberOfRegions = atoi(token);
            Splitter_destruct(splitter);
            lineFound = true;
            break;
        }
    }
    if (lineFound == false) {
        fprintf(stderr, "[%s] Error: No '#region:len:' found in the header. region len should be at least 1.\n",
                get_timestamp());
        exit(EXIT_FAILURE);
    }
    if (header->numberOfRegions <= 0 && lineFound) {
        fprintf(stderr, "[%s] Error: The value of '#region:len:' in the header should be at least 1.\n",
                get_timestamp());
        exit(EXIT_FAILURE);
    }
}

void CoverageHeader_updateRegionCoverages(CoverageHeader *header) {
    CoverageHeader_updateNumberOfRegions(header);
    if (header->numberOfRegions == 0) {
        header->regionCoverages = NULL;
        return;
    }
    header->regionCoverages = (int *) malloc(header->numberOfRegions * sizeof(int));
    stList *headerLines = header->headerLines;
    char *token;
    int numberOfParsedCoverages = 0;
    for (int i = 0; i < stList_length(headerLines); i++) {
        char *headerLine = stList_get(headerLines, i);
        if (strncmp("#region:coverage:", headerLine, strlen("#region:coverage:")) == 0) {
            Splitter *splitter = Splitter_construct(headerLine, ':');
            token = Splitter_getToken(splitter); //'#region'
            token = Splitter_getToken(splitter); //'coverage'
            token = Splitter_getToken(splitter); // index
            int regionIndex = atoi(token);
            token = Splitter_getToken(splitter); // coverage
            header->regionCoverages[regionIndex] = atoi(token);
            numberOfParsedCoverages += 1;
            Splitter_destruct(splitter);
        }
    }
    if (numberOfParsedCoverages != header->numberOfRegions) {
        fprintf(stderr,
                "[%s] Error: Number of parsed region coverages (%d) does not match '#region:len:%d' in the header line.\n",
                get_timestamp(), numberOfParsedCoverages, header->numberOfRegions);
        exit(EXIT_FAILURE);
    }
}

// make sure to set truth and prediction availability before calling this function
void CoverageHeader_updateNumberOfLabels(CoverageHeader *header) {
    // set it to zero in case no label line existed in header
    header->numberOfLabels = 0;
    stList *headerLines = header->headerLines;
    char *token;
    bool lineFound = false;
    for (int i = 0; i < stList_length(headerLines); i++) {
        char *headerLine = stList_get(headerLines, i);
        if (strncmp("#label:len", headerLine, strlen("#label:len")) == 0) {
            Splitter *splitter = Splitter_construct(headerLine, ':');
            token = Splitter_getToken(splitter); //#label
            token = Splitter_getToken(splitter); //len
            token = Splitter_getToken(splitter); //number
            header->numberOfLabels = atoi(token);
            Splitter_destruct(splitter);
            lineFound = true;
            break;
        }
    }
    if (lineFound == false && (header->isTruthAvailable || header->isPredictionAvailable)) {
        fprintf(stderr,
                "[%s] Error: '#label:len' should be set to a non-zero number if at least one of truth or prediction tags is set to true in the header.\n",
                get_timestamp());
        exit(EXIT_FAILURE);
    }
}

void CoverageHeader_updateStartOnlyMode(CoverageHeader *header) {
    // set it to false in case no label line existed in header
    header->startOnlyMode = false;
    stList *headerLines = header->headerLines;
    for (int i = 0; i < stList_length(headerLines); i++) {
        char *headerLine = stList_get(headerLines, i);
        if (strncmp("#start-only:true", headerLine, strlen("#start-only:true")) == 0) {
            header->startOnlyMode = true;
            break;
        }
    }
    if (header->startOnlyMode == true) {
        fprintf(stderr,
                "[%s] Start-only mode is activated.\n",
                get_timestamp());
    }
}


void CoverageHeader_updateTruthAvailability(CoverageHeader *header) {
    // set it to false in case no label line existed in header
    header->isTruthAvailable = false;
    stList *headerLines = header->headerLines;
    for (int i = 0; i < stList_length(headerLines); i++) {
        char *headerLine = stList_get(headerLines, i);
        if (strncmp("#truth:true", headerLine, strlen("#truth:true")) == 0) {
            header->isTruthAvailable = true;
            break;
        }
    }
    if (header->isTruthAvailable == false) {
        fprintf(stderr,
                "[%s] Truth tag was set to false (or not defined) so truth labels will not be parsed from file.\n",
                get_timestamp());
    }
}

void CoverageHeader_updatePredictionAvailability(CoverageHeader *header) {
    // set it to false in case no label line existed in header
    header->isPredictionAvailable = false;
    stList *headerLines = header->headerLines;
    for (int i = 0; i < stList_length(headerLines); i++) {
        char *headerLine = stList_get(headerLines, i);
        if (strncmp("#prediction:true", headerLine, strlen("#prediction:true")) == 0) {
            header->isPredictionAvailable = true;
            break;
        }
    }
    if (header->isPredictionAvailable == false) {
        fprintf(stderr,
                "[%s] Prediction tag was set to false (or not defined) so prediction labels will not be parsed from file.\n",
                get_timestamp());
    }
}


void *TrackReader_openFile(char *filePath, TrackFileFormat format) {
    void *fileReaderPtr = NULL;
    if (format == TRACK_FILE_FORMAT_COV || format == TRACK_FILE_FORMAT_BED) {
        fileReaderPtr = fopen(filePath, "r");
        if (fileReaderPtr == NULL) {
            fprintf(stderr, "[Error] Unable to open %s\n", filePath);
            exit(EXIT_FAILURE);
        }
    } else if (format == TRACK_FILE_FORMAT_COV_GZ || format == TRACK_FILE_FORMAT_BED_GZ) {
        gzFile gzReader = gzopen(filePath, "r");
        if (gzReader == Z_NULL) {
            fprintf(stderr, "[Error] Unable to open %s\n", filePath);
            exit(EXIT_FAILURE);
        }
        gzFile *gzReaderPtr = malloc(sizeof(gzFile));
        gzReaderPtr[0] = gzReader;
        fileReaderPtr = gzReaderPtr;
    } else {
        fprintf(stderr, "ERROR: FORMAT should be either BED, COV, COV_GZ or BED_GZ!\n");
        exit(EXIT_FAILURE);
    }
    return fileReaderPtr;
}

int TrackReader_readLine(TrackReader *trackReader, char **linePtr, int maxSize) {
    size_t len = 0;
    ssize_t read = 0;
    // file can be either gz-compressed or not
    if (trackReader->trackFileFormat == TRACK_FILE_FORMAT_BED ||
        trackReader->trackFileFormat == TRACK_FILE_FORMAT_COV) {
        FILE *fileReaderPtr = (FILE *) trackReader->fileReaderPtr;
        read = getline(linePtr, &len, fileReaderPtr);
        char *line = *linePtr;
        if (0 < read && line[read - 1] == '\n') line[read - 1] = '\0'; // replace \n with \0
        // getline will return -1 if the file is ended
        return read;
    } else if (trackReader->trackFileFormat == TRACK_FILE_FORMAT_BED_GZ ||
               trackReader->trackFileFormat == TRACK_FILE_FORMAT_COV_GZ) {
        gzFile *fileReaderPtr = (gzFile *) trackReader->fileReaderPtr;
        gzFile gzReader = fileReaderPtr[0];
        char *line = gzgets(gzReader, *linePtr, maxSize);
        if (gzeof(gzReader)) {
            return -1;
        }
        if (line == Z_NULL) {
            int errnum;
            fprintf(stderr, "[Error] gz-comperssed file in TrackReader cannot be read properly: %s\n",
                    gzerror(gzReader, &errnum));
            exit(EXIT_FAILURE);
        }
        read = strlen(line);
        if (0 < read && line[read - 1] == '\n') line[read - 1] = '\0'; // replace \n with \0
        return read;
    } else {
        return -1;
    }
}

TrackReader *TrackReader_constructFromTableInMemory(stHash *coverageBlockTable, char *faiPath, bool zeroBasedCoors) {
    TrackReader *trackReader = malloc(sizeof(TrackReader));
    trackReader->trackFileFormat = TRACK_MEMORY_COV;
    trackReader->fileReaderPtr = NULL;
    if (faiPath != NULL) {
        trackReader->contigLengthTable = ptBlock_get_contig_length_stHash_from_fai(faiPath);
    } else {
        fprintf(stderr, "Error: fai path cannot be empty!");
        exit(EXIT_FAILURE);
    }
    trackReader->s = -1;
    trackReader->e = -1;
    trackReader->attrbs = NULL;
    trackReader->attrbsLen = 0;
    trackReader->zeroBasedCoors = zeroBasedCoors;
    trackReader->contigList = ptBlock_get_sorted_contig_list(blockTable);
    trackReader->coverageBlockTable = coverageBlockTable;
    trackReader->nextContigIndexToRead = 0;
    trackReader->nextBlockIndexToRead = 0;
    strcpy(trackReader->ctg, (char *)stList_get(trackReader->contigList, trackReader->nextContigIndexToRead));
    trackReader->coverageBlockListBeingIterated = (stList *) stHash_search(trackReader->coverageBlockTable, trackReader->ctg);
    int *ctgLenPtr = stHash_search(trackReader->contigLengthTable, trackReader->ctg);
    trackReader->ctgLen = *ctgLenPtr;
    return trackReader;
}

TrackReader *TrackReader_construct(char *filePath, char *faiPath, bool zeroBasedCoors) {
    TrackReader *trackReader = malloc(sizeof(TrackReader));
    trackReader->trackFileFormat = TrackReader_getTrackFileFormat(filePath);
    trackReader->fileReaderPtr = TrackReader_openFile(filePath, trackReader->trackFileFormat);
    if (faiPath != NULL) {
        trackReader->contigLengthTable = ptBlock_get_contig_length_stHash_from_fai(faiPath);
    } else {
        trackReader->contigLengthTable = NULL;
    }
    trackReader->ctgLen = -1;
    trackReader->s = -1;
    trackReader->e = -1;
    trackReader->attrbs = NULL;
    trackReader->attrbsLen = 0;
    trackReader->zeroBasedCoors = zeroBasedCoors;
    return trackReader;
}

int64_t TrackReader_getFilePosition(TrackReader *trackReader) {
    if (trackReader->trackFileFormat == TRACK_FILE_FORMAT_COV ||
        trackReader->trackFileFormat == TRACK_FILE_FORMAT_BED) {
        return ftell(trackReader->fileReaderPtr);
    } else if (trackReader->trackFileFormat == TRACK_FILE_FORMAT_COV_GZ ||
               trackReader->trackFileFormat == TRACK_FILE_FORMAT_BED_GZ) {
        gzFile *gzReaderPtr = trackReader->fileReaderPtr;
        gzFile gzReader = gzReaderPtr[0];
        return gztell(gzReader);
    } else {
        return -1;
    }
}

void TrackReader_setFilePosition(TrackReader *trackReader, int64_t filePosition) {
    if (trackReader->trackFileFormat == TRACK_FILE_FORMAT_COV ||
        trackReader->trackFileFormat == TRACK_FILE_FORMAT_BED) {
        fseek(trackReader->fileReaderPtr, filePosition, SEEK_SET);
    } else if (trackReader->trackFileFormat == TRACK_FILE_FORMAT_COV_GZ ||
               trackReader->trackFileFormat == TRACK_FILE_FORMAT_BED_GZ) {
        gzFile *gzReaderPtr = trackReader->fileReaderPtr;
        gzFile gzReader = gzReaderPtr[0];
        gzseek(gzReader, filePosition, SEEK_SET);
    }
}

void TrackReader_destruct(TrackReader *trackReader) {
    // free trackReader attrbs
    for (int i = 0; i < trackReader->attrbsLen; i++) {
        free(trackReader->attrbs[i]);
    }
    free(trackReader->attrbs);
    trackReader->attrbs = NULL;
    // close file
    if (trackReader->trackFileFormat == TRACK_FILE_FORMAT_COV ||
        trackReader->trackFileFormat == TRACK_FILE_FORMAT_BED) {
        fclose(trackReader->fileReaderPtr);
    } else if (trackReader->trackFileFormat == TRACK_FILE_FORMAT_COV_GZ ||
               trackReader->trackFileFormat == TRACK_FILE_FORMAT_BED_GZ) {
        gzFile *gzReaderPtr = trackReader->fileReaderPtr;
        gzclose(gzReaderPtr[0]);
        free(trackReader->fileReaderPtr);
    }
    if (trackReader->contigLengthTable != NULL) {
        stHash_destruct(trackReader->contigLengthTable);
    }
    if (trackReader->contigList != NULL){
        stList_destruct(trackReader->contigList);
    }
    free(trackReader);
}


int TrackReader_next(TrackReader *trackReader) {
    if (trackReader->trackFileFormat == TRACK_FILE_FORMAT_COV ||
        trackReader->trackFileFormat == TRACK_FILE_FORMAT_COV_GZ) {
        return TrackReader_readNextTrackCov(trackReader);
    } else if (trackReader->trackFileFormat == TRACK_FILE_FORMAT_BED ||
               trackReader->trackFileFormat == TRACK_FILE_FORMAT_BED_GZ) {
        return TrackReader_readNextTrackBed(trackReader);
    } else if (trackReader->trackFileFormat == TRACK_MEMORY_COV){
        return TrackReader_readNextFromMemory(trackReader);
    }
    else {
        fprintf(stderr, "ERROR: FORMAT should be either BED, COV, COV_GZ or BED_GZ or reading from memory!\n");
        exit(EXIT_FAILURE);
    }
}

int TrackReader_readNextFromMemory(TrackReader *trackReader){
    if( (trackReader->coverageBlockListBeingIterated != NULL) &&
        (stList_length(trackReader->coverageBlockListBeingIterated) <= trackReader->nextBlockIndexToRead)){
        // go to next contig
        trackReader->nextContigIndexToRead += 1;
        trackReader->nextBlockIndexToRead = 0;
        if (trackReader->nextContigIndexToRead < stList_length(trackReader->contigList)) {
            // update contig name
            strcpy(trackReader->ctg, (char *)stList_get(trackReader->contigList, trackReader->nextContigIndexToRead));
            // update coverage block list
            trackReader->coverageBlockListBeingIterated = (stList *) stHash_search(trackReader->coverageBlockTable, trackReader->ctg);
            // update contig length
            int *ctgLenPtr = stHash_search(trackReader->contigLengthTable, trackReader->ctg);
            trackReader->ctgLen = *ctgLenPtr;
        }else{ // no more contig is left
            trackReader->coverageBlockListBeingIterated = NULL;
            trackReader->nextBlockIndexToRead = -1;
            trackReader->nextContigIndexToRead = -1;
            trackReader->ctg[0] = '\0';
            trackReader->s = -1;
            trackReader->e = -1;
            trackReader->ctgLen = -1;
        }
    }
    if (trackReader->coverageBlockListBeingIterated == NULL) {
        return -1;
    } else {
        ptBlock *block = stList_get(trackReader->coverageBlockListBeingIterated, trackReader->nextBlockIndexToRead);
        trackReader->s = trackReader->zeroBasedCoors ? block->rfs : block->rfs + 1;
        trackReader->e = trackReader->zeroBasedCoors ? block->rfe : block->rfe + 1;
        CoverageInfo *coverageInfo = (CoverageInfo *) block->data;
        if (trackReader->attrbs == NULL) {
            trackReader->attrbs = malloc(5 * sizeof(char *));
            for (int i = 0; i < 5; i++) {
                trackReader->attrbs[i] = malloc(200 * sizeof(char));
            }
            trackReader->attrbsLen = 5;
        }
        int len = 0;
        int *annotation_indices = CoverageInfo_getAnnotationIndices(coverageInfo, &len);
        char *annotation_entry_str = String_joinIntArray(annotation_indices, len, ',');
        sprintf(trackReader->attrbs[0], "%.0f", coverageInfo->coverage);
        sprintf(trackReader->attrbs[1], "%.0f", coverageInfo->coverage_high_mapq);
        sprintf(trackReader->attrbs[2], "%.0f", coverageInfo->coverage_high_clip);
        sprintf(trackReader->attrbs[3], "%s", annotation_entry_str); // annotation string
        sprintf(trackReader->attrbs[4], "%d", CoverageInfo_getRegionIndex(coverageInfo)); // region index
        free(annotation_entry_str);
        free(annotation_indices);
        trackReader->nextBlockIndexToRead += 1;
        return 1; // just to be greater than 0
    }
}

int TrackReader_readNextTrackBed(TrackReader *trackReader) {
    char *line = malloc(LINE_MAX_SIZE);
    ssize_t read = TrackReader_readLine(trackReader, &line, LINE_MAX_SIZE);
    // if this is the end of the file
    if (read == -1) {
        trackReader->ctg[0] = '\0';
        trackReader->s = -1;
        trackReader->e = -1;
        free(line);
        return read;
    }
    char *token;

    if (read == 0) {
        fprintf(stderr, "[Warning] TrackReader was empty. Go to the next line!\n");
        return TrackReader_readNextTrackBed(trackReader);
    }


    if (line[0] == '#') { // skip header lines
        return TrackReader_readNextTrackBed(trackReader);
    }
    // free trackReader attrbs to fill it with new ones
    for (int i = 0; i < trackReader->attrbsLen; i++) {
        free(trackReader->attrbs[i]);
    }
    free(trackReader->attrbs);
    trackReader->attrbs = NULL;
    trackReader->attrbsLen = 0;
    trackReader->ctgLen = -1;

    Splitter *splitter = Splitter_construct(line, '\t');
    token = Splitter_getToken(splitter);
    strcpy(trackReader->ctg, token);

    // get contig length of fai was available
    if (trackReader->contigLengthTable != NULL) {
        int *ctg_len_ptr = stHash_search(trackReader->contigLengthTable, trackReader->ctg);
        if (ctg_len_ptr == NULL) {
            fprintf(stderr,
                    "[%s] Warning: fai file does not contain this contig name %s. The contig length will be set to -1.\n",
                    trackReader->ctg);
        } else {
            trackReader->ctgLen = *ctg_len_ptr;
        }
    }
    token = Splitter_getToken(splitter);
    trackReader->s = trackReader->zeroBasedCoors ? atoi(token) : atoi(token) + 1;
    token = Splitter_getToken(splitter);
    trackReader->e = trackReader->zeroBasedCoors ? atoi(token) - 1 : atoi(token);
    while ((token = Splitter_getToken(splitter)) != NULL) {
        trackReader->attrbsLen += 1;
        if (trackReader->attrbsLen == 1) {
            trackReader->attrbs = malloc(1 * sizeof(char *));
        } else { // increase the size of the attrbs if there are more attrbs
            trackReader->attrbs = realloc(trackReader->attrbs, trackReader->attrbsLen * sizeof(char *));
        }
        // save the currect attrb
        trackReader->attrbs[trackReader->attrbsLen - 1] = malloc(strlen(token) + 1);
        strcpy(trackReader->attrbs[trackReader->attrbsLen - 1], token);
    }
    Splitter_destruct(splitter);
    return read;
}

int TrackReader_readNextTrackCov(TrackReader *trackReader) {

    char *line = malloc(LINE_MAX_SIZE);

    ssize_t read = TrackReader_readLine(trackReader, &line, LINE_MAX_SIZE);
    // if this is the end of the file
    if (read == -1) {
        trackReader->ctg[0] = '\0';
        trackReader->ctgLen = 0;
        trackReader->s = -1;
        trackReader->e = -1;
        free(line);
        return read;
    }

    //fprintf(stderr,"%s\n", line);
    char *token;
    // free trackReader attrbs to fill it with new ones
    for (int i = 0; i < trackReader->attrbsLen; i++) {
        free(trackReader->attrbs[i]);
    }
    free(trackReader->attrbs);
    trackReader->attrbs = NULL;
    trackReader->attrbsLen = 0;

    if (read == 0) {
        fprintf(stderr, "[Warning] line read by TrackReader was empty. Go to the next line!\n");
        free(line);
        return TrackReader_readNextTrackCov(trackReader);
    }

    if (line[0] == '#') { // skip header lines
	free(line);
        return TrackReader_readNextTrackCov(trackReader);
    }

    if (line[0] == '>') { // contig name and size is after '>'
        Splitter *splitter = Splitter_construct(line, ' ');
        token = Splitter_getToken(splitter);
        strcpy(trackReader->ctg, token + 1); // skip '>' and copy
        token = Splitter_getToken(splitter);
        trackReader->ctgLen = atoi(token);
        Splitter_destruct(splitter);
	free(line);
        return TrackReader_readNextTrackCov(trackReader);
    } else {
        Splitter *splitter = Splitter_construct(line, '\t');
        token = Splitter_getToken(splitter);
        trackReader->s = trackReader->zeroBasedCoors ? atoi(token) - 1 : atoi(token);
        token = Splitter_getToken(splitter);
        trackReader->e = trackReader->zeroBasedCoors ? atoi(token) - 1 : atoi(token);
        while ((token = Splitter_getToken(splitter)) != NULL) {
            //fprintf(stderr, "%s\n",token);
            trackReader->attrbsLen += 1;
            if (trackReader->attrbsLen == 1) {
                trackReader->attrbs = malloc(1 * sizeof(char *));
            } else {// increase the size of the attrbs if there is more attrbs
                trackReader->attrbs = realloc(trackReader->attrbs, trackReader->attrbsLen * sizeof(char *));
            }
            // save the current attrb
            trackReader->attrbs[trackReader->attrbsLen - 1] = malloc(strlen(token) + 1);
            strcpy(trackReader->attrbs[trackReader->attrbsLen - 1], token);
        }
        Splitter_destruct(splitter);
	free(line);
	return read;
    }
}
