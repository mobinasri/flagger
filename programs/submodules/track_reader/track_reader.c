#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "sonLib.h"
#include "common.h"
#include "track_reader.h"
#include <zlib.h>

# define LINE_MAX_SIZE 8192   /* line length maximum */

TrackFileFormat TrackReader_getTrackFileFormat(char *filePath){
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

	if (trackFileFormat == TRACK_FILE_FORMAT_UNDEFINED){
		fprintf(stderr, "[Error] %s should be either cov, cov.gz, bed or bed.gz! \n", filePath);
		exit(EXIT_FAILURE);
	}
	free(extension);
	return trackFileFormat;
}

void *TrackReader_openFile(char *filePath, TrackFileFormat format){
	void *fileReaderPtr = NULL;
	if (format == TRACK_FILE_FORMAT_COV || format == TRACK_FILE_FORMAT_BED){
                fileReaderPtr = fopen(filePath, "r");
                if (fileReaderPtr == NULL){
                        fprintf(stderr, "[Error] Unable to open %s\n", filePath);
                        exit(EXIT_FAILURE);
                }
        } else if (format == TRACK_FILE_FORMAT_COV_GZ || format == TRACK_FILE_FORMAT_BED_GZ){
                gzFile gzReader = gzopen(filePath, "r");
                if (gzReader == Z_NULL){
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

int TrackReader_readLine(TrackReader *trackReader, char **linePtr, int maxSize){
	size_t len = 0;
	ssize_t read = 0;
	// file can be either gz-compressed or not
	if (trackReader->trackFileFormat == TRACK_FILE_FORMAT_BED || trackReader->trackFileFormat == TRACK_FILE_FORMAT_COV){
		FILE *fileReaderPtr = (FILE *) trackReader->fileReaderPtr;
		read = getline(linePtr, &len, fileReaderPtr);
		char *line = *linePtr;
		if (0 < read && line[read-1] == '\n') line[read-1] = '\0'; // replace \n with \0
		// getline will return -1 if the file is ended
		return read;
	}
	else if (trackReader->trackFileFormat == TRACK_FILE_FORMAT_BED_GZ || trackReader->trackFileFormat == TRACK_FILE_FORMAT_COV_GZ){
		gzFile *fileReaderPtr = (gzFile *) trackReader->fileReaderPtr;
		gzFile gzReader = fileReaderPtr[0];
		char *line = gzgets(gzReader, *linePtr, maxSize);
		if (gzeof(gzReader)){
			return -1;
		}
		if(line == Z_NULL){
			int errnum;
			fprintf(stderr, "[Error] gz-comperssed file in TrackReader cannot be read properly: %s\n", gzerror(gzReader, &errnum));
			exit(EXIT_FAILURE);
		}
		read = strlen(line);
		if (0 < read && line[read-1] == '\n') line[read-1] = '\0'; // replace \n with \0
		return read;
	}else{
		return -1;
	}
}

void TrackReader_updateHeaderLines(TrackReader *trackReader){
	char *line = malloc(LINE_MAX_SIZE);
	ssize_t read;
	// set pointer to the start of the file
	TrackReader_setFilePosition(trackReader, 0);
	int64_t filePosition = TrackReader_getFilePosition(trackReader);
	while(0 < (read = TrackReader_readLine(trackReader, &line, LINE_MAX_SIZE))){
		if (line[0] == '#'){
			stList_append(trackReader->headerLines, copyString(line));
			// update file position
			filePosition = TrackReader_getFilePosition(trackReader);
		}
		else{
			// set file position to the beginning of this line
			// since this line has information about the first contig
			TrackReader_setFilePosition(trackReader, filePosition);
			break;
		}
	}
}

TrackReader* TrackReader_construct(char *filePath, char *faiPath, bool zeroBasedCoors){
	TrackReader *trackReader = malloc(sizeof(TrackReader));
	trackReader->trackFileFormat = TrackReader_getTrackFileFormat(filePath);
	trackReader->fileReaderPtr = TrackReader_openFile(filePath, trackReader->trackFileFormat);
	trackReader->headerLines = stList_construct3(0,free);
	TrackReader_updateHeaderLines(trackReader);
	if(faiPath != NULL){
		trackReader->contigLengthTable = ptBlock_get_contig_length_stHash_from_fai(faiPath);
	}else{
		trackReader->contigLengthTable = NULL;
	}
	trackReader->ctgLen = -1;
	trackReader->s=-1; 
	trackReader->e=-1;
	trackReader->attrbs = NULL;
	trackReader->attrbsLen = 0;
	trackReader->zeroBasedCoors = zeroBasedCoors;
	return trackReader;
}

int64_t TrackReader_getFilePosition(TrackReader *trackReader){
	if (trackReader->trackFileFormat == TRACK_FILE_FORMAT_COV || trackReader->trackFileFormat == TRACK_FILE_FORMAT_BED){
                return ftell(trackReader->fileReaderPtr);
        } else if (trackReader->trackFileFormat == TRACK_FILE_FORMAT_COV_GZ || trackReader->trackFileFormat == TRACK_FILE_FORMAT_BED_GZ){
                gzFile *gzReaderPtr = trackReader->fileReaderPtr;
		gzFile gzReader = gzReaderPtr[0];
                return gztell(gzReader);
        } else{
		return -1;
	}
}

void TrackReader_setFilePosition(TrackReader *trackReader, int64_t filePosition){
	if (trackReader->trackFileFormat == TRACK_FILE_FORMAT_COV || trackReader->trackFileFormat == TRACK_FILE_FORMAT_BED){
		fseek(trackReader->fileReaderPtr, filePosition, SEEK_SET);
        } else if (trackReader->trackFileFormat == TRACK_FILE_FORMAT_COV_GZ || trackReader->trackFileFormat == TRACK_FILE_FORMAT_BED_GZ){
                gzFile *gzReaderPtr = trackReader->fileReaderPtr;
                gzFile gzReader = gzReaderPtr[0];
		gzseek(gzReader, filePosition, SEEK_SET);
        }
}


void TrackReader_destruct(TrackReader *trackReader){
	// free trackReader attrbs
    	for(int i = 0; i < trackReader->attrbsLen; i++){
                 free(trackReader->attrbs[i]);
	}
    	free(trackReader->attrbs);
    	trackReader->attrbs = NULL;
	// close file
	if (trackReader->trackFileFormat == TRACK_FILE_FORMAT_COV || trackReader->trackFileFormat == TRACK_FILE_FORMAT_BED){
                fclose(trackReader->fileReaderPtr);
        } else if (trackReader->trackFileFormat == TRACK_FILE_FORMAT_COV_GZ || trackReader->trackFileFormat == TRACK_FILE_FORMAT_BED_GZ){
		gzFile *gzReaderPtr = trackReader->fileReaderPtr;
                gzclose(gzReaderPtr[0]);
		free(trackReader->fileReaderPtr);
        }
	if (trackReader->contigLengthTable != NULL){
		stHash_destruct(trackReader->contigLengthTable);
	}
	if (trackReader->headerLines != NULL){
                stList_destruct(trackReader->headerLines);
        }
	free(trackReader);
}


int TrackReader_next(TrackReader* trackReader){
	if (trackReader->trackFileFormat == TRACK_FILE_FORMAT_COV || trackReader->trackFileFormat == TRACK_FILE_FORMAT_COV_GZ){
		return TrackReader_readNextTrackCov(trackReader);
	}
	else if(trackReader->trackFileFormat == TRACK_FILE_FORMAT_BED || trackReader->trackFileFormat == TRACK_FILE_FORMAT_BED_GZ){
		return TrackReader_readNextTrackBed(trackReader);
	}
	else{
		fprintf(stderr, "ERROR: FORMAT should be either BED, COV, COV_GZ or BED_GZ!\n");
		exit(EXIT_FAILURE);
	}
}

int TrackReader_readNextTrackBed(TrackReader* trackReader){
    char *line = malloc(LINE_MAX_SIZE);
    ssize_t read = TrackReader_readLine(trackReader, &line, LINE_MAX_SIZE);
    // if this is the end of the file
    if (read == -1){
         trackReader->ctg[0] = '\0';
         trackReader->s = -1;
         trackReader->e = -1;
	 free(line);
         return read;
    }
    char* token;
    if(read == 0){
            fprintf(stderr, "[Warning] TrackReader was empty. Go to the next line!\n");
            return TrackReader_readNextTrackBed(trackReader);
    }
    // free trackReader attrbs to fill it with new ones
    for(int i = 0; i < trackReader->attrbsLen; i++){
                 free(trackReader->attrbs[i]);
    }
    free(trackReader->attrbs);
    trackReader->attrbs = NULL;
    trackReader->attrbsLen = 0;
    trackReader->ctgLen = -1;

    Splitter* splitter = Splitter_construct(line, '\t');
    token = Splitter_getToken(splitter);
    strcpy(trackReader->ctg, token);

    // get contig length of fai was available
    if(trackReader->contigLengthTable != NULL){
	    int *ctg_len_ptr = stHash_search(trackReader->contigLengthTable, trackReader->ctg);
	    if (ctg_len_ptr == NULL){
		    fprintf(stderr, "[%s]Warning: fai file does not contain this contig name %s. The contig length will be set to -1.\n", trackReader->ctg);
	    }
	    else{
		    trackReader->ctgLen = *ctg_len_ptr;
	    }
    }
    token = Splitter_getToken(splitter);
    trackReader->s = trackReader->zeroBasedCoors ? atoi(token) : atoi(token) + 1;
    token = Splitter_getToken(splitter);
    trackReader->e = trackReader->zeroBasedCoors ? atoi(token) - 1 : atoi(token);
    while ((token = Splitter_getToken(splitter)) != NULL){
	    trackReader->attrbsLen += 1;
	    if(trackReader->attrbsLen == 1){
		    trackReader->attrbs = malloc(1 * sizeof(char*));
	    }
	    else{ // increase the size of the attrbs if there are more attrbs
		    trackReader->attrbs = realloc(trackReader->attrbs, trackReader->attrbsLen * sizeof(char*));
	    }
	    // save the currect attrb
	    trackReader->attrbs[trackReader->attrbsLen - 1] = malloc(strlen(token) + 1);
	    strcpy(trackReader->attrbs[trackReader->attrbsLen - 1], token);
    }
    Splitter_destruct(splitter);
    return read;
}

int TrackReader_readNextTrackCov(TrackReader *trackReader){
    char *line = malloc(LINE_MAX_SIZE);
    ssize_t read = TrackReader_readLine(trackReader, &line, LINE_MAX_SIZE);
    // if this is the end of the file
    if (read == -1){
         trackReader->ctg[0] = '\0';
	 trackReader->ctgLen = 0;
         trackReader->s = -1;
         trackReader->e = -1;
	 free(line);
         return read;
    }

    //fprintf(stderr,"%s\n", line);
    char* token;
    // free trackReader attrbs to fill it with new ones
    for(int i = 0; i < trackReader->attrbsLen; i++){
                 free(trackReader->attrbs[i]);
    }
    free(trackReader->attrbs);
    trackReader->attrbs = NULL;
    trackReader->attrbsLen = 0;
    if(read == 0){
            fprintf(stderr, "[Warning] line read by TrackReader was empty. Go to the next line!\n");
            free(line);
            return TrackReader_readNextTrackCov(trackReader);
    }
    if (line[0] == '>') { // contig name and size is after '>'
	Splitter* splitter = Splitter_construct(line, ' ');
	token = Splitter_getToken(splitter);
        strcpy(trackReader->ctg, token + 1); // skip '>' and copy
	token = Splitter_getToken(splitter);
	trackReader->ctgLen = atoi(token);
	Splitter_destruct(splitter);
        TrackReader_readNextTrackCov(trackReader);
    }
    else {
	Splitter* splitter = Splitter_construct(line, '\t');
	token = Splitter_getToken(splitter);
        trackReader->s = trackReader->zeroBasedCoors ? atoi(token) - 1 : atoi(token);
	token = Splitter_getToken(splitter);
	trackReader->e = trackReader->zeroBasedCoors ? atoi(token) - 1 : atoi(token);
	while ((token = Splitter_getToken(splitter)) != NULL){
	    //fprintf(stderr, "%s\n",token);
            trackReader->attrbsLen += 1;
	    if(trackReader->attrbsLen == 1){
                    trackReader->attrbs = malloc(1 * sizeof(char*));
            }
            else{// increase the size of the attrbs if there is more attrbs
                    trackReader->attrbs = realloc(trackReader->attrbs, trackReader->attrbsLen * sizeof(char*));
            }
	    // save the current attrb
            trackReader->attrbs[trackReader->attrbsLen - 1] = malloc(strlen(token) + 1);
            strcpy(trackReader->attrbs[trackReader->attrbsLen - 1], token);
	}
	Splitter_destruct(splitter);
    }
    free(line);
    return read;
}