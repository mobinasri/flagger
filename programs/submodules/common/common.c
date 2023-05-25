#include "common.h"


bool file_exists (char *filename) {
    struct stat   buffer;
    return (stat (filename, &buffer) == 0);
}

int min(int a, int b){
	return a < b ? a : b;
}

int max(int a, int b){
        return b < a ? a : b;
}

uint8_t maxCharArray(uint8_t* a, int len){
	uint8_t m = 0;
	for(int i=0; i<len; i++){
		m = a[i] < m ? m : a[i];
	}
	return m;
}

uint8_t minCharArray(uint8_t* a, int len){
        uint8_t m = 255;
        for(int i=0; i<len; i++){
                m = m < a[i] ? m : a[i];
        }
        return m;
}

Splitter* Splitter_construct(char* str, char delimiter){
        Splitter* splitter = malloc(sizeof(Splitter));
        splitter->str = malloc((strlen(str) + 1) * sizeof(char));
        strcpy(splitter->str, str);
        splitter->token = malloc((strlen(str) + 1) * sizeof(char));
	splitter->delimiter = delimiter;
        splitter->offset = 0;
}

void Splitter_destruct(Splitter* splitter){
        free(splitter->str);
        free(splitter->token);
        free(splitter);
}

char* Splitter_getToken(Splitter* splitter){
	int i = splitter->offset;
	int j = 0;
	while(splitter->str[i] != '\0' && splitter->str[i] != splitter->delimiter){
		splitter->token[j] = splitter->str[i];
		j++; i++;
	}
	splitter->offset = splitter->str[i] == splitter->delimiter ? i + 1 : i;
	splitter->token[j] = '\0';
	if (j == 0){ // end of the string
		free(splitter->token);
		splitter->token = NULL;
	}
	return splitter->token;
}
