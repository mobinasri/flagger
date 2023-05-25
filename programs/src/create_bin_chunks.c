#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include "sonLib.h"
#include "common.h"
#include "hmm.h"
#include "tpool.h"
#include <pthread.h>
#include <time.h>

static struct option long_options[] =
        {
                {"inputCov",  required_argument, NULL, 'i'},
                {"threads",   required_argument, NULL, 't'},
                {"chunkLen",  required_argument, NULL, 'l'},
                {"windowLen", required_argument, NULL, 'w'},
                {NULL,        0,                 NULL, 0}
        };


int main(int argc, char *argv[]) {
    int c;
    char covPath[200];
    int nThreads;
    int chunkLen;
    int windowLen = 1000;
    char *program;
    (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
    while (~(c = getopt_long(argc, argv, "i:t:l:w:h", long_options, NULL))) {
        switch (c) {
            case 'i':
                strcpy(covPath, optarg);
                break;
            case 't':
                nThreads = atoi(optarg);
                break;
            case 'l':
                chunkLen = atoi(optarg);
                break;
            case 'w':
                windowLen = atoi(optarg);
                break;
            default:
                if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
            help:
                fprintf(stderr, "\nUsage: %s\n", program);
                fprintf(stderr, "Options:\n");
                fprintf(stderr, "         --inputCov, -i         path to .cov file\n");
                fprintf(stderr, "         --threads, -t         number of threads\n");
                fprintf(stderr, "         --chunkLen, -l         chunk length\n");
                fprintf(stderr, "         --windowLen, -w         window length (default = 100)\n");
                return 1;
        }
    }
    sprintf(outputPath, "%s.chunks.l_%d.w_%d.bin", covPath, chunkLen, windowLen);
    FILE *fp = fopen(outputPath, "wb+");
    Batch *batch = Batch_construct(covPath, chunkLen, nThreads, 0, windowLen);
    batch->templateChunkIdx = 0;
    batch->nThreadChunks = 0;
    fwrite(&(batch->chunkLen), sizeof(int32_t), 1, fp); // first 4 bytes
    fwrite(&(batch->windowLen), sizeof(int32_t), 1, fp); // second 4 bytes
    while (Batch_readThreadChunks(batch)) {
        for (int t = 0; t < batch->nThreadChunks; t++) {
            Chunk *chunk = batch->threadChunks[t];
            fprintf(stderr, "Writing the chunk %s:%d-%d\n", chunk->ctg, chunk->s, chunk->e);
            // write the length of contig name + null character
            int32_t ctgNameLen = strlen(chunk->ctg) + 1;
            fwrite(&ctgNameLen, sizeof(int32_t), 1, fp);
            // write the contig name
            fwrite(chunk->ctg, sizeof(char), ctgNameLen, fp);
            // write start and e
            fwrite(&(chunk->s), sizeof(int32_t), 1, fp);
            fwrite(&(chunk->e), sizeof(int32_t), 1, fp);
            //write actual size of chunk
            fwrite(&(chunk->seqLen), sizeof(int32_t), 1, fp);
            int32_t *cov1Array = malloc(chunk->seqLen * sizeof(int32_t));
            int32_t *cov2Array = malloc(chunk->seqLen * sizeof(int32_t));
            for (int i = 0; i < chunk->seqLen; i++) {
                cov1Array[i] = (int) chunk->seqEmit[i]->data[0];
                cov2Array[i] = (int) chunk->seqEmit[i]->data[1];
            }
            // write cov1 and then cov2 array
            fwrite(cov1Array, sizeof(int32_t), chunk->seqLen, fp);
            fwrite(cov2Array, sizeof(int32_t), chunk->seqLen, fp);
        }
    }
    fflush(fp);
    fclose(fp);
    fprintf(stderr, "Done writing %s!\n", outputPath);
    Batch_destruct(batch);
}
