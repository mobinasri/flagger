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
#include <pthread.h>

typedef enum StatType {FORWARD, BACKWARD, POSTERIOR, EMISSION, TRANSITION} StatType;
const char* const COMP_COLORS[] = {"162,0,37", "250,104,0", "0,138,0", "0,138,0", "0,138,0", "170,0,255", "128,128,128"};
const char* const COMP_NAMES[] = {"Err", "Dup", "Hap1", "Hap2", "Hap3", "Col", "Unk"};

struct stat st = {0};

pthread_mutex_t mutex;

typedef struct Arguments{
        HMM* model;
        EM* em;
	char name[200];
        char path[200];
        pthread_mutex_t* mutexPtr;
} Arguments;

Arguments* Arguments_construct(HMM* model, EM* em, char* name, char* path, pthread_mutex_t* mutexPtr){
        Arguments* args = malloc(sizeof(Arguments));
        args->model = model;
        args->em = em;
	args->name[0] = '\0';
	args->path[0] = '\0';
	if (name != NULL)
        	strcpy(args->name, name);
	if (path != NULL)
        	strcpy(args->path, path);
        args->mutexPtr = mutexPtr;
}



void saveHMMStats(StatType statType, char* path, HMM* model){
        FILE * fp;

        fp = fopen(path, "w+");
        if (fp == NULL){
                printf("Couldn't open %s\n", path);
                exit(EXIT_FAILURE);
        }

        // Print transition or emission probs
        for(int r=0; r < model->nClasses; r++){
		fprintf(fp, "#r=%d\n------------\n",r);
		if (statType == TRANSITION){
			char* matrixStr = MatrixDouble_toString(model->trans[r]);
			fprintf(fp, "%s\n\n", matrixStr);
		}
		else if (statType == EMISSION){
			for(int c=0; c < model->nComps; c++){
				fprintf(fp, "##c=%d\n",c);
				Gaussian* gaussian = model->emit[r][c];
				for(int m = 0; m < gaussian->n; m++){
					fprintf(fp, "##m=%d\n",m);
					char* vecStr = VectorDouble_toString(gaussian->mu[m]);
					fprintf(fp, "mu = \n%s\n\n", vecStr);
					char* matrixStr = MatrixDouble_toString(gaussian->cov[m]);
                                	fprintf(fp, "cov = \n%s\n\n", matrixStr);
				}
			}
		}
	}
	fclose(fp);
}


void saveEMStats(StatType statType, char* path, EM* em){
	FILE * fp;

        fp = fopen(path, "w+");
        if (fp == NULL){
                printf("Couldn't open %s\n", path);
                exit(EXIT_FAILURE);
        }

	// Print header
        fprintf(fp, "#i\t");
	for(int c =0; c < em->nComps; c++){
		fprintf(fp, "comp_%d\t", c);
	}
	fprintf(fp, "Scale\n");

	double* p;
	// Print numbers one position per each line
	for(int i=0; i < em->seqLength; i++){
		switch(statType){
			case FORWARD:
                		p = getForward(em, i);
				break;
			case BACKWARD:
				p = getBackward(em, i);
				break;
			case POSTERIOR:
				p = getPosterior(em, i);
				break;
		}
                fprintf(fp, "%d\t", i);
                for(int c =0; c < em->nComps; c++){
                        fprintf(fp, "%.3E\t", p[c]);
                }
		fprintf(fp, "%.3E\t\n", em->scales[i]);
		free(p);
        }
	fclose(fp);
}

void initMuSixComps(VectorDouble*** mu, int coverage, int* nMixtures){
	double errCov = coverage > 20 ? (double) coverage / 20 : 1 ;
	double dupCov = (double) coverage / 2.0;
	double hapCov = (double) coverage;
	double colCov = (double) coverage * 2.0;

	// Erroneous
	mu[0][0]->data[0] = 1;
	mu[0][0]->data[1] = 1;
	// Duplicated
	mu[1][0]->data[0] = dupCov;
	mu[1][0]->data[1] = 0;
	// Haploid Hom
	mu[2][0]->data[0] = hapCov;
        mu[2][0]->data[1] = 0;
	// Haploid Het
	mu[3][0]->data[0] = hapCov;
	mu[3][0]->data[1] = hapCov;
	// Collapsed
	mu[4][0]->data[0] = colCov;
	mu[4][0]->data[1] = colCov;
}

HMM* makeAndInitModel(int* coverages, int nClasses, int nComps, int nEmit, int* nMixtures){

	VectorDouble**** mu = malloc(nClasses * sizeof(VectorDouble***));
	for(int r = 0; r < nClasses; r++){
		mu[r] = malloc(nComps * sizeof(VectorDouble**));
		for(int c = 0; c < nComps; c++){
			mu[r][c] = VectorDouble_constructArray1D(nMixtures[c], nEmit);
		}
	}
	for (int r = 0; r < nClasses; r++){
		initMuSixComps(mu[r], coverages[r], nMixtures);
	}
   	HMM* model = HMM_construct(nClasses, nComps, nEmit, nMixtures, mu);
	//model->emit[1][1]->cov->data[1][1] = model->emit[1][2]->mu->data[0] / 8.0;
	//model->emit[1][2]->cov->data[1][1] = 4.0 * model->emit[1][2]->mu->data[0];
	//model->emit[1][3]->cov->data[1][1] = 8.0 * model->emit[1][2]->mu->data[0];
	for(int r = 0; r < nClasses; r++){
                for(int c = 0; c < nComps; c++){
                        VectorDouble_destructArray1D(mu[r][c], nMixtures[c]);
                }
		free(mu[r]);
        }
	free(mu);
	return model;
}

EM** Batch_buildEmArray(Batch* batch, HMM* model){
	EM** emArray = (EM**) malloc(batch->nThreadChunks * sizeof(EM*));
	Chunk* chunk;
	for(int t=0; t < batch->nThreadChunks; t++){
		chunk = batch->threadChunks[t];
		emArray[t] = EM_construct(chunk->seqEmit, chunk->seqClass, chunk->seqLen, model);
	}
	return emArray;
}

void* trainModelSaveStats(void* args_){
	Arguments* args = (Arguments*) args_;
	HMM* model = args->model;
	EM* em = args->em;
	char* name = args->name;
	char* dir = args->path;
	pthread_mutex_t* mutexPtr = args->mutexPtr;
        char path[200];
        //for(int itr=1; itr <= nItr; itr++){
		//fprintf(stderr,"\t\tIteration %d\n", itr);
                //fprintf(stderr,"\t\t\tRun forward\n");
		runForward(model, em);
		//fprintf(stderr,"\t\t\tRun borward\n");
                runBackward(model, em);

                //sprintf(path, "%s/forward_%s.txt", dir, name);
                //saveEMStats(FORWARD, path, em);
                //sprintf(path, "%s/backward_%s.txt", dir, name);
                //saveEMStats(BACKWARD, path, em);
                //sprintf(path, "%s/posterior_%d.txt", dir, itr);
                //saveEMStats(POSTERIOR, path, em);

		//fprintf(stderr, "\t\t\tUpdate sufficient Stats ...\n");
		pthread_mutex_lock(mutexPtr);
                updateSufficientStats(model, em);
		pthread_mutex_unlock(mutexPtr);

                //estimateParameters(model);
                //resetSufficientStats(model);

                //sprintf(path, "%s/transition_%d.txt", dir, itr);
                //saveHMMStats(TRANSITION, path, model);
                //sprintf(path, "%s/emission_%d.txt", dir, itr);
                //saveHMMStats(EMISSION, path, model);
        //}
}

void* infer(void* args_){
        Arguments* args = (Arguments*) args_;
        HMM* model = args->model;
        EM* em = args->em;
        char path[200];
        runForward(model, em);
        runBackward(model, em);

                //sprintf(path, "%s/forward_%d.txt", dir, itr);
                //saveEMStats(FORWARD, path, em);
                //sprintf(path, "%s/backward_%d.txt", dir, itr);
                //saveEMStats(BACKWARD, path, em);
                //sprintf(path, "%s/posterior_%d.txt", dir, itr);
                //saveEMStats(POSTERIOR, path, em);
        //}
}


void Batch_trainModelSaveStats(Batch* batch, int batchIdx, HMM* model, char* dir){
	EM** emArray =  Batch_buildEmArray(batch, model);
	char dir1[100];
	char dir2[100];
	char name[200];
	EM* em;
	pthread_t* tids = malloc(batch->nThreadChunks * sizeof(pthread_t));
	Arguments** args = malloc(batch->nThreadChunks * sizeof(Arguments*));
	for(int t=0; t < batch->nThreadChunks; t++){
		em = emArray[t];
		/*sprintf(dir1, "%s/batch_%d", dir, batchIdx);
		sprintf(dir2, "%s/batch_%d/chunk_%d", dir, batchIdx, t);
		if (stat(dir1, &st) == -1) {
			mkdir(dir1, 0700);
		}
		if (stat(dir2, &st) == -1) {
                        mkdir(dir2, 0700);
                }*/
		Chunk* chunk = batch->threadChunks[t];
		fprintf(stderr, "Batch %d, Chunk %d, %s: [%d-%d] \n", batchIdx, t, chunk->ctg, chunk->s, chunk->e);
		sprintf(name, "ba%d_ch%d", batchIdx, t);
		args[t] = Arguments_construct(model, em, name, dir, &mutex);
		fprintf(stderr, "Args is created\n");
		pthread_create(&tids[t], NULL, trainModelSaveStats, (void*) args[t]);
		fprintf(stderr, "Thread %d is running\n", t);
	}
	for(int t=0; t < batch->nThreadChunks; t++){
		assert(pthread_join(tids[t], NULL) == 0 );
		fprintf(stderr, "Thread %d is finished\n", t);
		free(args[t]);
		EM_destruct(emArray[t]);
	}
	free(args);
	free(emArray);
	free(tids);
}

uint8_t getCompIdx(EM* em, int loc){
	double* p = getPosterior(em, loc);
	uint8_t idx = 0;
	for(int c =1; c < em->nComps; c++){
		idx = p[idx] < p[c] ? c : idx;
        }
        free(p);
	return idx;
}

void Batch_inferSaveOutput(Batch* batch, int batchIdx, HMM* model, FILE* outputFile){
	EM** emArray =  Batch_buildEmArray(batch, model);
        char dir1[100];
        char dir2[100];
        EM* em;
        pthread_t* tids = malloc(batch->nThreadChunks * sizeof(pthread_t));
        for(int t=0; t < batch->nThreadChunks; t++){
                em = emArray[t];
                Chunk* chunk = batch->threadChunks[t];
                fprintf(stderr, "Batch %d, Chunk %d, %s: [%d-%d] \n", batchIdx, t, chunk->ctg, chunk->s, chunk->e);
                Arguments* args = Arguments_construct(model, em, "infer", NULL, NULL);
                pthread_create(&tids[t], NULL, infer, (void*) args);
                fprintf(stderr, "Thread %d is running\n", t);
        }
	int s = -1;
	int e = -1;
	int compIdx = -1;
	int preCompIdx = -1;
	int windowLen = batch->windowLen;
        for(int t=0; t < batch->nThreadChunks; t++){
                assert(pthread_join(tids[t], NULL) == 0 );
                fprintf(stderr, "Thread %d is finished\n", t);
		fprintf(stderr, "Writing the output...\n");
		Chunk* chunk = batch->threadChunks[t];
		s = chunk->s;
		e = chunk->s;
		preCompIdx = -1;
		for (int i=0; i < chunk->seqLen; i++){
			compIdx = getCompIdx(emArray[t], i);
			// if component changed write the block
			if (preCompIdx != -1 && preCompIdx != compIdx){
				fprintf(outputFile, "%s\t%d\t%d\t%s\t0\t+\t%d\t%d\t%s\n", 
						    chunk->ctg, s, e, COMP_NAMES[preCompIdx],
						    s, e, COMP_COLORS[preCompIdx]);
				s = e;
			}
			e = i == chunk->seqLen - 1 ? chunk->e + 1: e + windowLen; // maybe the last window is not complete
			preCompIdx = compIdx;
		}
		fprintf(outputFile, "%s\t%d\t%d\t%s\t0\t+\t%d\t%d\t%s\n", 
                                    chunk->ctg, s, e, COMP_NAMES[preCompIdx],
                                    s, e, COMP_COLORS[preCompIdx]); // there will be an unwritten component finally
		fprintf(stderr, "Done writing!\n");
                EM_destruct(emArray[t]);
        }
        free(emArray);
        free(tids);
}

int main(int argc, char *argv[]) {
   int c;
   char covPath[200];
   char covIndexPath[200];
   char outputDir[200];
   char trackName[200];
   int nThreads;
   int chunkLen;
   int maxNBatch=10;
   int windowLen=100;
   int nIteration=1;
   char *program;
   (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
   while (~(c=getopt(argc, argv, "c:o:t:n:w:b:l:r:h"))) {
		switch (c) {
			case 'c':
                                strcpy(covPath, optarg);
                                break;
			case 't':
                                nThreads = atoi(optarg);
                                break;
			case 'o':
                                strcpy(outputDir, optarg);
                                break;
			case 'l':
				chunkLen = atoi(optarg);
				break;
			case 'w':
				windowLen = atoi(optarg);
				break;
			case 'b':
                                maxNBatch = atoi(optarg);
                                break;
			case 'n':
				nIteration = atoi(optarg);
				break;
			case 'r':
				strcpy(trackName, optarg);
				break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
			help:	
				fprintf(stderr, "\nUsage: %s  -a <COV/BED_FILE> -b <COV/BED_FILE> -o <OUT_COV_FILE> \n", program);
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -c         path to .cov file\n");
				fprintf(stderr, "         -t         number of threads\n");
				fprintf(stderr, "         -l         chunk length\n");
				fprintf(stderr, "         -b         max number of batches\n");
				fprintf(stderr, "         -n         number of iterations\n");
				fprintf(stderr, "         -w         window length (default = 100)\n");
				fprintf(stderr, "         -r         track name\n");
				fprintf(stderr, "         -o         output dir\n");
				return 1;	
		}		
   }
   int nClasses = 3;
   int nComps = 5;
   int nEmit = 2 ;
   int* coverages = malloc(nClasses * sizeof(int));
   coverages[0] = 15;
   coverages[1] = 20;
   coverages[2] = 30;
   int* nMixtures = malloc(nComps * sizeof(int));
   nMixtures[0] = 1;
   nMixtures[1] = 1;
   nMixtures[2] = 1;
   nMixtures[3] = 1;
   nMixtures[4] = 1;
   HMM* model = makeAndInitModel(coverages, nClasses, nComps, nEmit, nMixtures);

   for (int r = 0; r < nClasses; r++){
	   for(int c = 0; c < nComps; c++){
		   Gaussian* gaussian = model->emit[r][c];
                   for(int m = 0; m < gaussian->n; m++){
		   	char* muStr = VectorDouble_toString(model->emit[r][c]->mu[m]);
		   	char* covStr = MatrixDouble_toString(model->emit[r][c]->cov[m]);
		   	fprintf(stderr, "r=%d, c=%d, m=%d\n%s\n%s\n\n", r, c, m, muStr, covStr);
		   }
	   }
   }

   char outputPath[200];
   int i = 0;
   Batch* batch = Batch_construct(covPath, chunkLen, nThreads, nEmit, windowLen);
   int batchIdx = 0 ;
   for(int itr=0; itr < nIteration; itr ++){
   	batch->templateChunkIdx = 0;
	batch->nThreadChunks = 0;
   	fprintf(stderr, "Batch is built\n");
   	batchIdx = 0;
	while(Batch_readThreadChunks(batch)){
		fprintf(stderr, "Run EM for batch %d\n", batchIdx);
                Batch_trainModelSaveStats(batch, batchIdx, model, outputDir);
		batchIdx += 1;
   	}

	fprintf(stderr, "Estimating parameters\n");
   	estimateParameters(model);
   	resetSufficientStats(model);

	fprintf(stderr, "Saving HMM stats\n");
	sprintf(outputPath, "%s/transition_%d.txt", outputDir, itr);
	fprintf(stderr, "%s\n", outputPath);
   	saveHMMStats(TRANSITION, outputPath, model);
	sprintf(outputPath, "%s/emission_%d.txt", outputDir,itr);
   	saveHMMStats(EMISSION, outputPath, model);
   }

   //Run inference
   sprintf(outputPath, "%s/%s.flagger.bed", outputDir, trackName);
   FILE* fp = fopen(outputPath, "w+");
   fprintf(fp, "track name=%s visibility=1 itemRgb=\"On\"\n", trackName);
   batchIdx = 0;
   batch->templateChunkIdx = 0;
   batch->nThreadChunks = 0;
   while(Batch_readThreadChunks(batch)){
	   fprintf(stderr, "[Inference] Running EM for batch %d\n", i);
	   Batch_inferSaveOutput(batch, batchIdx, model, fp);
	   batchIdx += 1;
   }
   fflush(fp);
   fclose(fp);
   Batch_destruct(batch);
   HMM_destruct(model);
}
