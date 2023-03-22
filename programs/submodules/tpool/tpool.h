#include <stdbool.h>
#include <stddef.h>
#include <pthread.h>
#include "hmm.h"

#ifndef __THREAD_POOL_H__
#define __THREAD_POOL_H__

typedef struct tpool tpool_t;
typedef struct tpool_batch tpool_batch_t;

typedef void (*thread_func_t)(void *arg);

struct tpool_work {
    thread_func_t      func;
    void              *arg;
    struct tpool_work *next;
};
typedef struct tpool_work tpool_work_t;

typedef struct work_arg_t {
	Chunk* templateChunkIdx;
	HMM* model;
	Batch* batch;
	char dir[200];
	char name[200];
}work_arg_t;

struct tpool {
    tpool_work_t    *work_first;
    tpool_work_t    *work_last;
    pthread_mutex_t  work_mutex;
    pthread_cond_t   work_cond;
    pthread_cond_t   working_cond;
    size_t           working_cnt;
    size_t           thread_cnt;
    bool             stop;
    int              idx;
};


typedef struct tpool_batch{
	Batch* batch;
	tpool_t* tm;
} tpool_batch;


tpool_t *tpool_create(size_t num, Batch** batches);
void tpool_destroy(tpool_t *tm);

bool tpool_add_work(tpool_t *tm, thread_func_t func, void *arg);
void tpool_wait(tpool_t *tm);

#endif /* __THREAD_POOL_H__ */