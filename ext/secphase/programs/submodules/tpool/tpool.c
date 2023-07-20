#include "tpool.h"
#include <pthread.h>
#include <string.h>
#include "stdlib.h"

// The tpool code is taken from
// https://nachtimwald.com/2019/04/12/thread-pool-in-c/

work_arg_t *tpool_shallow_copy_work_arg(work_arg_t * arg_src){
    work_arg_t *arg_dest = malloc(sizeof(work_arg_t));
    arg_dest->baq_flag = arg_src->baq_flag;
    arg_dest->consensus = arg_src->consensus;
    arg_dest->indel_threshold = arg_src->indel_threshold;
    arg_dest->min_q = arg_src->min_q;
    arg_dest->min_score = arg_src->min_score;
    arg_dest->prim_margin_score = arg_src->prim_margin_score;
    arg_dest->prim_margin_random = arg_src->prim_margin_random;
    arg_dest->set_q = arg_src->set_q;
    arg_dest->conf_d = arg_src->conf_d;
    arg_dest->conf_e = arg_src->conf_e;
    arg_dest->conf_b = arg_src->conf_b;
    strcpy(arg_dest->inputPath, arg_src->inputPath);
    strcpy(arg_dest->fastaPath, arg_src->fastaPath);
    arg_dest->marker_mode = arg_src->marker_mode;
    arg_dest->variant_ref_blocks_per_contig = arg_src->variant_ref_blocks_per_contig;
    arg_dest->modified_blocks_by_vars_per_contig = arg_src->modified_blocks_by_vars_per_contig;
    arg_dest->modified_blocks_by_marker_per_contig = arg_src->modified_blocks_by_marker_per_contig;
    arg_dest->variant_blocks_all_haps_per_contig = arg_src->variant_blocks_all_haps_per_contig;
    arg_dest->marker_blocks_all_haps_per_contig = arg_src->marker_blocks_all_haps_per_contig;
    arg_dest->reads_modified_by_vars = arg_src->reads_modified_by_vars;
    arg_dest->reads_modified_by_marker = arg_src->reads_modified_by_marker;
    arg_dest->output_log_file = arg_src->output_log_file;
    arg_dest->bam_fo = arg_src->bam_fo;
    arg_dest->mutexPtr = arg_src->mutexPtr;
    arg_dest->alignments = arg_src->alignments;
    arg_dest->alignments_len = arg_src->alignments_len;
    arg_dest->flank_margin = arg_src->flank_margin;
    arg_dest->sam_hdr = arg_src->sam_hdr;
    return arg_dest;
}

work_arg_t *tpool_create_work_arg(bool baq_flag, bool consensus, int indel_threshold, int min_q,
                                  int min_score, double prim_margin_score, double prim_margin_random, int set_q,
                                  double conf_d, double conf_e, double conf_b, char *inputPath,
                                  char *fastaPath, bool marker_mode,
                                  stHash *variant_ref_blocks_per_contig,
                                  stHash *modified_blocks_by_vars_per_contig,
                                  stHash *modified_blocks_by_marker_per_contig,
                                  stHash *variant_blocks_all_haps_per_contig,
                                  stHash *marker_blocks_all_haps_per_contig,
                                  int *reads_modified_by_vars,
                                  int *reads_modified_by_marker,
                                  FILE *output_log_file,
                                  samFile* bam_fo,
                                  pthread_mutex_t *mutexPtr, ptAlignment** alignments, int alignments_len,
                                  int flank_margin, sam_hdr_t *sam_hdr) {
    work_arg_t *arg = malloc(sizeof(work_arg_t));
    arg->baq_flag = baq_flag;
    arg->consensus = consensus;
    arg->indel_threshold = indel_threshold;
    arg->min_q = min_q;
    arg->min_score = min_score;
    arg->prim_margin_score = prim_margin_score;
    arg->prim_margin_random = prim_margin_random;
    arg->set_q = set_q;
    arg->conf_d = conf_d;
    arg->conf_e = conf_e;
    arg->conf_b = conf_b;
    strcpy(arg->inputPath, inputPath);
    strcpy(arg->fastaPath, fastaPath);
    arg->marker_mode = marker_mode;
    arg->variant_ref_blocks_per_contig =variant_ref_blocks_per_contig;
    arg->modified_blocks_by_vars_per_contig = modified_blocks_by_vars_per_contig;
    arg->modified_blocks_by_marker_per_contig = modified_blocks_by_marker_per_contig;
    arg->variant_blocks_all_haps_per_contig = variant_blocks_all_haps_per_contig;
    arg->marker_blocks_all_haps_per_contig = marker_blocks_all_haps_per_contig;
    arg->reads_modified_by_vars = reads_modified_by_vars;
    arg->reads_modified_by_marker = reads_modified_by_marker;
    arg->output_log_file = output_log_file;
    arg->bam_fo = bam_fo;
    arg->mutexPtr = mutexPtr;
    arg->alignments = alignments;
    arg->alignments_len = alignments_len;
    arg->flank_margin = flank_margin;
    arg->sam_hdr = sam_hdr;
    return arg;
}


static tpool_work_t *tpool_work_create(thread_func_t func, void *arg) {
    tpool_work_t *work;

    if (func == NULL)
        return NULL;

    work = malloc(sizeof(*work));
    work->func = func;
    work->arg = arg;
    work->next = NULL;
    return work;
}

static void tpool_work_destroy(tpool_work_t *work) {
    if (work == NULL)
        return;
    free(work);
}

static tpool_work_t *tpool_work_get(tpool_t *tm) {
    tpool_work_t *work;

    if (tm == NULL)
        return NULL;

    work = tm->work_first;
    if (work == NULL)
        return NULL;

    if (work->next == NULL) {
        tm->work_first = NULL;
        tm->work_last = NULL;
    } else {
        tm->work_first = work->next;
    }

    return work;
}

static void *tpool_worker(void *arg) {
    tpool_t *tm = arg;
    tpool_work_t *work;
    while (1) {
        pthread_mutex_lock(&(tm->work_mutex));

        while (tm->work_first == NULL && !tm->stop)
            pthread_cond_wait(&(tm->work_cond), &(tm->work_mutex));

        if (tm->stop)
            break;

        work = tpool_work_get(tm);
        tm->working_cnt++;
        tm->remaining_cnt--;
        pthread_mutex_unlock(&(tm->work_mutex));

        if (work != NULL) {
            work_arg_t *work_arg = work->arg;
            work->func(work_arg);
            tpool_work_destroy(work);
        }

        pthread_mutex_lock(&(tm->work_mutex));
        tm->working_cnt--;
        if (!tm->stop && tm->working_cnt == 0 && tm->work_first == NULL)
            pthread_cond_signal(&(tm->working_cond));
        pthread_mutex_unlock(&(tm->work_mutex));
    }

    tm->thread_cnt--;
    pthread_cond_signal(&(tm->working_cond));
    pthread_mutex_unlock(&(tm->work_mutex));
    return NULL;
}


tpool_t *tpool_create(size_t num) {
    tpool_t *tm;
    pthread_t thread;
    size_t i;

    if (num == 0)
        num = 2;

    tm = calloc(1, sizeof(*tm));
    tm->thread_cnt = num;

    tm->working_cnt = 0;
    tm->remaining_cnt = 0;

    pthread_mutex_init(&(tm->work_mutex), NULL);
    pthread_cond_init(&(tm->work_cond), NULL);
    pthread_cond_init(&(tm->working_cond), NULL);

    tm->work_first = NULL;
    tm->work_last = NULL;

    for (i = 0; i < num; i++) {
        pthread_create(&thread, NULL, tpool_worker, tm);
        pthread_detach(thread);
    }

    return tm;
}

void tpool_destroy(tpool_t *tm) {
    tpool_work_t *work;
    tpool_work_t *work2;

    if (tm == NULL)
        return;

    pthread_mutex_lock(&(tm->work_mutex));
    work = tm->work_first;
    while (work != NULL) {
        work2 = work->next;
        tpool_work_destroy(work);
        work = work2;
    }
    tm->stop = true;
    pthread_cond_broadcast(&(tm->work_cond));
    pthread_mutex_unlock(&(tm->work_mutex));

    tpool_wait(tm);

    pthread_mutex_destroy(&(tm->work_mutex));
    pthread_cond_destroy(&(tm->work_cond));
    pthread_cond_destroy(&(tm->working_cond));

    free(tm);
}

bool tpool_add_work(tpool_t *tm, thread_func_t func, void *arg) {
    tpool_work_t *work;

    if (tm == NULL)
        return false;

    work = tpool_work_create(func, arg);
    if (work == NULL)
        return false;

    pthread_mutex_lock(&(tm->work_mutex));
    if (tm->work_first == NULL) {
        tm->work_first = work;
        tm->work_last = tm->work_first;
    } else {
        tm->work_last->next = work;
        tm->work_last = work;
    }

    tm->remaining_cnt++;
    pthread_cond_broadcast(&(tm->work_cond));
    pthread_mutex_unlock(&(tm->work_mutex));

    return true;
}


void tpool_wait(tpool_t *tm) {
    if (tm == NULL)
        return;

    pthread_mutex_lock(&(tm->work_mutex));
    while (1) {
        if ((!tm->stop && (tm->working_cnt != 0 || tm->work_first != NULL)) || (tm->stop && tm->thread_cnt != 0)) {
            pthread_cond_wait(&(tm->working_cond), &(tm->work_mutex));
        } else {
            break;
        }
    }
    pthread_mutex_unlock(&(tm->work_mutex));
}

