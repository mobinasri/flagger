#include <stdbool.h>
#include <stddef.h>
#include <pthread.h>
#include "sonLib.h"
#include <sys/types.h>
#include "ptVariant.h"
#include "ptBlock.h"
#include "ptAlignment.h"
#include "sam.h"


#ifndef __THREAD_POOL_H__
#define __THREAD_POOL_H__

typedef struct tpool tpool_t;

typedef void (*thread_func_t)(void *arg);

struct tpool_work {
    thread_func_t func;
    void *arg;
    struct tpool_work *next;
};
typedef struct tpool_work tpool_work_t;

typedef struct work_arg_t {
    bool baq_flag;
    bool consensus;
    int indel_threshold;
    int min_q;
    int min_score;
    double prim_margin_score;
    double prim_margin_random;
    int set_q;
    double conf_d;
    double conf_e;
    double conf_b;
    char inputPath[1000];
    char fastaPath[1000];
    bool marker_mode;
    stHash *variant_ref_blocks_per_contig;
    stHash *modified_blocks_by_vars_per_contig;
    stHash *modified_blocks_by_marker_per_contig;
    stHash *variant_blocks_all_haps_per_contig;
    stHash *marker_blocks_all_haps_per_contig;
    int *reads_modified_by_vars;
    int *reads_modified_by_marker;
    FILE *output_log_file;
    samFile* bam_fo;
    pthread_mutex_t *mutexPtr;
    ptAlignment** alignments;
    int alignments_len;
    int flank_margin;
    sam_hdr_t *sam_hdr;
} work_arg_t;

struct tpool {
    tpool_work_t *work_first;
    tpool_work_t *work_last;
    pthread_mutex_t work_mutex;
    pthread_cond_t work_cond;
    pthread_cond_t working_cond;
    size_t working_cnt;
    size_t remaining_cnt;
    size_t thread_cnt;
    bool stop;
    int idx;
};


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
                                  pthread_mutex_t *mutexPtr,
                                  ptAlignment** alignments, int alignments_len,
                                  int flank_margin, sam_hdr_t *sam_hdr);

work_arg_t *tpool_shallow_copy_work_arg(work_arg_t * arg_src);

tpool_t *tpool_create(size_t num);

void tpool_destroy(tpool_t *tm);

bool tpool_add_work(tpool_t *tm, thread_func_t func, void *arg);

void tpool_wait(tpool_t *tm);

#endif /* __THREAD_POOL_H__ */
