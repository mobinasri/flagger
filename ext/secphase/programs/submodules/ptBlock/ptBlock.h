#ifndef PT_BLOCK_H
#define PT_BLOCK_H

#include <assert.h>
#include <math.h>
#include <float.h>
#include "common.h"
#include "sonLib.h"
#include "stdlib.h"
#include "ptAlignment.h"

/*
 * @abstract Structure for saving a block
 * @field rfs           Start coordinate on reference (0-based inclusive)
 * @field rfe           End coordinate on reference (0-based inclusive)
 * @field sqs           Start coordinate on sequence (0-based inclusive)
 * @field sqe           End coordinate on sequence (0-based inclusive)
 * @field rds_f         Start coordinate on read's forward strand (0-based inclusive)
 * @field rde_f         End coordinate on read's forward strand (0-based inclusive)
 * @field data		Pointer to some data attached to the block
 * @field free_data	Function for freeing the memory allocated to data
 */
typedef struct {
    int rfs; // ref start
    int rfe; // ref end
    int sqs; // seq start
    int sqe; // seq end
    int rds_f;
    int rde_f;
    void *data;

    void (*destruct_data)(void *);

    void *(*copy_data)(void *);

    void (*extend_data)(void *, void *);
} ptBlock;


/* Construct a ptBlock structure
 * Note that this constructor only receives the coordinates
 * and related data should be added by another function; ptBlock_set_data()
 */
ptBlock *ptBlock_construct(int rfs, int rfe, int sqs, int sqe, int rds_f, int rde_f);


/* Set data and the related functions in a ptBlock structure.
 * It also receives three functions
 *  1- free_data for freeing the memory allocated to the data.
 *  2- copy_data for copying data
 *  3- append_data for appending new data to what is already saved in the block
 */
void ptBlock_set_data(ptBlock *block, void *data, void (*destruct_data)(void *), void *(*copy_data)(void *),
                      void (*extend_data)(void *, void *));


void ptBlock_extend_data(ptBlock *block, void *data);


void ptBlock_destruct_data(ptBlock *block);


void *ptBlock_copy_data(ptBlock *block);


/* Make a copy of a ptBlock structure
 */
ptBlock *ptBlock_copy(ptBlock *block);


/* Destruct a ptBlock structure
 * It first frees the data augmented to the block
 * then frees the ptBlock struct
*/
void ptBlock_destruct(ptBlock *block);


/// Get Functions ////

int ptBlock_get_rfs(ptBlock *block);

int ptBlock_get_rfe(ptBlock *block);

int ptBlock_get_sqs(ptBlock *block);

int ptBlock_get_sqe(ptBlock *block);

int ptBlock_get_rds_f(ptBlock *block);

int ptBlock_get_rde_f(ptBlock *block);


//// Set Functions ////

void ptBlock_set_rfs(ptBlock *block, int rfs);

void ptBlock_set_rfe(ptBlock *block, int rfe);

void ptBlock_set_sqs(ptBlock *block, int sqs);

void ptBlock_set_sqe(ptBlock *block, int sqe);

void ptBlock_set_rds_f(ptBlock *block, int rds_f);

void ptBlock_set_rde_f(ptBlock *block, int rde_f);


// Compare two blocks by rfs
int ptBlock_cmp_rfs(const void *a, const void *b);

// Compare two blocks by rds_f
int ptBlock_cmp_rds_f(const void *a, const void *b);

// Compare two blocks by sqs
int ptBlock_cmp_sqs(const void *a, const void *b);


void ptBlock_sort_stHash_by_rfs(stHash *blocks_per_contig);

/**
 * Parse a bed file and save the tracks in a hash table
 * Each contig/chromosome is saved as a key
 * All blocks in each contig/chromosome are save in a single list as the corresponding value
 *
 * @param bed_path  Path to a BED file
 * @return  stHash table with contigs as keys and tracks as values
 *
 */
stHash *ptBlock_parse_bed(char *bed_path);


/**
 * Merge blocks (should be sorted by stList_sort)
 * Since merging can be performed by either rfs/rfe or rds_f/rde_f this function
 * takes two functions as inputs for getting start and end of each block
 * the functions that can be used for this purpose are defined above (starting with ptBlock_get_)
 * This enables merging block by different criteria
 * it also takes set_end function for updating the end location of merged block if it had
 * overlap with any new block during the merging process
 *
 *
 * @param blocks  stList of blocks each of which saved as a ptBlock
 * @param  get_start   Function for getting the start coordinate of each block
 * @param  get_end   Function for getting the end coordinate of each block
 * @param set_end  Function for updating the end coordinate of each block
 *
 * @return  merged blocks
 *
 * @Note    input blocks should be sorted based on the start coordinates so the sorting criteria
 *          should match the merging criteria
 *          For example if blocks are sorted by ptBlock_cmp_rfs you can merge by ptBlock_get_rfs/rfe
 *          and if they are sorted by ptBlock_cmp_rds_f you can merge by ptBlock_get_rds_f/rde_f
 *
 */
stList *ptBlock_merge_blocks(stList *blocks,
                             int (*get_start)(ptBlock *), int (*get_end)(ptBlock *),
                             void (*set_end)(ptBlock *, int));


/**
 * Merge blocks for each contig (each list should be sorted by stList_sort. ptBlock_parse_bed
 * sorts blocks by default)
 * Since merging can be performed by either rfs/rfe or rds_f/rde_f this function
 * takes two functions as inputs for getting start and end of each block
 * the functions that can be used for this purpose are defined above (starting with ptBlock_get_)
 * This enables merging block by different criteria
 * it also takes set_end function for updating the end location of merged block if it had
 * overlap with any new block uring the merging process
 *
 *
 * @param blocks    stHash table of blocks (can be output of ptBlock_parse_bed).
 *                  Keys are contigs/chromosomes and values are stLists
 *                  of blocks each of which saved as a ptBlock
 * @param  get_start   Function for getting the start coordinate of each block
 * @param  get_end   Function for getting the end coordinate of each block
 * @param set_end  Function for updating the end coordinate of each block
 *
 * @return  merged blocks saved in a stHash table
 *
 * @Note    input blocks should be sorted based on the start coordinates so the sorting criteria
 *          should match the merging criteria
 *          For example if blocks are sorted by ptBlock_cmp_rfs you can merge by ptBlock_get_rfs/rfe
 *          and if they are sorted by ptBlock_cmp_rds_f you can merge by ptBlock_get_rds_f/rde_f
 *
 */
stHash *ptBlock_merge_blocks_per_contig(stHash *blocks_per_contig,
                                        int (*get_start)(ptBlock *), int (*get_end)(ptBlock *),
                                        void (*set_end)(ptBlock *, int));

stHash *ptBlock_merge_blocks_per_contig_by_rf(stHash *blocks_per_contig);

stHash *ptBlock_merge_blocks_per_contig_by_rd_f(stHash *blocks_per_contig);

stHash *ptBlock_merge_blocks_per_contig_by_sq(stHash *blocks_per_contig);

void ptBlock_add_alignment(stHash *blocks_per_contig, ptAlignment *alignment);

void ptBlock_save_in_bed(stHash *blocks_per_contig, char* bed_path);

int ptBlock_get_total_number(stHash *blocks_per_contig);

void ptBlock_add_blocks_by_contig(stHash *blocks_per_contig, char* contig, stList* blocks_to_add);


#endif /* PT_BLOCK_H */

