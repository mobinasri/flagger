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


/*
 * @abstract Structure for iterating over the blocks saved in a stHash table. Table should have
 *           contig names as keys and a stList of blocks as values
 * @field blocks_per_contig     The stHash table of blocks
 * @field ctg_list              A list of contig names created once "ptBlockItrPerContig_construct" is called
 * @field ctg_index             The index of the current contig in ctg_list
 * @field block_index           The index of the current block in the block list related to the current contig
 */
typedef struct {
    void* blocks_per_contig;

    stList* ctg_list;
    int ctg_index;
    int block_index;
} ptBlockItrPerContig;

/* Construct a ptBlock structure
 * Note that this constructor only receives the coordinates
 * and related data should be added by another function; ptBlock_set_data()
 */
ptBlock *ptBlock_construct(int rfs, int rfe, int sqs, int sqe, int rds_f, int rde_f);

/* Construct a ptBlock structure with count
 *
 */
ptBlock *ptBlock_construct_with_count(int rfs, int rfe, int sqs, int sqe, int rds_f, int rde_f, int count);

int ptBlock_get_count(ptBlock* block);

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

// the three functions for keeping count data in ptBlock structs

void extend_count_data(void* dest_, void *src_);


void destruct_count_data(void *src_);


void *copy_count_data(void *src_);

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


/* Compare two ptBlock structures
 *
 * @param block_1       the first block
 * @param block_2       the second block
 * @return              true if the ref/read coordinates of the two blocks were exactly the same
 */
bool ptBlock_is_equal(ptBlock* block_1, ptBlock* block_2);


/* Compare two lists of ptBlock structures
 *
 * @param blocks_1      the first list of blocks
 * @param blocks_2      the second list of blocks
 * @return              true if the ref/read coordinates of the two block lists were exactly the same
 */
bool ptBlock_is_equal_stList(stList* blocks_1, stList* blocks_2);


/* Compare two stHash tables of blocks
 *
 * @param blocks_1      the first table of blocks
 * @param blocks_2      the second table of blocks
 * @return              true if the contig names and ref/read coordinates of the blocks in the two tables were exactly the same
 */
bool ptBlock_is_equal_stHash(stHash* blocks_per_contig_1, stHash* blocks_per_contig_2);


// Functions for block iterator


/* Construct a ptBlockItrPerContig structure
 *
 * This function creates an iterator that can be used
 * for iterating over the blocks saved in a stHash table
 * having contig names as keys. It will keep a reference to the
 * table and also the indices related to the current contig and block
 * to keep track of the previously iterated blocks.
 *
 * @param blocks_per_contig     stHash table of blocks (each value is a stList of blocks)
 * @return block_iter           the block iterator
 */
ptBlockItrPerContig *ptBlockItrPerContig_construct(stHash *blocks_per_contig);



/* Return the next block and update the contig name
 *
 * @param block_iter        the block iterator
 * @param ctg_name_ptr      a reference to the contig name (to be able to update in place)
 * @return block            ptBlock of the next block
 */
ptBlock* ptBlockItrPerContig_next(ptBlockItrPerContig *block_iter, char** ctg_name_ptr);


/* Destruct a ptBlockItrPerContig structure
 *
 * This function destructs the block iterator. It will NOT
 * destruct the table of blocks that was used for creating the iterator
 */
void ptBlockItrPerContig_destruct(ptBlockItrPerContig *blockItr);


/* Add a block to the blocks table in which each key is the corresponding contig name
 * and each value is a stList of blocks
 *
 * @param blocks_per_contig     stHash table of blocks (each value is a stList of blocks)
 * @param block                 the block to add
 * @param ctg_name              the contig name
 */
void ptBlock_add_block_to_stList_table(stHash* blocks_per_contig, ptBlock* block, char* ctg_name);


/* This function recieves a stHash table of blocks and returns a stList of smaller tables. Each 
 * smaller table is batch that contains roughly 1/split_number of the whole blocks (in bases). 
 * It is useful for parsing alignments from an indexed bam file through multi-threading especially
 * when the process is IO bound.
 *
 * @param blocks_per_contig     stHash table of blocks (each value is a stList of blocks)
 * @param split_number		Number of batches to create (length of the output list)
 * @return a stList of stHash tables. Each table is a batch of regions including about 1/split_number
 *         of the given blocks
 */
stList* ptBlock_split_into_batches(stHash *blocks_per_contig, int split_number);

/* Print all block in the given table. The blocks are printed in BED format. start is 0-based and end is 1-based
 *
 * @param blocks_per_contig     stHash table of blocks (each value is a stList of blocks)
 */
void ptBlock_print_blocks_stHash(stHash* blocks_per_contig);

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
 * Merge blocks (should be sorted by stList_sort)
 * Since merging can be performed by either rfs/rfe or rds_f/rde_f this function
 * takes two functions as inputs for getting start and end of each block
 * the functions that can be used for this purpose are defined above (starting with ptBlock_get_)
 * This enables merging block by different criteria
 * it also takes set_end/set_start function for updating the end and start locations of merged block if it had
 * overlap with any new block during the merging process
 * This function is different from ptBlock_merge_blocks (without "_v2"). The current one
 *  will create a separate block wherever we have different overlapping blocks however the
 *  other one will create one block for all consecutive overlapping blocks. "_v2" is useful
 *  if we want to know the number of overlapping blocks in each merged block.
 *
 *
 * @param blocks  stList of blocks each of which saved as a ptBlock
 * @param  get_start   Function for getting the start coordinate of each block
 * @param  get_end   Function for getting the end coordinate of each block
 * @param  set_start  Function for updating the start coordinate of each block
 * @param  set_end  Function for updating the end coordinate of each block
 *
 * @return  merged blocks
 *
 * @Note    input blocks should be sorted based on the start coordinates so the sorting criteria
 *          should match the merging criteria
 *          For example if blocks are sorted by ptBlock_cmp_rfs you can merge by ptBlock_get_rfs/rfe
 *          and if they are sorted by ptBlock_cmp_rds_f you can merge by ptBlock_get_rds_f/rde_f
 *
 */

stList *ptBlock_merge_blocks_v2(stList *blocks,
                                int (*get_start)(ptBlock *), int (*get_end)(ptBlock *),
                                void (*set_start)(ptBlock *, int), void (*set_end)(ptBlock *, int));

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

/**
 * Similar to ptBlock_merge_blocks_per_contig
 * but it works with ptBlock_merge_blocks_v2() instead of ptBlock_merge_blocks()
 *
 */
stHash *ptBlock_merge_blocks_per_contig_v2(stHash *blocks_per_contig,
                                        int (*get_start)(ptBlock *), int (*get_end)(ptBlock *),
                                        void (*set_start)(ptBlock *, int), void (*set_end)(ptBlock *, int));

stHash *ptBlock_merge_blocks_per_contig_by_rf_v2(stHash *blocks_per_contig);

stHash *ptBlock_merge_blocks_per_contig_by_rd_f_v2(stHash *blocks_per_contig);

stHash *ptBlock_merge_blocks_per_contig_by_sq_v2(stHash *blocks_per_contig);



void ptBlock_add_alignment(stHash *blocks_per_contig, ptAlignment *alignment, bool init_count_data);

void ptBlock_save_in_bed(stHash *blocks_per_contig, char* bed_path, bool print_count_data);

int ptBlock_get_total_number(stHash *blocks_per_contig);

void ptBlock_add_blocks_by_contig(stHash *blocks_per_contig, char* contig, stList* blocks_to_add);

stList* ptBlock_copy_stList(stList* blocks);

#endif /* PT_BLOCK_H */

