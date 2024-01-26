#ifndef PT_BLOCK_H
#define PT_BLOCK_H

#include <assert.h>
#include <math.h>
#include <float.h>
#include "common.h"
#include "sonLib.h"
#include "stdlib.h"
#include "stdio.h"
#include "ptAlignment.h"
#include <zlib.h>
#include "tpool.h"
#include "ptBlock.h"
#include "cJSON.h"

/*! @typedef
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

/*! @typedef
 * @abstract Structure for keeping useful information about a block with the same coverage/annotation
 * @field annotation_flag     a 32-bit flag where each bit represents a single annotation. Therefore it can only
 *                           keep track of 32 distinct annotations. For example "0...00000100" means the 3rd annotation
 *                           or "0...00011100" means that this block is completely within three different annotations;
 *                           3rd, 4th and 5th
 * @field coverage              The total read depth of coverage in this block
 * @field coverage_high_mapq    The read depth of coverage for only the alignments with high mapqs (for example >20)
 * @field coverage_high_clip    The read depth of coverage for only the alignments that are highly clipped (for example >10%)
 */
typedef struct CoverageInfo {
    int32_t annotation_flag;
    u_int8_t coverage;
    u_int8_t coverage_high_mapq;
    u_int8_t coverage_high_clip;
} CoverageInfo;

/*! @typedef
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

/**
 * Set data and the related functions in a ptBlock structure.
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

// the function for representing the count data as a string
char *get_string_count_data(void* src_);

/**
 * Creates an instance of CoverageInfo given the required attributes
 *
 * @param annotation_flag       32-bit flag for representing at most 32 different annotations
 * @param coverage              coverage value for the related block
 * @param coverage_high_mapq    coverage of the alignments with high mapq for the related block
 * @param coverage_high_clip    coverage of the highly clipped alignments for the related block
 */
CoverageInfo *CoverageInfo_construct(int32_t annotation_flag,
                            u_int8_t coverage,
                            u_int8_t coverage_high_mapq,
                            u_int8_t coverage_high_clip);

CoverageInfo *CoverageInfo_copy(CoverageInfo *coverageInfo);
CoverageInfo **CoverageInfo_copy1DArray(CoverageInfo **coverageInfo, int len);
void CoverageInfo_destruct1DArray(CoverageInfo **coverageInfo1DArray, int len);
void CoverageInfo_destruct(CoverageInfo *coverageInfo);


/**
 * Receives an alignment, creates a CoverageInfo struct based on the given thresholds on mapq and clipping ratio
 *
 *
 * @param alignment             the alignment
 * @param min_mapq              the minimum mapq value for determining if the alignment is of high mapq or not
 * @param min_clipping_ratio    the minimum clipping ratio for determining if the alignment is highly clipped or not
 *
 */
CoverageInfo *CoverageInfo_construct_from_alignment(ptAlignment *alignment, int min_mapq, double min_clipping_ratio);

// the three functions for keeping CoverageInfo data in ptBlock structs

void extend_cov_info_data(void *dest_, void *src_);


void destruct_cov_info_data(void* src);


void *copy_cov_info_data(void* src_);

// the functions for converting coverage info into a single string

// format 1 is useful for debugging purposes
char *get_string_cov_info_data_format_1(void* src_);

// format 2 is compatible for HMM-Flagger
char *get_string_cov_info_data_format_2(void* src_);

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
 * All blocks in each contig/chromosome are saved in a single list as the corresponding value
 *
 * @param bed_path  Path to a BED file
 * @return  stHash table with contigs as keys and tracks as values
 *
 */
stHash *ptBlock_parse_bed(char *bed_path);


/**
 * Compare two ptBlock structures
 *
 * @param block_1       the first block
 * @param block_2       the second block
 * @return              true if the ref/read coordinates of the two blocks were exactly the same
 */
bool ptBlock_is_equal(ptBlock* block_1, ptBlock* block_2);


/**
 * Compare two lists of ptBlock structures
 *
 * @param blocks_1      the first list of blocks
 * @param blocks_2      the second list of blocks
 * @return              true if the ref/read coordinates of the two block lists were exactly the same
 */
bool ptBlock_is_equal_stList(stList* blocks_1, stList* blocks_2);


/**
 * Compare two stHash tables of blocks
 *
 * @param blocks_1      the first table of blocks
 * @param blocks_2      the second table of blocks
 * @return              true if the contig names and ref/read coordinates of the blocks in the two tables were exactly the same
 */
bool ptBlock_is_equal_stHash(stHash* blocks_per_contig_1, stHash* blocks_per_contig_2);


// Functions for block iterator


/**
 * Construct a ptBlockItrPerContig structure
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



/**
 * Return the next block and update the contig name
 *
 * @param block_iter        the block iterator
 * @param ctg_name          contig name (to be able to update in place)
 * @return block            ptBlock of the next block
 */
ptBlock* ptBlockItrPerContig_next(ptBlockItrPerContig *block_iter, char* ctg_name);


/**
 * Destruct a ptBlockItrPerContig structure
 *
 * This function destructs the block iterator. It will NOT
 * destruct the table of blocks that was used for creating the iterator
 */
void ptBlockItrPerContig_destruct(ptBlockItrPerContig *blockItr);


/**
 * Add a block to the blocks table in which each key is the corresponding contig name
 * and each value is a stList of blocks
 *
 * @param blocks_per_contig     stHash table of blocks (each value is a stList of blocks)
 * @param block                 the block to add
 * @param ctg_name              the contig name
 */
void ptBlock_add_block_to_stList_table(stHash* blocks_per_contig, ptBlock* block, char* ctg_name);


/**
 * This function recieves a stHash table of blocks and returns a stList of smaller tables. Each
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

/**
 * Print all block in the given table. The blocks are printed in BED format. start is 0-based and end is 1-based
 *
 * @param blocks_per_contig     stHash table of blocks (each value is a stList of blocks)
 * @param get_string_function   a function that takes the additional data ("data" attribute) saved in
 *                              the ptBlock struct and converts it to a string.
 *                              If it is set to NULL then only coordinates
 *                              will be printed with no additional data.
 * @param fp                    opened file to write the blocks in (can be also stdout/stderr)
 *                              it can also be a POINTER to the output of gzopen() for the compressed mode
 * @param is_compressed		    true if fp is of type gzFile*, false if stderr/stdout/FILE*
 */
void ptBlock_print_blocks_stHash_in_bed(stHash* blocks_per_contig,
                                        char * (*get_string_function)(void *),
                                        void* fp,
                                        bool is_compressed);

/**
 * Print all block in the given table. The blocks are printed in COV format. start is 1-based and end is 1-based
 *
 * @param blocks_per_contig     stHash table of blocks (each value is a stList of blocks)
 * @param get_string_function   a function that takes the additional data ("data" attribute) saved in
 *                              the ptBlock struct and converts it to a string.
 *                              If it is set to NULL then only coordinates
 *                              will be printed with no additional data.
 * @param fp                    opened file to write the blocks in (can be also stdout/stderr)
 *                              it can also be a POINTER to the output of gzopen() for the compressed mode
 * @param is_compressed		    true if fp is of type gzFile*, false if stderr/stdout/FILE*
 * @param ctg_to_len            stHash table to convert contig name to contig length (each value is of type int*)
 */
void ptBlock_print_blocks_stHash_in_cov(stHash* blocks_per_contig,
                                        char * (*get_string_function)(void *),
                                        void* fp,
                                        bool is_compressed,
                                        stHash* ctg_to_len);

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


/**
 * Receives an alignment, creates a ptBlock struct based on the reference start and end
 * coordinates and add this block to the given block table
 *
 * @param blocks_per_contig     stHash table of blocks (each value is a stList of blocks)
 * @param  alignment            the alignment
 * @param  init_count_data      if true then create a count data for new ptBlock and set its value to 1
 *
 */
void ptBlock_add_alignment(stHash *blocks_per_contig,
                           ptAlignment *alignment,
                           bool init_count_data);

/**
 * Receives an alignment, creates a ptBlock struct based on the reference start and end
 * coordinates and add this block to the given block table. The new ptBlock will have a
 * CoverageInfo struct as its data. The content of the data is determined based on "min_mapq"
 * and "min_clipping_ratio"
 *
 * @param blocks_per_contig     stHash table of blocks (each value is a stList of blocks)
 * @param alignment             the alignment
 * @param min_mapq              the minimum mapq value for determining if the alignment is of high mapq or not
 * @param min_clipping_ratio    the minimum clipping ratio for determining if the alignment is highly clipped or not
 *
 */
void ptBlock_add_alignment_as_CoverageInfo(stHash *blocks_per_contig,
                                           ptAlignment *alignment,
                                           int min_mapq,
                                           double min_clipping_ratio);

/**
 * Receives a table of blocks and adds the given data to all blocks. Note that the given data
 * will be copied once per block before adding.
 *
 *
 * @param blocks_per_contig     stHash table of blocks (each value is a stList of blocks)
 * @param free_data             function for freeing the memory allocated to the data.
 * @param copy_data             function for copying data
 * @param append_data           function for appending new data to what is already saved in the block
 *
 */
void ptBlock_add_data_to_all_blocks_stHash(stHash *blocks_per_contig,
                                           void* data,
                                           void (*destruct_data)(void *),
                                           void *(*copy_data)(void *),
                                           void (*extend_data)(void *, void *));


/**
 * copy and add all blocks from the source block table to the destination block table
 *
 *
 * @param blocks_per_contig_dest    the destination block table (stHash table with block lists as values)
 * @param blocks_per_contig_src     the source block table (stHash table with block lists as values)
 *
 */
void ptBlock_extend_block_tables(stHash *blocks_per_contig_dest, stHash *blocks_per_contig_src);



void ptBlock_save_in_bed(stHash *blocks_per_contig, char* bed_path, bool print_count_data);

int64_t ptBlock_get_total_number(stHash *blocks_per_contig);

int64_t ptBlock_get_total_length(stHash *blocks_per_contig, int (*get_start)(ptBlock *), int (*get_end)(ptBlock *));

int64_t ptBlock_get_total_length_by_rf(stHash *blocks_per_contig);

int64_t ptBlock_get_total_length_by_rd_f(stHash *blocks_per_contig);

int64_t ptBlock_get_total_length_by_sq(stHash *blocks_per_contig);

void ptBlock_add_blocks_by_contig(stHash *blocks_per_contig, char* contig, stList* blocks_to_add);

stList* ptBlock_copy_stList(stList* blocks);

// takes the path to a bam file
// returns a stHash table from contig name to length
stHash* ptBlock_get_contig_length_stHash_from_bam(char* bam_path);


// functions and structs for extracting coverage information from bam/sam file


typedef struct ArgumentsCovExt {
    stHash* coverage_blocks_per_contig;
    stHash* ref_blocks_per_contig_to_parse;
    char* bam_path;
    pthread_mutex_t *mutexPtr;
    int min_mapq;
    double min_clipping_ratio;
} ArgumentsCovExt;

// parse bam file and create a stHash table of blocks
// coverage information is saved in a CoverageInfo object
// available as the "data" attribute of each resulting block
stHash* ptBlock_multi_threaded_coverage_extraction(char* bam_path,
                                                   int threads,
                                                   int min_mapq,
                                                   double min_clipping_ratio);

// make a block table that covers the whole reference sequences
stHash* ptBlock_get_whole_genome_blocks_per_contig(char* bam_path);

int get_annotation_index(stList* annotation_names, char* annotation_name);
stList *parse_annotation_names_and_save_in_stList(char* json_path);

// parse bam file and create a stHash table of blocks
// this function calls "ptBlock_multi_threaded_coverage_extraction"
// and it adds the blocks with 0 coverage and also fills the annotation
// fields based on the bed files whose paths are given in a json file
stHash* ptBlock_multi_threaded_coverage_extraction_with_zero_coverage_and_annotation(char* bam_path,
                                                                                     char* json_path,
                                                                                     int threads,
                                                                                     int min_mapq,
                                                                                     double min_clipping_ratio);
#endif /* PT_BLOCK_H */

