#ifndef CIGAR_IT_H
#define CIGAR_IT_H

#include <stdbool.h>
#include "sam.h"

// cs and MD tags are split into operations by ptCigarIt_find_cs_token() and
// ptCigarIt_find_md_token(), which find the same tokens as these POSIX regexes
// (leftmost-longest match, searching forward from the current offset):
//   cs: (:[0-9]+)|([+-][a-z]+)|((\*[a-z]+)+)
//   MD: ([A-Z](0[A-Z])*)|([0-9]+)|(\^[A-Z]+)
// They are not matched with regcomp()/regexec() because the macOS regex library
// rejected or mis-matched valid tags (e.g. "regcomp() failed, returning nonzero (11)").

typedef struct {
        // the constant attributes
        uint32_t* cigar;
        int n; // number of operations
        bool is_rev;
        // the attributes below may change in each iteration
        uint8_t op; // current cigar operation
        int len; // length of the current operation
        int idx;
        // cs tag attributes
        bool use_cs;
        char* cs;
        int offset_cs;
	// md tag attributes
	bool use_md;
	char* md;
	int offset_md;
	// private attributes
	// This is used to sync the iteration over mis/matches in CIGAR and MD (or CS) tag
	int match_remain;
        //location attributes
        int rds_f; // read start forward
        int rde_f; // read end forward
        int sqs; // seq start
        int sqe; // seq end
        int rfs; // ref start
        int rfe; // ref end
}ptCigarIt;

void ptCigarIt_destruct(ptCigarIt* cigar_it);

ptCigarIt* ptCigarIt_construct(bam1_t* b, bool use_cs, bool use_md);

// Find the first cs (or MD) token in the string s. On success it returns true and sets
// [*start, *end) to the token's position in s; it returns false if s has no more tokens.
bool ptCigarIt_find_cs_token(const char* s, int* start, int* end);

bool ptCigarIt_find_md_token(const char* s, int* start, int* end);

int ptCigarIt_next_cs(ptCigarIt* cigar_it);

int ptCigarIt_next(ptCigarIt* cigar_it);

#endif

