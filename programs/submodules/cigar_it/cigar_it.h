#ifndef CIGAR_IT_H
#define CIGAR_IT_H

#include <stdbool.h>
#include "sam.h"
#include <regex.h>

// One match per cs/MD operation. Only the whole match is used, never the sub-groups.
// '+', '*' and '^' are escaped outside brackets so every libc (glibc, musl, macOS) reads them the same way.
#define CS_PATTERN "(:[0-9]+)|((\\+|-)[a-z]+)|((\\*[a-z]+)+)"
#define MD_PATTERN "([A-Z](0[A-Z])*)|([0-9]+)|(\\^[A-Z]+)"

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
        const regex_t *preg_cs; // shared, compiled once
	// md tag attributes
	bool use_md;
	char* md;
	int offset_md;
	const regex_t *preg_md; // shared, compiled once
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

int ptCigarIt_next_cs(ptCigarIt* cigar_it);

int ptCigarIt_next(ptCigarIt* cigar_it);

#endif

