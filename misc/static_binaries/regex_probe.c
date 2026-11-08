// Temporary diagnostic: finds cs/MD tags that the libc regex splits differently from the
// expected tokens, to see why bam2cov output differs on macOS.
//
// For every alignment in a BAM it walks the cs (or MD) tag the same way cigar_it.c does
// (regexec with REG_NOTBOL, then continue right after the match) and compares each match with
// a plain scanner that finds the expected token. On the first difference in a tag it reports
// the regexec return code, the offset, the tag length and the text around that point, and
// retries regexec on a short copy of the remaining text to tell apart length-dependent failures
// from pattern problems.
//
// The tags are checked once from one thread and once from 8 threads sharing the compiled
// patterns, as bam2cov does. Findings are printed as GitHub annotations.
//
//   regex_probe <bam> <cs|md>
#include <pthread.h>
#include <regex.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sam.h"
#include "cigar_it.h"

#define MAX_REPORTS 4
#define N_THREADS 8
#define SHORT_COPY 200

static regex_t preg;
static bool use_cs;
static char** tags;
static char** names;
static int n_tags;

static bool is_digit(char c) { return '0' <= c && c <= '9'; }
static bool is_lower(char c) { return 'a' <= c && c <= 'z'; }
static bool is_upper(char c) { return 'A' <= c && c <= 'Z'; }

// expected next token at or after s[0]: same tokens as CS_PATTERN / MD_PATTERN
static bool find_token(const char* s, int* start, int* end) {
    for (int i = 0; s[i] != '\0'; i++) {
        int j = i + 1;
        if (use_cs) {
            if (s[i] == ':' && is_digit(s[j])) {
                while (is_digit(s[j])) j++;
            } else if ((s[i] == '+' || s[i] == '-') && is_lower(s[j])) {
                while (is_lower(s[j])) j++;
            } else if (s[i] == '*' && is_lower(s[j])) {
                while (is_lower(s[j])) j++;
                while (s[j] == '*' && is_lower(s[j + 1])) {
                    j++;
                    while (is_lower(s[j])) j++;
                }
            } else continue;
        } else {
            if (is_upper(s[i])) {
                while (s[j] == '0' && is_upper(s[j + 1])) j += 2;
            } else if (is_digit(s[i])) {
                while (is_digit(s[j])) j++;
            } else if (s[i] == '^' && is_upper(s[j])) {
                while (is_upper(s[j])) j++;
            } else continue;
        }
        *start = i;
        *end = j;
        return true;
    }
    return false;
}

typedef struct {
    int first, step;  // tags first, first+step, ...
    int n_bad;
    int n_reports;
    char reports[MAX_REPORTS][1024];
} job_t;

static void check_tag(job_t* job, int t) {
    const char* tag = tags[t];
    int len = strlen(tag);
    int off = 0;
    while (1) {
        regmatch_t pm;
        int s = 0, e = 0;
        int rc = regexec(&preg, tag + off, 1, &pm, REG_NOTBOL);
        bool want = find_token(tag + off, &s, &e);
        bool same = rc == 0 ? (want && pm.rm_so == s && pm.rm_eo == e) : (rc == REG_NOMATCH && !want);
        if (!same) {
            job->n_bad++;
            if (job->n_reports < MAX_REPORTS) {
                char err[128];
                if (rc != 0) regerror(rc, &preg, err, sizeof(err));
                else snprintf(err, sizeof(err), "matched [%d,%d), expected [%d,%d)",
                              (int) pm.rm_so, (int) pm.rm_eo, s, e);
                // retry on a short copy of the remaining text
                char copy[SHORT_COPY + 1];
                int n = len - off < SHORT_COPY ? len - off : SHORT_COPY;
                memcpy(copy, tag + off, n);
                copy[n] = '\0';
                regmatch_t pm2;
                int rc2 = regexec(&preg, copy, 1, &pm2, REG_NOTBOL);
                int lo = off < 40 ? 0 : off - 40;
                int hi = off + 60 < len ? off + 60 : len;
                snprintf(job->reports[job->n_reports++], 1024,
                         "read %s: regexec rc=%d (%s) at offset %d of %d; "
                         "short copy of the rest gives rc=%d%s; text: ...%.*s|%.*s...",
                         names[t], rc, err, off, len, rc2,
                         rc2 == 0 && want && pm2.rm_so == s && pm2.rm_eo == e ? " (correct)" : "",
                         off - lo, tag + lo, hi - off, tag + off);
            }
            return;
        }
        if (rc != 0) return;  // reached the end of the tag, as expected
        off += pm.rm_eo;
    }
}

static void* run_job(void* arg) {
    job_t* job = arg;
    for (int t = job->first; t < n_tags; t += job->step) check_tag(job, t);
    return NULL;
}

static void run(int n_threads) {
    job_t* jobs = calloc(n_threads, sizeof(job_t));
    pthread_t* th = malloc(n_threads * sizeof(pthread_t));
    for (int i = 0; i < n_threads; i++) {
        jobs[i].first = i;
        jobs[i].step = n_threads;
        pthread_create(&th[i], NULL, run_job, &jobs[i]);
    }
    int n_bad = 0, shown = 0;
    for (int i = 0; i < n_threads; i++) {
        pthread_join(th[i], NULL);
        n_bad += jobs[i].n_bad;
    }
    printf("::notice::regex_probe %s, %d thread(s): %d of %d tags split differently from the expected tokens\n",
           use_cs ? "cs" : "MD", n_threads, n_bad, n_tags);
    for (int i = 0; i < n_threads && shown < MAX_REPORTS; i++)
        for (int r = 0; r < jobs[i].n_reports && shown < MAX_REPORTS; r++, shown++)
            printf("::warning::regex_probe %s, %d thread(s): %s\n", use_cs ? "cs" : "MD", n_threads,
                   jobs[i].reports[r]);
    free(jobs);
    free(th);
}

int main(int argc, char** argv) {
    if (argc != 3) {
        fprintf(stderr, "usage: regex_probe <bam> <cs|md>\n");
        return 1;
    }
    use_cs = strcmp(argv[2], "cs") == 0;
    int rc = regcomp(&preg, use_cs ? CS_PATTERN : MD_PATTERN, REG_EXTENDED);
    if (rc != 0) {
        char err[128];
        regerror(rc, &preg, err, sizeof(err));
        printf("::warning::regex_probe: regcomp failed with %d (%s)\n", rc, err);
        return 0;
    }
    samFile* fp = sam_open(argv[1], "r");
    bam_hdr_t* hdr = fp ? sam_hdr_read(fp) : NULL;
    if (!hdr) {
        printf("::warning::regex_probe: cannot read %s\n", argv[1]);
        return 0;
    }
    bam1_t* b = bam_init1();
    int cap = 1024;
    tags = malloc(cap * sizeof(char*));
    names = malloc(cap * sizeof(char*));
    while (sam_read1(fp, hdr, b) >= 0) {
        uint8_t* aux = bam_aux_get(b, use_cs ? "cs" : "MD");
        if (!aux) continue;
        if (n_tags == cap) {
            cap *= 2;
            tags = realloc(tags, cap * sizeof(char*));
            names = realloc(names, cap * sizeof(char*));
        }
        tags[n_tags] = strdup(bam_aux2Z(aux));
        names[n_tags] = strdup(bam_get_qname(b));
        n_tags++;
    }
    int max_len = 0;
    for (int t = 0; t < n_tags; t++) {
        int l = strlen(tags[t]);
        if (l > max_len) max_len = l;
    }
    printf("read %d %s tags from %s (longest %d characters)\n", n_tags, use_cs ? "cs" : "MD", argv[1], max_len);
    run(1);
    run(N_THREADS);
    return 0;
}
