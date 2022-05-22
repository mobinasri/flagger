#include <getopt.h>
#include "sam.h"
#include <assert.h>
#include "faidx.h"
#include <time.h>
#include "bgzf.h"
#include <regex.h>
#include "sonLib.h"
#include "cigar_it.h"
#include "thread_pool.h"

typedef struct {
	char	contig_name[50];
        int32_t start;
	uint8_t mapq;
}Location;

Location* location_construct(char* contig_name, int start, uint8_t mapq){
	Location* loc = malloc(sizeof(Location));
	strcpy(loc->contig_name, contig_name);
	loc->start = start;
	loc->mapq = mapq;
	return loc;
}

void location_destruct(void* location){
	Location* loc = location;
	free(loc);
}

stHash* get_phased_read_table(char* phased_reads_path){
	stHash* phased_read_table = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, location_destruct);
	if(phased_reads_path == NULL) return phased_read_table;
	FILE* fp = fopen(phased_reads_path, "r");
	size_t read;
	size_t len;
	char* line = NULL;
	char* read_name;
	char contig_name_new[50];
	char contig_name_old[50];
	int start_new;
	int start_old;
	char* token;
	Location* loc;
	while((read = getline(&line, &len, fp)) != -1){
		line[strlen(line)-1] = '\0';
		if(line[0] == '$'){
			token = strtok(line, "\t");
			token = strtok(NULL, "\t");
			read_name = malloc(strlen(token) + 1);
			strcpy(read_name, token);
			start_new = -1;
			start_old = -1;
		}
		else if (line[0] == '@'){
			token = strtok(line, "\t");
                        token = strtok(NULL, "\t");
			token = strtok(NULL, "\t");
			strcpy(contig_name_new, token);
			token = strtok(NULL, "\t");
			start_new = atoi(token); // 0-based
			loc = location_construct(contig_name_new, start_new, -1);
			// if pri (old) and sec (new) start exactly at the same place ignore it
			if (start_old != -1 
			    && (start_old != start_new 
			        || strcmp(contig_name_old, contig_name_new) != 0)){
				stHash_insert(phased_read_table, read_name, loc);
			}
		}
		else if(line[0] == '*'){
			token = strtok(line, "\t");
                        token = strtok(NULL, "\t");
                        token = strtok(NULL, "\t");
                        strcpy(contig_name_old, token);
                        token = strtok(NULL, "\t");
                        start_old = atoi(token); // 0-based
			// if pri (old) and sec (new) start exactly at the same place ignore it
                        if (start_new != -1 
                            && (start_old != start_new
                                || strcmp(contig_name_old, contig_name_new) != 0)){
                                stHash_insert(phased_read_table, read_name, loc);
                        }
		}

	}
	return phased_read_table;
}

bool is_prim(stHash* phased_read_table, bam1_t* b, sam_hdr_t* sam_hdr){
	int start = b->core.pos; // 0-based
	const char* contig_name = sam_hdr_tid2name(sam_hdr, b->core.tid);
	char* read_name = bam_get_qname(b);
	Location* loc = stHash_search(phased_read_table, read_name);
	if(loc != NULL){
		if(strcmp(loc->contig_name, contig_name) == 0 && 
		   loc->start == start){
			return true;
			printf("%s\t%s\t%d\n", read_name, loc->contig_name, loc->start);
		}
		else{
			return false;
		}
	}
	return (b->core.flag & BAM_FSECONDARY) == 0; 
}

stHash* get_mapq_table(char* mapq_table_path){
	stHash* mapq_table = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, (void (*)(void *)) stList_destruct);
	if (mapq_table_path == NULL) return mapq_table;
        FILE* fp = fopen(mapq_table_path, "r");
        size_t read;
        size_t len;
        char* line = NULL;
        char* read_name;
        char* contig_name;
        int start;
	int mapq;
        char* token;
        Location* loc;
	stList* locs;
        while((read = getline(&line, &len, fp)) != -1){
                line[strlen(line)-1] = '\0';
                token = strtok(line, "\t");
                read_name = malloc(strlen(token) + 1);
               	strcpy(read_name, token);
                token = strtok(NULL, "\t");
                contig_name = malloc(strlen(token) + 1);
                strcpy(contig_name, token);
                token = strtok(NULL, "\t");
                start = atoi(token) - 1; // 0-based
		token = strtok(NULL, "\t");
                mapq = atoi(token);
		//printf("%d\n",mapq);
                loc = location_construct(contig_name, start, mapq);
		locs = stHash_search(mapq_table, read_name);
		if (locs == NULL){
			locs = stList_construct3(0, location_destruct);
			stHash_insert(mapq_table, read_name, locs);
		}
		stList_append(locs, loc);
        }
        return mapq_table;
}


stSet* get_read_set(char* exclude_path){
        stSet* exclude_read_set = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free);
        if (exclude_path == NULL) return exclude_read_set;
        FILE* fp = fopen(exclude_path, "r");
        size_t read;
        size_t len;
        char* line = NULL;
        char* read_name;
        while((read = getline(&line, &len, fp)) != -1){
                line[strlen(line)-1] = '\0';
                read_name = malloc(strlen(line));
		strcpy(read_name, line);
                stSet_insert(exclude_read_set, read_name);
        }
        return exclude_read_set;
}

uint8_t get_mapq(stHash* mapq_table, bam1_t* b, sam_hdr_t* sam_hdr){
        int start = b->core.pos; // 0-based
        const char* contig_name = sam_hdr_tid2name(sam_hdr, b->core.tid);
        char* read_name = bam_get_qname(b);
       	stList* locs = stHash_search(mapq_table, read_name);
	Location* loc;
        if(locs != NULL){
		for(int i = 0; i < stList_length(locs); i++){
			loc = stList_get(locs, i);
			if(strcmp(loc->contig_name, contig_name) == 0 &&
			   loc->start == start){
				//printf("%s\t%d\n",loc->contig_name, loc->mapq);
				return loc->mapq;
			}
		}
        }
        return b->core.qual;
}

int get_read_length(bam1_t* b){
	ptCigarIt* cigar_it = ptCigarIt_construct(b, true, true);
	int read_len = 0;
	while(ptCigarIt_next(cigar_it)){
		if(cigar_it->op == BAM_CDIFF ||
		   cigar_it->op == BAM_CMATCH ||
		   cigar_it->op == BAM_CEQUAL ||
		   cigar_it->op == BAM_CINS ||
		   cigar_it->op == BAM_CHARD_CLIP  ||
		   cigar_it->op == BAM_CSOFT_CLIP){
			read_len += cigar_it->len;
		}
	}
	ptCigarIt_destruct(cigar_it);
	return read_len;
}

int get_alignment_length(bam1_t* b){
        ptCigarIt* cigar_it = ptCigarIt_construct(b, true, true);
        int alignment_len = 0;
        while(ptCigarIt_next(cigar_it)){
                if(cigar_it->op == BAM_CDIFF ||
                   cigar_it->op == BAM_CMATCH ||
                   cigar_it->op == BAM_CEQUAL) {
                        alignment_len += cigar_it->len;
                }
        }
	ptCigarIt_destruct(cigar_it);
        return alignment_len;
}


static struct option long_options[] =
{
    {"inputBam", required_argument, NULL, 'i'},
    {"outputBam", required_argument, NULL, 'o'},
    {"phasingLog", required_argument, NULL, 'P'},
    {"mapqTable", required_argument, NULL, 'M'},
    {"minReadLen", required_argument, NULL, 'm'},
    {"minAlignmentLen", required_argument, NULL, 'a'},
    {"primaryOnly", 0, NULL, 'p'},
    {"exclude", required_argument, NULL, 'e'},
    {"threads", required_argument, NULL, 'n'},
    {"noTag", 0, NULL, 't'},
    {"maxMapq", required_argument, NULL, 'x'},
    {"maxDiv", required_argument, NULL, 'd'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char *argv[]){
        int c;
        char* input_path;
	char* output_path;
	char* exclude_path=NULL;
	char* phasing_log_path=NULL;
	char* mapq_table_path=NULL;
	bool primary_only = false;
	bool no_tag = false;
	int min_read_length = 5000;
	int min_alignment_length = 5000;
	double max_divergence = 0.12;
	int nthreads = 2;
	int max_mapq = 60;
	char *program;
        (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
        while (~(c=getopt_long(argc, argv, "i:o:x:e:P:M:tpm:a:n:d:h", long_options, NULL))) {
                switch (c) {
                        case 'i':
                                input_path = optarg;
                                break;
                        case 'o':
                                output_path = optarg;
                                break;
			case 'x':
                                max_mapq = atoi(optarg);
                                break;
			case 'e':
				exclude_path = optarg;
				break;
			case 'P':
                                phasing_log_path = optarg;
                                break;
			case 'M':
				mapq_table_path = optarg;
				break;
			case 't':
				no_tag = true;
				break;
			case 'p':
				primary_only = true;
				break;
			case 'm':
                                min_read_length = atoi(optarg);
                                break;
			case 'a':
                                min_alignment_length = atoi(optarg);
                                break;
			case 'n':
				nthreads = atoi(optarg);
				break;
			case 'd':
				max_divergence = atof(optarg);
				break;
                        default:
                                if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
                        help:
                                fprintf(stderr, "\nUsage: %s  -i <INPUT_BAM> -o <OUTPUT_BAM> -p <PHASING_LOG> -m <MAPQ_TABLE>\n\tModify the input bam file:\n\t* Apply the phasing log by swapping the primary and secondary alignments whenever necessary(stdout log of ./phase_reads)\n\t* Set the MAPQs to the values given in the mapq table\n\t\tmapq table is a tab delimited text containing 4 columns:\n\t\t1. read name\n\t\t2. contig name\n\t\t3. left-most coordinate on contig (1-based)\n\t\t4. adjusted mapq\n\t* Filter secondary alignments (After applying the phasing log)\n\t* Skip outputing the optional fields (like cs and MD tags)\n\t* Filter the reads shorter than the given threshold\n\t* Filter the alignments shorter than the given threshold\n\n", program);
                                fprintf(stderr, "Options:\n");
                                fprintf(stderr, "         --inputBam,\t-i         input bam file\n");
                                fprintf(stderr, "         --outputBam,\t-o         output bam file\n");
				fprintf(stderr, "         --maxMapq,\t-x         maximum mapq [default:60]\n");
				fprintf(stderr, "         --phasingLog,\t-P         the phasing log path [optional]\n");
				fprintf(stderr, "         --mapqTable,\t-M         the adjusted mapq table path [optional]\n");
				fprintf(stderr, "         --exclude,\t-e         Path to a file containing the read names that have to be excluded [optional]\n");
                                fprintf(stderr, "         --noTag,\t-t         output no optional fields\n");
                                fprintf(stderr, "         --primaryOnly,\t-p         output only primary alignments\n");
				fprintf(stderr, "         --minReadLen,\t-m         min read length [default: 5k]\n");
				fprintf(stderr, "         --minAlignmentLen,\t-a         min alignment length [default: 5k]\n");
				fprintf(stderr, "         --maxDiv,\t-d         min gap-compressed divergence (\"de\" tag) [default: 0.12]\n");
				fprintf(stderr, "         --threads,\t-n         number of threads (for bam I/O)[default: 2]\n");
				return 1;
                }
        }
	//open input bam file and read the header
        samFile* fp = sam_open(input_path, "r");
        sam_hdr_t* sam_hdr = sam_hdr_read(fp);
	//open output bam file and write the header
	samFile* fo = sam_open(output_path, "wb");
	sam_hdr_write(fo, sam_hdr);
	// Make a multi threading pool
        htsThreadPool p = {NULL, 0};
        if (nthreads > 0) {
                p.pool = hts_tpool_init(nthreads);
                if (!p.pool) {
                        fprintf(stderr, "Error creating thread pool\n");
                }
                else{ // Add I/O streams to the threading pool
                        hts_set_opt(fp, HTS_OPT_THREAD_POOL, &p);
                        hts_set_opt(fo, HTS_OPT_THREAD_POOL, &p);
                }
        }
	//get the table of phased reads (the locations of their primary alignments are stored as values)
	stHash* phased_read_table = get_phased_read_table(phasing_log_path);
	//get the table of reads with adjusted mapq (the locations (+ corresponding mapqs) are stored as values)
	//stHashIterator* it = stHash_getIterator(phased_read_table);
	char* read_name;
	//Location* loc;
	stHash* mapq_table = get_mapq_table(mapq_table_path);
	stSet* exclude_read_set = get_read_set(exclude_path);
	/*stHashIterator* it = stHash_getIterator(mapq_table);
	stList* locs;
	Location* loc;
	while((read_name = stHash_getNext(it)) != NULL){
		locs = stHash_search(mapq_table, read_name);
		for(int i=0; i< stList_length(locs); i++){
			loc = stList_get(locs, i);
			printf("%s\t%s\t%d\t%d\n", read_name, loc->contig_name, loc->start, loc->mapq);
		}
	}*/
	//modify mapq and sec/pri tags based on given logs and tables
	bam1_t* b = bam_init1();
	double divergence = 0.0;
	while(sam_read1(fp, sam_hdr, b) > -1) {
		if (b->core.flag & BAM_FUNMAP) continue; // if unmapped
		if (stSet_search(exclude_read_set, bam_get_qname(b)) != NULL) continue; //if read should be excluded
		if (is_prim(phased_read_table, b, sam_hdr)){
			b->core.flag &= ~BAM_FSECONDARY; // make it primary
		}
		else{
			if(primary_only) continue;
			b->core.flag |= BAM_FSECONDARY; // make it secondary
		}
		// skip short reads or short alignments
		if (get_read_length(b) < min_read_length || 
		    get_alignment_length(b) < min_alignment_length){
			continue;
		}
		b->core.qual = get_mapq(mapq_table, b, sam_hdr);
		if(max_mapq < b->core.qual) continue;
		divergence = bam_aux2f(bam_aux_get(b, "de"));
		if( max_divergence < divergence) continue;
		if(no_tag) b->l_data -= bam_get_l_aux(b);
		if (sam_write1(fo, sam_hdr, b) == -1) {
			fprintf(stderr, "Couldn't write %s\n", bam_get_qname(b));
                }
	}
	sam_close(fo);
	sam_close(fp);
        bam_destroy1(b);
        sam_hdr_destroy(sam_hdr);
	stHash_destruct(mapq_table);
	stHash_destruct(phased_read_table);
	if (p.pool)
                hts_tpool_destroy(p.pool);
}
