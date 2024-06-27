#include "ptVariant.h"
#include "common.h"

ptVariant *ptVariant_construct(char *contig, int32_t pos, int8_t type, float vaf, int8_t gq, int64_t ps) {
    ptVariant *variant = (ptVariant *) malloc(sizeof(ptVariant));
    strcpy(variant->contig, contig);
    variant->pos = pos;
    variant->type = type;
    variant->vaf = vaf;
    variant->gq = gq;
    variant->ps = ps;
    // Initialize gt and alleles variables
    // They can be filled later with ptMarker_append_gt() and ptMarker_append_allele()
    variant->gt = NULL;
    variant->gt_len = 0;
    variant->alleles = stList_construct3(0, free);
    variant->longest_allele_len = 0;
    return variant;
}


void ptVariant_destruct(ptVariant *variant) {
    free(variant->gt);
    stList_destruct(variant->alleles);
    free(variant);
}


int ptVariant_cmp(const void *a, const void *b) {
    ptVariant *v1 = (ptVariant *) a;
    ptVariant *v2 = (ptVariant *) b;
    // first compare by contig then pos
    if (strcmp(v1->contig, v2->contig) != 0) {
        return strcmp(v1->contig, v2->contig);
    } else {
        return v1->pos - v2->pos;
    }
}


void ptVariant_append_allele(ptVariant *variant, char *allele) {
    char *allele_copy = (char *) malloc((strlen(allele) + 1) * sizeof(char));
    strcpy(allele_copy, allele);
    stList_append(variant->alleles, allele_copy);

    // update longest allele len
    if (variant->longest_allele_len < strlen(allele)) {
        variant->longest_allele_len = strlen(allele);
    }
}


void ptVariant_append_gt(ptVariant *variant, int8_t gt) {
    // Increase size of gt list
    if (variant->gt == NULL) {
        variant->gt = malloc(1 * sizeof(int8_t));
    } else {
        variant->gt = realloc(variant->gt, (variant->gt_len + 1) * sizeof(int8_t));
    }
    variant->gt[variant->gt_len] = gt;
    variant->gt_len += 1;
}

ptVariant *ptVariant_copy(ptVariant *src) {
    ptVariant *dest = (ptVariant *) malloc(sizeof(ptVariant));
    strcpy(dest->contig, src->contig);
    dest->pos = src->pos;
    dest->type = src->type;
    dest->vaf = src->vaf;
    dest->gq = src->gq;
    dest->ps = src->ps;
    dest->gt = (int8_t *) malloc(src->gt_len * sizeof(int8_t));
    dest->gt_len = src->gt_len;
    memcpy(dest->gt, src->gt, src->gt_len);
    dest->alleles = stList_construct3(0, free);
    for (int i = 0; i < stList_length(src->alleles); i++) {
        char *allele = stList_get(src->alleles, i);
        char *allele_copy = (char *) malloc((strlen(allele) + 1) * sizeof(char));
        strcpy(allele_copy, allele);
        stList_append(dest->alleles, allele_copy);
    }
    dest->longest_allele_len = src->longest_allele_len;
    return dest;
}


bool ptVariant_is_equal(ptVariant *var1, ptVariant *var2) {
    return (var1->pos == var2->pos) && (strcmp(var1->contig, var2->contig) == 0);
}


bool ptVariant_exist_in_list(ptVariant *var, stList *var_list) {
    for (int i = 0; i < stList_length(var_list); i++) {
        if (ptVariant_is_equal(var, stList_get(var_list, i))) {
            return true;
        }
    }
    return false;
}

void ptVariant_print(ptVariant *variant) {
    printf("\n## Variant##\n");
    switch (variant->type) {
        case VCF_REF:
            printf("REF\n");
            break;
        case VCF_SNP:
            printf("SNP\n");
            break;
        case VCF_DEL:
            printf("DEL\n");
            break;
        case VCF_INS:
            printf("INS\n");
            break;
        case VCF_INDEL:
            printf("INDEL\n");
            break;
        case VCF_ANY:
            printf("ANY\n");
            break;
        default:
            printf("OTHER\n");
            break;
    }
    printf("contig\tpos\tvaf\tgq\tps\n");
    printf("%s\t%d\t%.2f\t%d\t%ld\n", variant->contig, variant->pos, variant->vaf, variant->gq, variant->ps);
    for (int i = 0; i < variant->gt_len; i++) {
        printf("Genotype %d:%d %s\n", i, variant->gt[i], (char *) stList_get(variant->alleles, variant->gt[i]));
        //printf("Genotype %d\n",variant->gt[i]);
    }
}


void ptVariant_swap_gt(stList *variants, int start, int end) {
    for (int i = start; i <= end; i++) {
        ptVariant *variant = stList_get(variants, i);
        int8_t gt1 = variant->gt[1];
        variant->gt[1] = variant->gt[0];
        variant->gt[0] = gt1;
    }
}


stList *read_phased_variants(char *vcf_path, bool consistent_gt, int min_gq) {

    stList *variants = stList_construct3(0, ptVariant_destruct);

    //open vcf file
    htsFile *fp = hts_open(vcf_path, "rb");

    //read header
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    bcf1_t *rec = bcf_init();


    int64_t ps_pre = -1; // The id of the previous phase block
    int32_t n_ref_ps_variants = 0; // Number of variants with their first genotype as ref in the current phase block
    int32_t n_ps_variants = 0; // Number of total variants in the current phase block
    int32_t n_variants = 0;


    // Each phase block id can take variable number of bytes
    // in the vcf so we are defining ps variables with different
    // sizes for casting properly
    int8_t *ps_int8;
    int16_t *ps_int16;
    int32_t *ps_int32;
    int64_t *ps_int64;
    int64_t ps;

    printf("Reading phased variants\n\n##Phase_Block_ID\tN_Ref\tN_All\tRef_Ratio\tSwapped\n");
    //save each vcf record
    while (bcf_read(fp, hdr, rec) >= 0) {
        //unpack for read REF,ALT,INFO,etc
        bcf_unpack(rec, BCF_UN_STR);
        bcf_unpack(rec, BCF_UN_FMT);

        int32_t *gt_arr = NULL, ngt_arr = 0;
        int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);

        // get genotypes
        bcf_fmt_t *fmt_gt = bcf_get_fmt(hdr, rec, "GT");

        // skip if the number of genotypes is not equal to 2
        if (fmt_gt->n != 2) continue;
        int8_t gt1 = bcf_gt_allele(fmt_gt->p[0]);
        int8_t gt2 = bcf_gt_allele(fmt_gt->p[1]);

        // check if the variant is phased or not
        // an issue was opened here
        // https://github.com/samtools/htslib/issues/1113
        // saying that the 1st genotype never has the phase bit
        // so we look at the 2nd genotype to check if it's phased or not
        int is_phased = bcf_gt_is_phased(fmt_gt->p[1]);

        // skip unphased variants
        if (!is_phased) continue;

        // get contig
        char *contig = bcf_hdr_id2name(hdr, rec->rid);

        // get position
        int32_t pos = rec->pos;

        // get allele frequency
        bcf_fmt_t *fmt_vaf = bcf_get_fmt(hdr, rec, "VAF");
        float *vaf_float_ptr = (float *) fmt_vaf->p;
        float vaf = vaf_float_ptr[0];

        // get genotype quality
        bcf_fmt_t *fmt_gq = bcf_get_fmt(hdr, rec, "GQ");
        assert(fmt_gq->type == BCF_BT_INT8);
        int8_t gq = fmt_gq->p[0];

        // check genotype quality
        if (gq < min_gq) continue;

        // get phase block id
        bcf_fmt_t *fmt_ps = bcf_get_fmt(hdr, rec, "PS");
        switch (fmt_ps->type) {
            case BCF_BT_INT8:
                ps_int8 = (int8_t *) fmt_ps->p;
                ps = ps_int8[0];
                break;
            case BCF_BT_INT16:
                ps_int16 = (int16_t *) fmt_ps->p;
                ps = ps_int16[0];
                break;
            case BCF_BT_INT32:
                ps_int32 = (int32_t *) fmt_ps->p;
                ps = ps_int32[0];
                break;
            case BCF_BT_INT64:
                ps_int64 = (int64_t *) fmt_ps->p;
                ps = ps_int64[0];
                break;
        }

        // get var type
        int8_t var_type = bcf_get_variant_types(rec);

        if (ps_pre != ps && stList_length(variants) > 0) { // if phase block changed
            printf("%d\t%d\t%d\t%.2f", ps_pre, n_ref_ps_variants, n_ps_variants,
                   (double) n_ref_ps_variants / n_ps_variants);
            // if less than half of the variants in the previous
            // phase block (ps_pre) have reference allele as their first
            // genotype then iterate over the variants in that
            // phase block and swap the genotypes for each variant
            if (n_ref_ps_variants < (0.5 * n_ps_variants) && consistent_gt) {
                n_variants = stList_length(variants);
                //Swap genotypes
                ptVariant_swap_gt(variants, n_variants - n_ps_variants, n_variants - 1);
                printf("\tYES\n");
            } else {
                printf("\tNO\n");
            }
            // reset variables related to phase block
            n_ps_variants = 0;
            n_ref_ps_variants = 0;
        }
        // Make a ptVariant structure
        ptVariant *variant = ptVariant_construct(contig, pos, var_type, vaf, gq, ps);

        // add alleles to the variant structure
        // rec->d.allele[0] is REF
        for (int i = 0; i < rec->n_allele; ++i) {
            ptVariant_append_allele(variant, rec->d.allele[i]);
        }

        // add genotypes
        ptVariant_append_gt(variant, gt1);
        ptVariant_append_gt(variant, gt2);

        // add variant to the list
        stList_append(variants, variant);

        // update variables related to phase block
        ps_pre = ps;
        n_ps_variants += 1;
        n_ref_ps_variants += gt1 == 0 ? 1 : 0;
    }

    printf("%d\t%d\t%d\t%.2f", ps_pre, n_ref_ps_variants, n_ps_variants,
           (double) n_ref_ps_variants / n_ps_variants);
    // check the last phase block
    if (n_ref_ps_variants < (0.5 * n_ps_variants) && consistent_gt) {
        n_variants = stList_length(variants);
        //Swap genotypes
        ptVariant_swap_gt(variants, n_variants - n_ps_variants, n_variants - 1);
        printf("\tYES\n");
    } else {
        printf("\tNO\n");
    }

    // Free variant record and vcf header
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    int ret;
    if ((ret = hts_close(fp))) {
        fprintf(stderr, "hts_close(%s): non-zero status %d\n", vcf_path, ret);
        exit(ret);
    }
    return variants;
}

stList *filter_ref_variants(stList *variants) {
    stList *selected_variants = stList_construct3(0, ptVariant_destruct);
    for (int i = 0; i < stList_length(variants); i++) {
        ptVariant *variant = stList_get(variants, i);
        if (variant->gt[0] == 0) continue;
        ptVariant *variant_copy = ptVariant_copy(variant);
        stList_append(selected_variants, variant_copy);
    }
    return selected_variants;
}


stHash *extract_variant_ref_blocks(stList *variants, const faidx_t *fai, int min_margin) {
    stHash *variant_blocks = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL,
                                               (void (*)(void *)) stList_destruct);
    ptVariant *variant = NULL;
    ptVariant *pre_variant = NULL;
    stList *blocks;
    int rfs = -1;
    int rfe = -1;
    int contig_len;
    ptBlock *block = NULL;
    for (int i = 0; i < stList_length(variants); i++) {
        variant = stList_get(variants, i);
        //ptVariant_print(variant);
        // if this is the first variant
        // or the contig has changed
        if (i == 0 || strcmp(variant->contig, pre_variant->contig) != 0) {
            // make a list of blocks for the new contig
            blocks = stList_construct3(0, ptBlock_destruct);
            char *contig = malloc(200);
            strcpy(contig, variant->contig);
            contig_len = faidx_seq_len64(fai, contig);
            // add the blocks list to the hash table with the key equal to contig name
            stHash_insert(variant_blocks, contig, blocks);
        }
        char *allele_ref = stList_get(variant->alleles, 0);
        char *allele_alt = stList_get(variant->alleles, variant->gt[0]);
        // Adjust margin if variant is longer than min_margin
        int margin = max(min_margin, variant->longest_allele_len);
        if (strlen(allele_ref) <= strlen(allele_alt)) { // insertion or snp
            rfs = variant->pos - margin;
            rfe = variant->pos + margin;
        } else if (strlen(allele_alt) < strlen(allele_ref)) { // deletion
            rfs = variant->pos - margin;
            rfe = variant->pos + strlen(allele_ref) - strlen(allele_alt) + margin;
        }
        // correct rfs and rfe if exceeded feasible coordinates
        rfs = rfs < 0 ? 0 : rfs;
        rfe = contig_len < rfe ? contig_len - 1 : rfe;
        // make a ptBlock struct
        block = ptBlock_construct(rfs, rfe,
                                  -1, -1, // sqs and sqe
                                  -1, -1); // rfd_s and rfd_e
        // make a list of variants for the new block
        stList *vars = stList_construct3(0, ptVariant_destruct);
        stList_append(vars, ptVariant_copy(variant));
        // save variant records as the block data
        ptBlock_set_data(block, (void *) vars, ptVariant_destruct_stList, ptVariant_copy_stList,
                         ptVariant_extend_stList);
        stList_append(blocks, block);
        pre_variant = variant;
    }

    // sort variants in each block
    // sort blocks
    // merge overlapping blocks
    stHashIterator *it = stHash_getIterator(variant_blocks);
    char *key;
    while ((key = stHash_getNext(it)) != NULL) {
        // stHash_remove does two things:
        // remove the key from table
        // return the value
        // we will sort and merge blocks (the returned value)
        // then insert the same key but this time with merged blocks as its value
        blocks = stHash_remove(variant_blocks, key);
        for (int i = 0; i < stList_length(blocks); i++) {
            block = stList_get(blocks, i);
            //printf("##%d %d\n", i, stList_length(block->data));
            stList_sort((stList *) block->data, ptVariant_cmp);
        }
        // sort blocks by rfs
        stList_sort(blocks, ptBlock_cmp_rfs);
        // merge blocks by rfs and rfe
        stList *merged_blocks = ptBlock_merge_blocks(blocks, ptBlock_get_rfs, ptBlock_get_rfe, ptBlock_set_rfe);
        // insert the same key
        stHash_insert(variant_blocks, key, merged_blocks);
        // destroy old unmerged blocks
        stList_destruct(blocks);
    }
    return variant_blocks;
}


char *fetch_read_seq(ptAlignment *alignment, ptBlock *block) {
    if (block->sqe < block->sqs) {
        char *empty_seq = malloc(sizeof(char));
        empty_seq[0] = '\0';
        return empty_seq;
    }
    uint8_t *seq = bam_get_seq(alignment->record);
    //printf("block->sqs - block->sqe= %d-%d\t len=%d\n", block->sqs , block->sqe, block->sqe - block->sqs + 1);
    assert(block->sqe < alignment->record->core.l_qseq);
    assert(0 <= block->sqs);
    int block_len = block->sqe - block->sqs + 1;
    //allocate read sequence
    char *block_seq = (uint8_t *) malloc(block_len + 1);
    for (int k = 0; k < block_len; k++)
        block_seq[k] = seq_nt16_str[bam_seqi(seq, block->sqs + k)];
    block_seq[block_len] = '\0';
    return block_seq;
}


char *fetch_corrected_ref_seq(const faidx_t *fai, ptBlock *block, char *contig_name) {
    if (block->rfe < block->rfs) {
        char *empty_seq = malloc(sizeof(char));
        empty_seq[0] = '\0';
        return empty_seq;
    }
    int offset_orig = 0;
    int offset_corrected = 0;
    stList *vars = (stList *) block->data;
    ptVariant *var = NULL;
    char *reg = malloc(200);
    memset(reg, '\0', 200);
    sprintf(reg, "{%s}:%d-%d", contig_name, block->rfs + 1, block->rfe + 1);
    int len;
    //printf("%s\n",reg);
    char *seq_orig = fai_fetch(fai, reg, &len);
    char *seq_corrected = malloc((strlen(seq_orig) * 2) * sizeof(char));
    int size = strlen(seq_orig) * 2;
    memset(seq_corrected, '\0', size);
    // copy the bases before the first variant
    //memcpy(seq_corrected, seq_orig, var->pos - block->rfs);
    //offset_orig += var->pos - block->rfs;
    //offset_corrected += var->pos - block->rfs;
    // apply variations and copy the bases between variants
    for (int i = 0; i < stList_length(vars); i++) {
        var = stList_get(vars, i);
        if (block->rfs + offset_orig < var->pos) {
            if (size <= var->pos - block->rfs - offset_orig + offset_corrected) {
                seq_corrected = realloc(seq_corrected, size * 2);
                size *= 2;
            }
            memcpy(seq_corrected + offset_corrected, seq_orig + offset_orig, var->pos - block->rfs - offset_orig);
            //printf("pos=%d, rfs=%d\n", var->pos , block->rfs);
            offset_corrected += var->pos - block->rfs - offset_orig;
            offset_orig = var->pos - block->rfs;
            seq_corrected[offset_corrected] = '\0';
            //printf("orig=%d, corrected=%d\n", offset_orig, offset_corrected);
            //printf("@@%s\n", seq_corrected);
        }
        char *allele_ref = stList_get(var->alleles, 0);
        char *allele_alt = stList_get(var->alleles, var->gt[0]);
        if (size <= strlen(allele_alt) + strlen(seq_corrected) + offset_corrected) {
            seq_corrected = realloc(seq_corrected, size * 2);
            size *= 2;
        }
        memcpy(seq_corrected + offset_corrected, allele_alt, strlen(allele_alt));
        offset_orig += strlen(allele_ref);
        offset_corrected += strlen(allele_alt);
        //seq_corrected[offset_corrected] = '\0';
        //printf("@#%s\n", seq_corrected);
    }
    if (offset_orig < (block->rfe - block->rfs + 1)) { // last bases to be copied
        strcpy(seq_corrected + offset_corrected, seq_orig + offset_orig);
        offset_corrected += strlen(seq_orig + offset_orig);
    }
    seq_corrected[offset_corrected] = '\0';

    //printf("orig\t%s\n", seq_orig);
    //printf("corr\t%s\n", seq_corrected);
    return seq_corrected;
}


stList *project_ref_blocks_to_read(ptAlignment *alignment, stList *ref_variant_blocks) {
    bam1_t *b = alignment->record;
    // Find the block that starts after the given alignment
    // based on the reference coordinates
    int j = 0;
    ptBlock *ref_block = stList_get(ref_variant_blocks, j);
    while (ref_block->rfs < alignment->rfs) {
        if (stList_length(ref_variant_blocks) - 1 == j) {
            ref_block = NULL;
            break;
        }
        j += 1;
        ref_block = stList_get(ref_variant_blocks, j);
    }
    // If the alignment does not include any block
    if (ref_block && alignment->rfe < ref_block->rfe) {
        ref_block = NULL;
    }
    // Initiate a list for keeping the blocks with read coordinates
    stList *read_blocks = stList_construct3(0, ptBlock_destruct);

    // If the alignment does not include any block return empty list
    if (ref_block == NULL) {
        return read_blocks;
    }

    //printf("first block: %d\t%d\n", ref_block->rfs, ref_block->rfe);
    //printf("alignment interval: %d\t%d\n", alignment->rfs, alignment->rfe);
    ptBlock *read_block;
    int rds_f = -1;
    int rde_f = -1;
    int cigar_rds_f;
    int cigar_rde_f;
    // Iterate over cigar operations
    ptCigarIt *cigar_it = ptCigarIt_construct(b, true, true);
    // Since projection is from reference to read coordinates
    // there is no need to do anything about insertions
    // Insertions do not advance in reference coordinates
    int i = 0;
    while (ptCigarIt_next(cigar_it)) {
        if (bam_is_rev(b)) {
            cigar_rds_f = -1 * cigar_it->rde_f;
            cigar_rde_f = -1 * cigar_it->rds_f;
        } else {
            cigar_rds_f = cigar_it->rds_f;
            cigar_rde_f = cigar_it->rde_f;
        }
        if (cigar_it->op == BAM_CMATCH ||
            cigar_it->op == BAM_CEQUAL ||
            cigar_it->op == BAM_CDIFF ||
            cigar_it->op == BAM_CDEL) {
            /*if(i == 0){
                printf("first cigar op: %d\t%d\n", cigar_it->rfs, cigar_it->rfe);
                i=1;
            }*/
            //The loop below iterates the blocks' overlap with the current cigar operation
            //There may exist multiple blocks within the current op
            while (ref_block && ref_block->rfe <= cigar_it->rfe) {
                //printf("block end cigar op: %d\t%d\n", cigar_it->rfs, cigar_it->rfe);

                // update start locations of the projected coordinates
                if (cigar_it->rfs <= ref_block->rfs) {
                    /*printf("rds_f updated\n");
                    printf("\tcigar op: %d\t%d\n", cigar_it->rfs, cigar_it->rfe);
                    printf("\tblock: %d\t%d\n", ref_block->rfs, ref_block->rfe);*/
                    rds_f = cigar_it->op == BAM_CDEL ? cigar_rds_f : cigar_rds_f + (ref_block->rfs - cigar_it->rfs);
                    //printf("rds_f = %d\n",rds_f);
                    //sqs = cigar_it->op == BAM_CDEL ? cigar_it->sqs : cigar_it->sqs + (ref_block->rfs - cigar_it->rfs);
                }
                /*printf("rde_f updated\n");
                                printf("\tcigar op: %d\t%d\n", cigar_it->rfs, cigar_it->rfe);
                                printf("\tblock: %d\t%d\n", ref_block->rfs, ref_block->rfe);*/
                // update end locations of the projected coordinates
                rde_f = cigar_it->op == BAM_CDEL ? cigar_rde_f : cigar_rds_f + (ref_block->rfe - cigar_it->rfs);
                //printf("rde_f = %d\n",rde_f);
                //sqe = cigar_it->op == BAM_CDEL ? cigar_it->sqe : cigar_it->sqs + (ref_block->rfe - cigar_it->rfs);

                // Construct block with read coordinates
                if (bam_is_rev(b)) { //positive strand
                    read_block = ptBlock_construct(ref_block->rfs,
                                                   ref_block->rfe,
                                                   -1,
                                                   -1,
                                                   -1 * rde_f,
                                                   -1 * rds_f);
                } else { // negative strand
                    read_block = ptBlock_construct(ref_block->rfs,
                                                   ref_block->rfe,
                                                   -1,
                                                   -1,
                                                   rds_f,
                                                   rde_f);

                }

                // append new read block
                if (read_block->rds_f <= read_block->rde_f) {
                    // make a copy of all the variants in the ref block
                    stList *variants_copy = (stList *) ptVariant_copy_stList((stList *) ref_block->data);
                    // add all copied variants to the read block
                    ptBlock_set_data(read_block, (void *) variants_copy, ptVariant_destruct_stList,
                                     ptVariant_copy_stList,
                                     ptVariant_extend_stList);
                    //printf("read block added: %d\t%d\n", read_block->rds_f, read_block->rde_f);
                    // add read block
                    stList_append(read_blocks, read_block);
                } else { // delete block if the whole block is within a deletion
                    ptBlock_destruct(read_block);
                }

                // update block index; j
                if (j < stList_length(ref_variant_blocks) - 1) {
                    j++;
                    ref_block = stList_get(ref_variant_blocks, j);
                } else if (j == stList_length(ref_variant_blocks) - 1) {
                    ref_block = NULL;
                }
                if (ref_block == NULL) break;// no more block remaining so break the cigar iteration
            }//end while
            // if the start of the next block is within the current operation
            if (ref_block &&
                cigar_it->rfs <= ref_block->rfs &&
                ref_block->rfs <= cigar_it->rfe) {
                /*printf("rds_f updated\n");
                printf("\tcigar op: %d\t%d\n", cigar_it->rfs, cigar_it->rfe);
                printf("\tblock: %d\t%d\n", ref_block->rfs, ref_block->rfe);*/
                rds_f = cigar_it->op == BAM_CDEL ? cigar_rds_f : cigar_rds_f + (ref_block->rfs - cigar_it->rfs);
                //printf("rds_f = %d\n",rds_f);
                //sqs = cigar_it->op == BAM_CDEL ? cigar_it->sqs : cigar_it->sqs + (ref_block->rfs - cigar_it->rfs);
            }//end if
        }// end if
    }//end of iteration over cigar ops
    //printf("CIGAR DONE!\n\n\n");
    ptCigarIt_destruct(cigar_it);
    return read_blocks;
}

stList *merge_variant_read_blocks(stList *blocks) {
    // merge blocks by rds_f and rde_f
    return ptBlock_merge_blocks(blocks, ptBlock_get_rds_f, ptBlock_get_rde_f, ptBlock_set_rde_f);;
}

int get_edit_distance(ptAlignment *alignment, faidx_t *fai, ptBlock *block) {
    char *read_seq = fetch_read_seq(alignment, block);
    char *corrected_ref_seq = fetch_corrected_ref_seq(fai, block, alignment->contig);
    //printf("read\t%s\nref\t%s\n", read_seq, corrected_ref_seq);
    int edit_distance = 0;
    if ((0 < strlen(read_seq)) && (0 < strlen(corrected_ref_seq))) {
        EdlibAlignResult result = edlibAlign(read_seq, strlen(read_seq), corrected_ref_seq, strlen(corrected_ref_seq),
                                             edlibDefaultAlignConfig());
        if (!(result.status == EDLIB_STATUS_OK)) {
            fprintf(stderr, "Edlib didn't work!\n");
            exit(EXIT_FAILURE);
        }
        edit_distance = result.editDistance;
        //printf("edit=%d\n\n", edit_distance);
        edlibFreeAlignResult(result);
    } else {
        // at least one of read_seq and corrected_ref_seq is '\0'
        edit_distance = max(strlen(read_seq), strlen(corrected_ref_seq));
    }
    free(corrected_ref_seq);
    free(read_seq);
    return edit_distance;
}


int get_total_edit_distance(ptAlignment *alignment, const faidx_t *fai, char *contig_name,
                            stList *variant_read_blocks, stList *projected_blocks) {
    bam1_t *b = alignment->record;
    int j = bam_is_rev(b) ? stList_length(variant_read_blocks) - 1 : 0;
    int edit_distance = 0;
    ptCigarIt *cigar_it = ptCigarIt_construct(b, true, true);
    int block_rds_f, block_rde_f, cigar_rds_f, cigar_rde_f, rfs, rfe, sqs, sqe;
    ptBlock *read_block = stList_get(variant_read_blocks, j);
    // reverse each block interval to make it work with the projecting algorithm below
    if (bam_is_rev(b)) {
        block_rds_f = -1 * read_block->rde_f;
        block_rde_f = -1 * read_block->rds_f;
    } else {
        block_rds_f = read_block->rds_f;
        block_rde_f = read_block->rde_f;
    }
    //iterate over cigar operations from left to right (w.r.t reference)
    while (ptCigarIt_next(cigar_it)) {
        // reverse each cigar operation interval to make it work with the projecting algorithm below
        if (bam_is_rev(b)) {
            cigar_rds_f = -1 * cigar_it->rde_f;
            cigar_rde_f = -1 * cigar_it->rds_f;
        } else {
            cigar_rds_f = cigar_it->rds_f;
            cigar_rde_f = cigar_it->rde_f;
        }
        // match/ mismatch/ insertion
        // for M and I the main algorithm is the same except
        // some minor parts that are corrected by
        // the conditional statement "cigar_it->op == BAM_CINS ? : "
        if (cigar_it->op == BAM_CMATCH ||
            cigar_it->op == BAM_CEQUAL ||
            cigar_it->op == BAM_CDIFF ||
            cigar_it->op == BAM_CINS) {
            //The loop below iterate the blocks overlap with the current cigar operation
            //There may exist multiple blocks within the current op
            while (read_block && block_rde_f <= cigar_rde_f) {
                // update start locations of the projected coordinates
                if (cigar_rds_f <= block_rds_f) {
                    rfs = cigar_it->op == BAM_CINS ? cigar_it->rfs : cigar_it->rfs + (block_rds_f - cigar_rds_f);
                    sqs = cigar_it->sqs + (block_rds_f - cigar_rds_f);
                }
                // update end locations of the projected coordinates
                rfe = cigar_it->op == BAM_CINS ? cigar_it->rfe : cigar_it->rfs + (block_rde_f - cigar_rds_f);
                sqe = cigar_it->sqs + (block_rde_f - cigar_rds_f);
                /*if(1000 < (rfe-rfs)){
                    printf("####%s\t%d\t%d\n",contig_name,rfs,rfe);
                    printf("%d\t%d\n",read_block->rds_f, read_block->rde_f);
                }*/
                // get block
                ptBlock *block = ptBlock_construct(rfs,
                                                   rfe,
                                                   sqs,
                                                   sqe,
                                                   read_block->rds_f,
                                                   read_block->rde_f);
                stList_append(projected_blocks, block);

                //printf("%d\t%d\n", rfs, rfe);
                // all_variants contains all variants from all regions in the genome
                // that have overlap with alignment of this read block
                stList *all_variants = (stList *) read_block->data;
                //printf("All variants %d\n", stList_length(all_variants));
                //for(int k=0; k < stList_length(all_variants);k++){
                //        ptVariant_print(stList_get(all_variants,k));
                // }
                // kept_variants will contain only the variants from this contig
                // and located within the block
                stList *kept_variants = stList_construct3(0, ptVariant_destruct);
                for (int i = 0; i < stList_length(all_variants); i++) {
                    ptVariant *variant = stList_get(all_variants, i);
                    if (strcmp(variant->contig, contig_name) == 0 &&
                        block->rfs <= variant->pos &&
                        variant->pos <= block->rfe) {
                        stList_append(kept_variants, ptVariant_copy(variant));
                    }
                }
                ptBlock_set_data(block, (void *) kept_variants, ptVariant_destruct_stList, ptVariant_copy_stList,
                                 ptVariant_extend_stList);
                //printf("Kept variants %d\n", stList_length(kept_variants));
                //for(int k=0; k < stList_length(kept_variants);k++){
                //	ptVariant_print(stList_get(kept_variants,k));
                //}
                edit_distance += get_edit_distance(alignment, fai, block);
                // update block index; j
                if (bam_is_rev(b) && j > 0) {
                    j--;
                    read_block = stList_get(variant_read_blocks, j);
                    block_rds_f = -1 * read_block->rde_f;
                    block_rde_f = -1 * read_block->rds_f;
                } else if (!bam_is_rev(b) && j < stList_length(variant_read_blocks) - 1) {
                    j++;
                    read_block = stList_get(variant_read_blocks, j);
                    block_rds_f = read_block->rds_f;
                    block_rde_f = read_block->rde_f;
                } else if (j == 0 || j == stList_length(variant_read_blocks) - 1) {
                    read_block = NULL;
                }
            }
            if (read_block == NULL) break;// no more block remaining so break the cigar iteration
            // if the start of the next block is within the current operation
            if (cigar_rds_f <= block_rds_f &&
                block_rds_f <= cigar_rde_f) {
                rfs = cigar_it->op == BAM_CINS ? cigar_it->rfs : cigar_it->rfs + (block_rds_f - cigar_rds_f);
                sqs = cigar_it->sqs + (block_rds_f - cigar_rds_f);
            }
        }
    }
    ptCigarIt_destruct(cigar_it);
    return edit_distance;
}


stHash *
ptVariant_parse_variants_and_extract_blocks(char *vcf_path, char *bed_path, faidx_t *fai, int min_margin, int min_gq) {
    stList *phased_variants = read_phased_variants(vcf_path, true, min_gq);
    fprintf(stdout, "[%s] Number of parsed phased variants = %d\n", get_timestamp(), stList_length(phased_variants));

    stHash *blocks_per_contig;
    stHash *merged_blocks_per_contig;
    stList *intersected_variants;
    if (bed_path != NULL && bed_path[0] != NULL) {
        blocks_per_contig = ptBlock_parse_bed(bed_path);
        fprintf(stdout, "[%s] Total length of parsed bed tracks = %d\n", get_timestamp(),
                ptBlock_get_total_length_by_rf(blocks_per_contig));

        merged_blocks_per_contig = ptBlock_merge_blocks_per_contig_by_rf(blocks_per_contig);
        fprintf(stdout, "[%s] Total length of merged bed tracks = %d\n", get_timestamp(),
                ptBlock_get_total_length_by_rf(merged_blocks_per_contig));

        intersected_variants = ptVariant_subset_stList(phased_variants, merged_blocks_per_contig);
        fprintf(stdout, "[%s] Number of intersected variants = %d\n", get_timestamp(),
                stList_length(intersected_variants));
    } else {
        intersected_variants = phased_variants;
        fprintf(stdout, "[%s] Number of intersected variants (No BED file was given) = %d\n", get_timestamp(),
                stList_length(intersected_variants));
    }

    stList *selected_variants = filter_ref_variants(intersected_variants);
    fprintf(stdout, "[%s] Number of selected (No REF) variants = %d\n", get_timestamp(),
            stList_length(selected_variants));

    stHash *variant_ref_blocks = extract_variant_ref_blocks(selected_variants, fai, min_margin);
    fprintf(stdout, "[%s] Variant blocks are created on the reference coordinates (min_margin = %d).\n",
            get_timestamp(), min_margin);
    fprintf(stdout, "[%s] Total length of variant blocks: %d\n",
            get_timestamp(), ptBlock_get_total_length_by_rf(variant_ref_blocks));

    // We can free these lists because variant_ref_blocks has all the necessary variants
    stList_destruct(phased_variants);
    stList_destruct(selected_variants);
    if (bed_path != NULL) {
        stList_destruct(intersected_variants);
        stHash_destruct(blocks_per_contig);
        stHash_destruct(merged_blocks_per_contig);
    }

    return variant_ref_blocks;

}


void ptVariant_save_variant_ref_blocks(stHash *variant_ref_blocks, char *bed_path) {

    char *contig;
    stList *contig_blocks;
    ptBlock *block;
    stList *vars;
    FILE *fp = fopen(bed_path, "w");
    if (fp == NULL) {
        fprintf(stderr, "[%s] Failed to open file %s.\n", get_timestamp(), bed_path);
    }

    // iterate over contigs with variant blocks and print them in the bed file
    stHashIterator *it = stHash_getIterator(variant_ref_blocks);
    while ((contig = stHash_getNext(it)) != NULL) {
        contig_blocks = stHash_search(variant_ref_blocks, contig);
        for (int i = 0; i < stList_length(contig_blocks); i++) {
            block = stList_get(contig_blocks, i);
            vars = (stList *) block->data;
            fprintf(fp, "%s\t%d\t%d\tN_VARS=%d\n", contig, block->rfs, block->rfe + 1, stList_length(vars));
        }
    }
    fclose(fp);
}


stList *ptVariant_get_merged_variant_read_blocks(stHash *variant_ref_blocks_per_contig, ptAlignment **alignments,
                                                 int alignments_len) {

    // Make a list of lists for saving read blocks; one list per alignment
    // Each internal list contains the blocks in the read coordinates for one alignment
    stList **read_blocks_per_alignment = (stList *) malloc(alignments_len * sizeof(stList * ));
    for (int i = 0; i < alignments_len; i++) {
        stList *blocks_contig = stHash_search(variant_ref_blocks_per_contig, alignments[i]->contig);
        if (blocks_contig == NULL) {
            read_blocks_per_alignment[i] = stList_construct3(0, ptBlock_destruct);
        } else {
            read_blocks_per_alignment[i] = project_ref_blocks_to_read(alignments[i], blocks_contig);
        }
    }
    // Find the maximum start point and minimum end point for alignments
    // Alignments of the same read may not cover the same regions of the read
    int max_rds_f = alignments[0]->rds_f;
    int min_rde_f = alignments[0]->rde_f;
    for (int i = 0; i < alignments_len; i++) {
        max_rds_f = max_rds_f < alignments[i]->rds_f ? alignments[i]->rds_f : max_rds_f;
        min_rde_f = alignments[i]->rde_f < min_rde_f ? alignments[i]->rde_f : min_rde_f;
    }
    // Put all read blocks in a single list
    stList *all_read_blocks = stList_construct3(0, ptBlock_destruct);
    for (int i = 0; i < alignments_len; i++) {
        int len = stList_length(read_blocks_per_alignment[i]);
        for (int j = 0; j < len; j++) {
            ptBlock *b = stList_get(read_blocks_per_alignment[i], j);
            // add block only if it is within the region which has alignments
            // to all haplotypes
            if ((max_rds_f <= b->rds_f) && (b->rde_f <= min_rde_f)) {
                stList_append(all_read_blocks, ptBlock_copy(b));
            }
        }
    }
    // Sort all read blocks by start pos
    stList_sort(all_read_blocks, ptBlock_cmp_rds_f);
    // Merge read blocks
    stList *read_blocks_merged = merge_variant_read_blocks(all_read_blocks);

    // free read_blocks_per_alignment
    for (int i = 0; i < alignments_len; i++) {
        stList_destruct(read_blocks_per_alignment[i]);
    }
    free(read_blocks_per_alignment);

    //free all_read_blocks
    stList_destruct(all_read_blocks);


    return read_blocks_merged;
}


stList **
set_scores_as_edit_distances(stList *read_blocks_merged, ptAlignment **alignments, int alignments_len, faidx_t *fai) {
    stList **projected_blocks_per_alignment_idx = (stList **) malloc(alignments_len * sizeof(stList * ));
    // get edit distances of the variant blocks for each alignment
    for (int i = 0; i < alignments_len; i++) {
        projected_blocks_per_alignment_idx[i] = stList_construct3(0, ptBlock_destruct);
        alignments[i]->score =
                -1 *
                get_total_edit_distance(alignments[i], fai, alignments[i]->contig, read_blocks_merged,
                                        projected_blocks_per_alignment_idx[i]);
    }
    return projected_blocks_per_alignment_idx;
}


bool overlap_variant_ref_blocks(stHash *variant_ref_blocks_per_contig, ptAlignment **alignments, int alignments_len) {
    // Iterate over alignments
    for (int i = 0; i < alignments_len; i++) {
        stList *blocks_contig = stHash_search(variant_ref_blocks_per_contig, alignments[i]->contig);
        // If there is no block then continue
        if (blocks_contig == NULL || stList_length(blocks_contig) == 0) continue;
        ptBlock *first_block = stList_get(blocks_contig, 0);
        ptBlock *last_block = stList_get(blocks_contig, stList_length(blocks_contig) - 1);
        // If the alignment ends after the first block or starts after the last block then continue
        if ((last_block->rfe < alignments[i]->rfs) || (alignments[i]->rfe < first_block->rfs)) continue;
        // Check all blocks
        for (int j = 0; j < stList_length(blocks_contig); j++) {
            ptBlock *block = stList_get(blocks_contig, j);
            // If a block was completely within an alignment then return true
            if ((alignments[i]->rfs <= block->rfs) && (block->rfe <= alignments[i]->rfe)) return true;
        }
    }
    return false;
}


void ptVariant_extend_stList(void *curr_vars_, void *new_vars_) {
    stList *curr_vars = (stList *) curr_vars_;
    stList *new_vars = (stList *) new_vars_;
    for (int i = 0; i < stList_length(new_vars); i++) {
        ptVariant *new_var = stList_get(new_vars, i);
        // if the new variant does not exist make a copy and add
        // to the first variants list
        if (!ptVariant_exist_in_list(new_var, curr_vars)) {
            stList_append(curr_vars, ptVariant_copy(new_var));
        }
    }
}


void *ptVariant_copy_stList(void *vars_) {
    stList *vars = (stList *) vars_;
    stList *vars_copy = stList_construct3(0, ptVariant_destruct);
    for (int i = 0; i < stList_length(vars); i++) {
        stList_append(vars_copy, ptVariant_copy(stList_get(vars, i)));
    }
    return (void *) vars_copy;
}

void *ptVariant_destruct_stList(void *vars_) {
    stList *vars = (stList *) vars_;
    stList_destruct(vars);
}


// variants should be sorted by contig and pos
// blocks should be sorted by rfs
stList *ptVariant_subset_stList(stList *variants, stHash *blocks_per_contig) {
    char contig[200];
    contig[0] = '\0';
    int j = 0;
    ptBlock *block;
    stList *blocks;
    ptVariant *variant;
    stList *subset_variants = stList_construct3(0, ptVariant_destruct);
    for (int i = 0; i < stList_length(variants); i++) {
        variant = stList_get(variants, i);
        if (strcmp(contig, variant->contig) != 0) {
            strcpy(contig, variant->contig);
            blocks = stHash_search(blocks_per_contig, contig);
            // handle the case when the whole contig is not in the bed file
            if (blocks != NULL) {
                j = 0;
                block = stList_get(blocks, j);
            }
        }
        if (blocks == NULL) continue;
        while (block->rfe < variant->pos && j < stList_length(blocks) - 1) {
            j += 1;
            block = stList_get(blocks, j);
        }
        if (block->rfe < variant->pos) continue;
        if (block->rfs <= variant->pos && variant->pos <= block->rfe) {
            stList_append(subset_variants, ptVariant_copy(variant));
        }
    }
    return subset_variants;
}

stHash *ptVariant_copy_stHash_blocks_per_contig(stHash *blocks_per_contig) {
    stHash *copy_blocks_per_contig = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL,
                                                       (void (*)(void *)) stList_destruct);
    stList *blocks;
    char *contig_name;
    stHashIterator *it = stHash_getIterator(blocks_per_contig);
    while ((contig_name = stHash_getNext(it)) != NULL) {
        blocks = stHash_search(blocks_per_contig, contig_name);
        stHash_insert(copy_blocks_per_contig, contig_name, ptBlock_copy_stList(blocks));
    }
    return copy_blocks_per_contig;
}