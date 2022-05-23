## Flagger-based evaluation pipeline

### Overview
Here is a description of a read-based pipeline that can detect different types of mis-assemblies in a draft diploid assembly. One core component of this pipeline is another pipeline named [**Flagger**](https://github.com/mobinasri/flagger/tree/main/docs/coverage). Flagger recieves the read alignments to the draft assembly and partition the assembly into 4 main components; erroneous, (falsely) duplicated, haploid and collapsed.

This evaluation has 5 steps:
- Align long reads to the diploid assembly
- Phase the ambiguous alignments using [secphase](https://github.com/mobinasri/secphase) (Optional)
- Remove the alignments with alternative alleles
- Run Flagger on the assembly using the alignments with no alternative alleles
- Combine the Flagger outputs

### 1. Align Reads
The ONT and HiFi reads can be aligned to a diploid assembly (~ 6Gbases long) with a long read aligner like winnowmap. Since the assembly is diploid the expected base-level coverage should be half of the sequencing coverage.
Here are the commands for producing the alignments (taken from the [winnowmap docs](https://github.com/marbl/Winnowmap)):
```` 
  # making the k-mer table with meryl
  meryl count k=15 output merylDB asm.fa
  meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt
  
  # alignment with winnowmap
  winnowmap -W repetitive_k15.txt -ax [map-ont | map-pb] -Y -L --eqx --cs -I8g <(cat pat_asm.fa mat_asm.fa) reads.fq.gz | \
    samtools view -hb > read_alignment.bam
````
### 2. Relocalize wrongly phased reads
In this step we use Secphase to rephase the reads with multiple alignments. In other words all the secondary and primary
alignments of the same read are scored based on marker consistency and 
the alignment with the highest score is selected as the primary alignment. The output of this section is 
a corrected version of the input bam file, in which the primary and secondary alignments are swapped 
whenever neccessary.

More information about Secphase is available [here](https://github.com/mobinasri/secphase)

### 3. Call variants 
By calling variants it is possible to detect the regions that needs polishing or the regions that have alignments from the wrong haplotype. It is recommeneded to use [Deepvariant](https://github.com/google/deepvariant) for calling variants from HiFi alignments and [Pepper-Margin-Deepvariant](https://github.com/kishwarshafin/pepper) for ONT. 
````
## For HiFi
## Taken from deepvariant doc
BIN_VERSION="1.3.0"
docker run \
  -v ${INPUT_DIR}:/input \
  -v ${OUTPUT_DIR}:/output \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type="PACBIO" \
  --ref="/input/${ASSEMBLY_FASTA}" \
  --reads="/input/${INPUT_BAM}" \
  --output_vcf="/output/${OUTPUT_VCF}" \
  --make_examples_extra_args="keep_supplementary_alignments=true, min_mapping_quality=0" \
  --call_variants_extra_args="use_openvino=true" \
  --num_shards=$(nproc) \
  --dry_run=false 
  
## For ONT
## Taken from pepper-margin-deepvariant doc
sudo docker run \
  -v ${INPUT_DIR}:/input \
  -v ${OUTPUT_DIR}:/output \
  kishwars/pepper_deepvariant:r0.6 \
  run_pepper_margin_deepvariant call_variant \
  -b "/input/${INPUT_BAM}" \
  -f "/input/${ASSEMBLY_FASTA}" \
  -o "/output" \
  -t $(nproc) \
  --ont_r9_guppy5_sup \
  --pepper_include_supplementary \
  --dv_min_mapping_quality 0 \
  --pepper_min_mapping_quality 0 \
  
 
# --ont_r9_guppy5_sup is preset for ONT R9.4.1 Guppy 5 "Sup" basecaller
# for ONT R10.4 Q20 reads: --ont_r10_q20
````

Note that for both variant callers, the minimum mapping quality is set to 0 which is neccessary to do if the assembly under evaluation is diploid.

### 4. Remove the alignments with alternative alleles
The called variants are then filtered to include only the biallelic snps.
````
## Get the biallelic snps
bcftools view -Ov -f PASS -m2 -M2 -v snps ${OUTPUT_VCF} > ${SNPS_VCF}
````
By having the biallelic snps it is possible to find the alignments with alternative alleles, remove them from the bam file and produce a new bam file.
`filter_alt_reads` is a program that can be used for this aim.
```
## Run filter_alt_reads to get a bam file with no alternative-contained alignments
docker run \
 -v ${INPUT_DIR}:/input \
 -v ${OUTPUT_DIR}:/output \
 quay.io/masri2019/hpp_coverage:latest \
 filter_alt_reads \
 -i "/input/${INPUT_BAM}" \
 -o "/output/${ALT_FILTERED_BAM}"
 -f "/output/${ALT_BAM}"
 -v "${SNPS_VCF}"
 -t $(nproc)
```
`${ALT_FILTERED_BAM}` is the bam file with no alignments that contain alternative alleles and `${ALT_BAM}` includes the removed alignments.

### 5. Run Flagger on the alignments with no alternative allele
`${ALT_FILTERED_BAM}` is then used as input to Flagger. Flagger outputs a bed file for each of the 4 components; 
erroneous, duplicated, haploid and collapsed. Any component other than the haploid one is pointing to unreliable blocks in
assembly. The 4 components are explained in detail [here](https://github.com/mobinasri/flagger/tree/main/docs/coverage#2-coverage-distribution-and-fitting-the-mixture-model). 

More information about Flagger is available [here]()

### Components

|Component| Status| Color |Description|
|:--------|:-----|:-----|:----------|
|Err  |**Erroneous** |Red| This block has low read coverage. If it is located in the middle of a contig it could be either a misjoin or a region that needs polishing|
|Dup  |**Duplicated** |Orange| This block is potentially a false duplication of another block. It should mainly include low-MAPQ alignments with half of the expected coverage. Probably one of the copies has to be polished to fix this issue|
|Hap  | **Haploid** |Green| This block is correctly assembled and has the expected read coverage |
|Col |**Collapsed** |Purple| Two or more highly similar haplotypes are collapsed into this block |

Each of these components has their own color when they are shown in the IGV or the UCSC Genome Browser.


### Running Flagger on HPRC assemblies

The haplotype-resolved assemblies of the HPRC-Y1 samples and their corresponding data sets are available in

https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=working/HPRC/

For more details read this github page:

https://github.com/human-pangenomics/HPP_Year1_Assemblies

We have used the Genbank version of the HPRC-Y1 assemblies.

The results are available in 
https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/e9ad8022-1b30-11ec-ab04-0a13c5208311--COVERAGE_ANALYSIS_Y1_GENBANK/FLAGGER/

