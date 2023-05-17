## Evaluating dual/diploid assemblies with Flagger

### Overview
Here is a description of a read-based pipeline that can detect different types of mis-assemblies in a draft dual/diploid assembly. (*What is a dual assembly? Read [this page](https://lh3.github.io/2021/10/10/introducing-dual-assembly)*). One core component of this pipeline is another pipeline named [**Flagger**](https://github.com/mobinasri/flagger/tree/main/docs/flagger). Flagger recieves the read alignments to a draft diploid assembly, detects the anomalies in the read coverage along the assembly and partition the assembly into 4 main components; erroneous, (falsely) duplicated, haploid and collapsed.



This evaluation has 5 steps:
- Align long reads to the diploid assembly
- Phase and relocalize the reads with secondary alignments using [secphase](https://github.com/mobinasri/secphase) (Optional)
- Call and filter variants (Optional)
- Use biallelic SNVs to remove the alignments with alternative alleles (Optional)
- Run Flagger in two modes
  - Using all alignments
  - Using the alignments with no alternative alleles (optinal but 3rd and 4th steps are required, recommended only when base-level accuracy is concerned)

### 1. Align long reads
The ONT and HiFi reads can be aligned to a dual/diploid assembly (~ 6Gbases long in human) with a long read aligner like winnowmap and minimap2. Since the assembly is dual/diploid the expected base-level coverage should be half of the sequencing coverage.
Here are the commands for producing the alignments (taken from the [winnowmap docs](https://github.com/marbl/Winnowmap)):
```` 
  # making the k-mer table with meryl
  meryl count k=15 output merylDB asm.fa
  meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt
  
  # alignment with winnowmap
  winnowmap -W repetitive_k15.txt -ax [map-ont | map-pb] -Y -L --eqx --cs -I8g <(cat pat_asm.fa mat_asm.fa) reads.fq.gz | \
    samtools view -hb > read_alignment.bam
````
Any other appropriate long read aligner can also be employed in this step.

### 2. Relocalize wrongly phased reads (Optional)
In this step we use Secphase to phase and relocalize the reads with multiple alignments. To be more precise all the secondary and primary
alignments of the same read are scored based on marker consistency and 
the alignment with the highest score is selected as the primary alignment. The output of this section is 
a corrected version of the input bam file, in which the primary and secondary alignments are swapped 
whenever neccessary. Secphase can work only if the secondary alignments are available with the full sequence and base quality array.

More information about Secphase is available [here](https://github.com/mobinasri/secphase)

### 3. Call and filter variants (Optional)
By calling variants it is possible to detect the regions that need polishing or the regions with alignments from the wrong haplotype. It is recommeneded to use [Deepvariant](https://github.com/google/deepvariant) for calling variants with HiFi alignments and [Pepper-Margin-Deepvariant](https://github.com/kishwarshafin/pepper) for ONT. 
````
## For HiFi
## Taken from deepvariant doc
BIN_VERSION="1.4.0"
docker run \
  -v ${INPUT_DIR}:/input \
  -v ${OUTPUT_DIR}:/output \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type="PACBIO" \
  --ref="/input/${ASSEMBLY_FASTA}" \
  --reads="/input/${INPUT_BAM}" \
  --output_vcf="/output/${OUTPUT_VCF}" \
  --make_examples_extra_args="keep_supplementary_alignments=true,min_mapping_quality=0" \
  --call_variants_extra_args="use_openvino=false" \
  --num_shards=$(nproc) \
  --dry_run=false 
  
## For ONT
## Taken from pepper-margin-deepvariant doc
sudo docker run \
  -v ${INPUT_DIR}:/input \
  -v ${OUTPUT_DIR}:/output \
  kishwars/pepper_deepvariant:r0.8 \
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

Note that for both variant callers, the minimum mapping quality is set to 0 which is neccessary to do if the assembly under evaluation is dual/diploid.

The called variants are then filtered to include only the biallelic snps with high quality and frequency.
````
## Get the biallelic snps
bcftools view -Ov -f PASS -m2 -M2 -v snps -e 'FORMAT/VAF<~{vafCutoff} | FORMAT/GQ<~{qCutoff}' ${OUTPUT_VCF} > ${SNPS_VCF}
````

### 4. Remove the alignments with alternative alleles (Optional)
By having the biallelic snps it is possible to find the alignments with alternative alleles, remove them from the bam file.
`filter_alt_reads` is a program that can be used for this aim.
```
## Run filter_alt_reads to get a bam file with no alternative-contained alignments
docker run \
 -v ${INPUT_DIR}:/input \
 -v ${OUTPUT_DIR}:/output \
 mobinasri/flagger:v0.3.0 \
 filter_alt_reads \
 -i "/input/${INPUT_BAM}" \
 -o "/output/${ALT_FILTERED_BAM}"
 -f "/output/${ALT_BAM}"
 -v "${SNPS_VCF}"
 -t $(nproc)
 -m 1000 
 -r 0.4
```

For each alignment `filter_alt_reads` iterates over the CIGAR string and clusters the snps closer than the number given to the `-m` parameter. That alignment will be removed if it encompasses a snp cluster in which more than `-r` ratio of the snps have alternative alleles. 
`${ALT_FILTERED_BAM}` is the cleaned bam file and `${ALT_BAM}` includes the removed alignments.

### 5. Run Flagger
#### 5.1 (Mode 1) Using all alignments
The produced alignment file (`${INPUT_BAM}` prior to step 4) can be used as the input to Flagger. Flagger outputs a bed file with 5 labels; 
erroneous (Err), duplicated (Dup), haploid (Hap), collapsed (Col) and unkown (Unk). Any component other than the haploid one is pointing to unreliable blocks in assembly and unkown label is for the bases couldn't be assigned confidently. The 4 components are explained in detail [here](https://github.com/mobinasri/flagger/tree/main/docs/coverage#2-coverage-distribution-and-fitting-the-mixture-model). 

#### 5.2 (Mode 2) Using alignments with no alternative alleles
`${ALT_FILTERED_BAM}` (after step 4) can also be used as the input to Flagger. Some of the regions flagged as collapsed in `${INPUT_BAM}` may turn out to be haploid in `${ALT_FILTERED_BAM}`, which means in the original alignment they had mismapped reads from other regions. There may be some regions flagged as haploid originally but turn out to be erroneous in `${ALT_FILTERED_BAM}`, which indicates regions that needs polishing.

More information about Flagger is available [here](https://github.com/mobinasri/flagger/tree/main/docs/flagger)

#### Running pipeline with WDL

It is easier to run the pipeline using the WDLs described below. A WDL file can be run locally using Cromwell, which is an open-source Workflow Management System for bioinformatics. The latest releases of Cromwell are available [here](https://github.com/broadinstitute/cromwell/releases) and the documentation is available [here](https://cromwell.readthedocs.io/en/stable/CommandLine/).

It is recommended to run the whole pipeline using [flagger_end_to_end.wdl](https://github.com/mobinasri/flagger/blob/main/wdls/workflows/flagger_end_to_end.wdl). It will run Flagger in both modes mentioned above. `flagger_end_to_end.wdl` requires the alignment of each haplotype assembly to a reference like [chm13v2.0](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz). It will use these alignments for extracting the potentially coverage-biased regions and also stratifying the final results based on CenSat, SD and sex annotations.

Recommended values for the parameters of flagger_end_to_end.wdl:


|Parameter| Value|
|:--------|:-----|
|FlaggerEndToEnd.maxReadDivergence | 0.02 for HiFi and 0.09 for ONT|
|FlaggerEndToEnd.variantCaller | "dv" for HiFi and "pmdv" for ONT|
|preprocess.deepVariantModelType | Should be set based on the latest version of deepvariant ("PACBIO" for v1.4.0)|
|preprocess.pepperModelType | Should be set based on the ONT guppy version  (read https://github.com/kishwarshafin/pepper) "--ont_r9_guppy5_sup" for R9-guppy5 |
|preprocess.moreOptions | "-m 1000 -r 0.4" |
|preprocess.qCutoff | 10 |
|preprocess.vafCutoff | 0.3|
|FlaggerEndToEnd.refBiasedBlocksBedArray | [ "gs://masri/flagger/v0.3.0/chm13v1.1_hifi_r1_high_biased.bed", "gs://masri/flagger/v0.3.0/chm13v1.1_hifi_r2_low_biased.bed" ] for HiFi and  [ "gs://masri/flagger/v0.3.0/chm13v1.1_ont_r2_low_biased.bed"] for ONT |
|FlaggerEndToEnd.refBiasedRegionFactorArray | [ 1.25, 0.75 ] for HiFi and [0.75] for ONT |
|FlaggerEndToEnd.refBiasedRegionNameArray | [ "hifi_biased_high", "hifi_biased_low" ] for HiFi and  ["ont_biased_low" ] for ONT |
|FlaggerEndToEnd.refCntrBed | "gs://masri/flagger/v0.3.0/chm13v2.0.censat.bed" |
|FlaggerEndToEnd.refCntrCtBed | "gs://masri/flagger/v0.3.0/chm13v2.0.ct.bed" |
|FlaggerEndToEnd.refSDBed| "gs://masri/flagger/v0.3.0/chm13v2.0.sd.bed" |
|FlaggerEndToEnd.refSexBed| "gs://masri/flagger/v0.3.0/chm13v2.0.sex.bed" |
|FlaggerEndToEnd.refName | "chm13v2.0"|
|FlaggerEndToEnd.secphaseOptions | "--hifi" for HiFi and "--ont" for ONT |


All files with gs urls are publicly accessible so if you are running the WDL on Terra you can use the same urls. They are also available in the directories `misc/annotations` and `misc/biased_regions` of this repository for those who want to run locally. This WDL also needs the alignment of each haplotype to the reference (like chm13v2.0). Those alignments can be produced using [asm2asm_aligner.wdl](https://github.com/mobinasri/flagger/blob/main/wdls/tasks/alignment/asm2asm_aligner.wdl). Here is the list of recommended parametes for this workflow:
|Parameter| Value|
|:--------|:-----|
|asm2asmAlignment.aligner|"minimap2" |
|asm2asmAlignment.alignmentBam.options |"-L --eqx --cs"|
|asm2asmAlignment.refAssemblyFastaGz | "gs://masri/flagger/v0.3.0/chm13v2.0.fa.gz" |
|asm2asmAlignment.alignmentBam.threadCount |32|
|asm2asmAlignment.preset | asm5|
|asm2asmAlignment.suffix | "chm13_v2.0" |
|asm2asmAlignment.alignmentBam.memSize | 48|

Below are the main commands for running Flagger locally using Cromwell.
```
wget https://github.com/broadinstitute/cromwell/releases/download/85/cromwell-85.jar
wget https://github.com/broadinstitute/cromwell/releases/download/85/womtool-85.jar

# Get version 0.3.0 of Flagger
wget https://github.com/mobinasri/flagger/archive/refs/tags/v0.3.0.zip

unzip v0.3.0.zip

# make a directory for saving outputs and json files
mkdir workdir 

cd workdir

java -jar ../womtool-58.jar inputs ../flagger-0.3.0/wdls/workflows/flagger_end_to_end.wdl > inputs.json
```

After modifying `inputs.json` based on the recommended parameters and the paths to input files; `assemblyFastaGz`, `fai`, `hap1ToRefBam`, `hap2ToRefBam`. and removing any other parameter from the json file you can run the command below:

```
# run flagger workflow
java -jar ../cromwell-58.jar run ../flagger-0.3.0/wdls/workflows/flagger_end_to_end.wdl -i inputs.json -m outputs.json
```
The paths to output files will be saved in `outputs.json`. The instructions for running any other WDL is similar.

It is also possible to run the pipeline in only the first mode using [flagger_end_to_end_no_variant_calling.wdl](https://github.com/mobinasri/flagger/blob/main/wdls/workflows/flagger_end_to_end_no_variant_calling.wdl), which ignores variant calling and filtering alignments.

If the assembly is related to a species without any reliable annotated reference [flagger_end_to_end_no_variant_calling_no_ref.wdl](https://github.com/mobinasri/flagger/blob/main/wdls/workflows/flagger_end_to_end_no_variant_calling_no_ref.wdl) can be used. This WDL does not need reference annotation files and the alignments to the reference assembly. It operates in the first mode which ignores variant calling and filtering alignments.

#### Dockstore links

All WDLs are uploaded to Dockstore for easier import into platforms like Terra or AnVIL.
- [Dockstore link for flagger_end_to_end.wdl](https://dockstore.org/workflows/github.com/mobinasri/flagger/FlaggerEndToEnd:v0.3.0?tab=info)

- [Dockstore link for flagger_end_to_end_no_variant_calling.wdl](https://dockstore.org/workflows/github.com/mobinasri/flagger/FlaggerEndToEndNoVariantCalling:v0.3.0?tab=info)

- [Dockstore link for flagger_end_to_end_no_variant_calling_no_ref.wdl](https://dockstore.org/workflows/github.com/mobinasri/flagger/FlaggerEndToEndNoVariantCallingNoRef:v0.3.0?tab=info)



### Components

|Component| Status| Color |Description|
|:--------|:-----|:-----|:----------|
|Err  |**Erroneous** |Red| This block has low read coverage. If it is located in the middle of a contig it could be either a misjoin or a region that needs polishing|
|Dup  |**Duplicated** |Orange| This block is potentially a false duplication of another block. It should mainly include low-MAPQ alignments with half of the expected coverage. Probably one of the copies has to be polished to fix this issue|
|Hap  | **Haploid** |Green| This block is correctly assembled and has the expected read coverage |
|Col |**Collapsed** |Purple| Two or more highly similar haplotypes are collapsed into this block |
|Unk |**Unknown** |Gray| These blocks could not be assigned confidently (usually on the edges of other components)|

Each of these components has their own color when they are shown in the IGV or the UCSC Genome Browser.


### Running Flagger on HPRC assemblies

The haplotype-resolved assemblies of the HPRC-Y1 samples and their corresponding data sets are available in

https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=working/HPRC/

For more details read this github page:

https://github.com/human-pangenomics/HPP_Year1_Assemblies

We have used the Genbank version of the HPRC-Y1 assemblies.

The v0.1 results are available in 
https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/e9ad8022-1b30-11ec-ab04-0a13c5208311--COVERAGE_ANALYSIS_Y1_GENBANK/FLAGGER/

### Publications
Liao, Wen-Wei, Mobin Asri, Jana Ebler, Daniel Doerr, Marina Haukness, Glenn Hickey, Shuangjia Lu et al. "[A draft human pangenome reference.](https://www.nature.com/articles/s41586-023-05896-x)" Nature 617, no. 7960 (2023): 312-324.
