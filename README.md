## Evaluating dual/diploid assemblies with Flagger

### Overview
Here is a description of a read-based pipeline that can detect different types of mis-assemblies in a draft dual/diploid assembly. (*What is a dual assembly? Read [this page](https://lh3.github.io/2021/10/10/introducing-dual-assembly)*). The core component of this pipeline is [**Flagger**](https://github.com/mobinasri/flagger/tree/main/docs/flagger). Flagger recieves the read alignments to a draft diploid assembly, detects the anomalies in the read coverage along the assembly and partition the assembly into 4 main components; erroneous, (falsely) duplicated, haploid (correctly assembled) and collapsed.



This evaluation has 5 steps:
- Align long reads to the diploid assembly
- Phase and relocalize the reads with secondary alignments using [secphase](https://github.com/mobinasri/secphase) (Optional)
- Run Flagger

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

### 3. Run Flagger
The produced alignment file (`${INPUT_BAM}` in step 1) can be used as the input to Flagger. Flagger outputs a bed file with 5 labels; 
erroneous (Err), duplicated (Dup), haploid (Hap), collapsed (Col) and unkown (Unk). Any component other than the haploid one is pointing to unreliable blocks in assembly and unkown label is for the bases that couldn't be assigned confidently. The 4 components are explained in detail [here](https://github.com/mobinasri/flagger/tree/main/docs/coverage#2-coverage-distribution-and-fitting-the-mixture-model). 


More information about Flagger is available [here](https://github.com/mobinasri/flagger/tree/main/docs/flagger)

#### Running pipeline with WDL

It is easier to run the pipeline using the WDLs described below. A WDL file can be run locally using Cromwell, which is an open-source Workflow Management System for bioinformatics. The latest releases of Cromwell are available [here](https://github.com/broadinstitute/cromwell/releases) and the documentation is available [here](https://cromwell.readthedocs.io/en/stable/CommandLine/).

It is recommended to run the whole pipeline using [flagger_end_to_end_with_mapping.wdl](https://github.com/mobinasri/flagger/blob/dev-0.3.0/wdls/workflows/flagger_end_to_end_with_mapping.wdl).

Recommended values for the parameters of flagger_end_to_end_with_mapping.wdl:


| Parameter | Description | Type | Default | 
| --- | --- | --- | ------------ |
|sampleName| Sample name; for example 'HG002'| String | No Default |
|suffixForFlagger| Suffix string that contains information about this analysis; for example 'hifi_winnowmap_flagger_for_hprc'| String | No Default |
|suffixForMapping| Suffix string that contains information about this alignment. It will be appended to the name of the final alignment. For example 'hifi_winnowmap_v2.03_hprc_y2'| String | No Default |
|hap1AssemblyFasta| Path to uncompressed or gzip-compressed fasta file of the 1st haplotype.| File | No Default |
|hap2AssemblyFasta| Path to uncompressed or gzip-compressed fasta file of the 2nd haplotype.| File | No Default |
|readfiles| An array of read files. Their format can be either fastq, fq, fastq.gz, fq.gz, bam or cram. For cram format referenceFastaForReadExtraction should also be passed.| Array[File] | No Default |
|aligner| Name of the aligner. It can be either minimap2, winnowmap or veritymap.| String | winnowmap |
|preset| Paremeter preset should be selected based on aligner and sequencing platform. Common presets are map-pb/map-hifi/map-ont for minimap2, map-pb/map-ont for winnowmap and hifi-haploid/hifi-haploid-complete/hifi-diploid/ont-haploid-complete for veritymap| String | No Default |
|kmerSize| The kmer size for using minimap2 or winnowmap. With winnowmap kmer size should be 15 and with minimap2 kmer size should be 17 and 19 for using the presets map-ont and map-hifi/map-pb respectively.| Int | 15 |
|alignerOptions| Aligner options. It can be something like '--eqx --cs -Y -L -y' for minimap2/winnowmap. Note that if assembly is diploid and aligner is either minimap2 or winnowmap '-I8g' is necessary. If the reads contain modification tags and these tags are supposed to be present in the final alignment file, alignerOptions should contain '-y' and the aligner should be either minimap2 or winnowmap. If running secphase is enabled it is recommended to add '-p0.5' to alignerOptions; it will keep more secondary alignments so secphase will have more candidates per read. For veritymap '--careful' can be used but not recommended for whole-genome assembly since it increases the runtime dramatically.| String | --eqx -Y -L -y |
|readExtractionOptions| The options to be used while converting bam to fastq with samtools fastq. If the reads contain epigentic modification tags it is necessary to use '-TMm,Ml,MM,ML'| String | -TMm,Ml,MM,ML |
|referenceFastaForReadExtraction| If reads are in CRAM format then the related reference should be passed to this parameter for read extraction.| File | No Default (Optional) |
|enableAddingMDTag| If true it will call samtools calmd to add MD tags to the final bam file.| Boolean | true |
|splitNumber| The number of chunks which the input reads should be equally split into. Note that enableSplittingReadsEqually should be set to true if user wants to split reads into equally sized chunks.| Int | 16 |
|enableSplittingReadsEqually| If true it will merge all reads together and then split them into multiple chunks of roughly equal size. Each chunk will then be aligned via a separate task. This feature is useful for running alignment on cloud/slurm systems where there are  multiple nodes available with enough computing power and having alignments tasks distributed among small nodes is more efficient or cheaper than running a single alignment task in a large node. If the  whole workflow is being on a single node it is not recommened to use this feature since mergin and splitting reads takes its own time. | Boolean | false |
|minReadLengthForMapping| If it is greater than zero, a task will be executed for filtering reads shorter than this value before alignment.| Int | 0 |
|alignerThreadCount | The number of threads for mapping in each alignment task | Int | 16 |
|alignerMemSize | The size of the memory in Gb for mapping in each alignment task | Int | 48 |
|alignerDockerImage | The mapping docker image | String | mobinasri/long_read_aligner:v0.4.0 | 
|correctBamOptions| Options for the correct_bam program that can filters short/highly divergent alignments  | String | --primaryOnly --minReadLen 5000 --minAlignment 5000 --maxDiv 0.1 | 
|preemptible| Number of retries to use preemptible nodes on Terra/GCP. | Int | 2 |
|zones| Name of the zone for taking nodes on Terra/GCP. | String | us-west2-a| 
|maxReadDivergenceForFlagger| Alignments with gap-compressed ratio higher than this will be filtered in the pre-process step of flagger. | Float | 0.1 |
|potentialBiasesBedArray| Array of bed files each of which contains regions with potential coverage bias for example one bed file can contain HSat2 regions in haplotype 1. | Array[File] | No Default (Optional) |
|sexBed| Optional bed file containing regions assigned to X/Y chromosomes. (can be either in ref or asm coordinates)| File | No Default (Optional) |
|SDBed| Optional Bed file containing Segmental Duplications. (can be either in ref or asm coordinates)| File | No Default (Optional) |
|cntrBed| Optional Bed file containing peri/centromeric satellites (ASat, HSat, bSat, gSat) without 'ct' blocks.| File | No Default (Optional) |
|cntrCtBed| Optional Bed file containing centromere transition 'ct' blocks.| File | No Default (Optional) | 
|additionalStratificationBedArray| Array of additional stratification bed files for final stats tsv file. | Array[File] | No Default (Optional) |
|additionalStratificationNameArray| Array of names for the stratifications provided in the argument additionalStratificationBedArray. |  Array[File] | No Default (Optional) |
|enableProjectingBedsFromRef2Asm| If True it means that the given bed files are in ref coors (e.g. chm13v2) and they have to be projected to asm coors. | Boolean | false |
|projectionReferenceFastaGz| The given bed files are in the coordinates of this reference. A reference should be passed if enableProjectingBedsFromRef2Asm is true. | File | No Default (Optional) |
|enableRunningSecphase | If true it will run secphase in the marker mode using the wdl parameters starting with 'secphase' otherwise skip it. | Boolean | false |
|secphaseDockerImage| Docker image for running Secphase | String | mobinasri/secphase:v0.4.3 |
|secphaseOptions| String containing secphase options (can be either --hifi or --ont). | String | --hifi |
|secphaseVersion| Secphase version.| String | v0.4.3
|enableOutputtingWig| If true it will make wig files from cov files and output them. wig files can be easily imported into IGV sessions | Boolean | true |
|enableOutputtingBam| If true it will output read alignment bam file and its related index | Boolean | false |
|windowSize| The size of the window flagger uses for finding coverage distrubutions (Default = 5Mb)| Int | 5000000 | 
|sortPdfPagesByHaplotype| Sort the coverage distributions plotted in the output pdf by haplotype | Boolean | false |
|hap1ContigPattern| The pattern that will be used for finding the names of the contigs belonging to haplotype1. It will be skipped if sortPdfPagesByHaplotype is false. | String | hap1 |
|hap2ContigPattern| The pattern that will be used for finding the names of the contigs belonging to haplotype2. It will be skipped if sortPdfPagesByHaplotype is false. | String | hap2 |


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
wget https://github.com/mobinasri/flagger/archive/refs/tags/v0.4.0.zip

unzip v0.4.0.zip

# make a directory for saving outputs and json files
mkdir workdir 

cd workdir

java -jar ../womtool-58.jar inputs ../flagger-0.4.0/wdls/workflows/flagger_end_to_end.wdl > inputs.json
```

After modifying `inputs.json` based on the recommended parameters and the paths to input files; `assemblyFastaGz`, `fai`, `hap1ToRefBam`, `hap2ToRefBam`. and removing any other parameter from the json file you can run the command below:

```
# run flagger workflow
java -jar ../cromwell-58.jar run ../flagger-0.4.0/wdls/workflows/flagger_end_to_end.wdl -i inputs.json -m outputs.json
```
The paths to output files will be saved in `outputs.json`. The instructions for running any other WDL is similar.


#### Dockstore links

All WDLs are uploaded to Dockstore for easier import into platforms like Terra or AnVIL.
- [Dockstore link for flagger_end_to_end.wdl](https://dockstore.org/workflows/github.com/mobinasri/flagger/FlaggerEndToEndWithMapping:v0.4.0?tab=info)



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
