## Evaluating dual/diploid assemblies with HMM-Flagger

### ***Note: HMM-Flagger is an HMM-based version of Flagger, offering improved functionality and performance. It replaces the previous version, Flagger v0.4.0, and is recommended for all users from now on.***
### Overview

HMM-Flagger is a read-mapping-based tool that can detect different types of mis-assemblies in a dual or diploid genome assembly. HMM-Flagger recieves the read alignments to a genome assembly, uses Hidden Markov Model to detect anomalies in the read coverage along the assembly and finally partitions the assembly into four classes; erroneous, falsely duplicated, haploid (structurally correct) and collapsed.

## Quick Start In Three Steps (Needs a BAM and FASTA file)
### 1. Create a whole-genome BED file

```
# Go to the working directory
# Put fasta, bam and bam index files in this directory
cd ${WORKING_DIR} 


# Create fasta index assuming that file is not gz-compressed
samtools faidx ${FASTA_FILE}


# Create a bed file covering whole genome
cat ${FASTA_FILE}.fai | \
    awk '{print $1"\t0\t"$2}' > whole_genome.bed
```


### 2. Convert BAM to COV

```
# Put path to the whole-genome bed file in a json file
echo "{" > annotations_path.json
echo \"whole_genome\" : \"${PWD}/whole_genome.bed\" >> annotations_path.json
echo "}" >> annotations_path.json


# Convert bam to cov.gz with bam2cov program
docker run -it --rm -v${WORKING_DIR}:${WORKING_DIR} mobinasri/flagger:v1.0.0-prerelease \
  bam2cov --bam ${WORKING_DIR}/${BAM_FILE} \
                  --output ${WORKING_DIR}/coverage_file.cov.gz \
                  --annotationJson ${WORKING_DIR}/annotations_path.json \
                  --threads 16 \
                  --baselineAnnotation whole_genome
```

### 3. Run HMM-Flagger

```
mkdir -p ${WORKING_DIR}/hmm_flagger_outputs
docker run -it --rm -v${WORKING_DIR}:${WORKING_DIR} mobinasri/flagger:v1.0.0-prerelease \
        hmm_flagger \
            --input ${WORKING_DIR}/coverage_file.cov.gz \
            --outputDir ${WORKING_DIR}/hmm_flagger_outputs  \
            --alphaTsv /home/programs/config/alpha_optimum_trunc_exp_gaussian_w_4000_n_50.tsv \
            --labelNames Err,Dup,Hap,Col \
            --threads 16
```

The bed file contating HMM-Flagger predictions will be located in `hmm_flagger_outputs/final_flagger_prediction.bed`. A summary tsv file (counts will be based on windows rather than bases) will be located in `hmm_flagger_outputs/prediction_summary_final.tsv`

For example this table shows that `95.99%` of the assembly is flagged as correctly assembled:
```
cat hmm_flagger_outputs/prediction_summary_final.tsv | grep base_level | grep whole_genome | grep percentage
PREDICTION	base_level	percentage	annotation	whole_genome	ALL_SIZES	ALL_LABELS	1.24	1.37	95.99	1.40	0.00
```

### Note about coverage-biased regions
It has been observed that peri/centromeric satellite regions are prone to have coverage biases (in PacBio HiFi and ONT data). Users can add satellite repeat arrays as separate annotations to the input json for the `bam2cov` program (with enabling `--runBiasDetection`). `bam2cov` will detect if there is coverage biases in any annotation provided in the json file and add the neccessary information to the created coverage file to help HMM-Flagger predict the labels more accurately. It is recommended to provide each satellite in a separate bed file for example putting HSat1A and Hsat2 in separate files.

## Running pipeline with WDL

If users have a set of unalinged reads it is easier to run the pipeline (read mapping + HMM-Flagger) using the WDLs described below. Using WDLs it is also easier to incorporate coverage biases that may exist in the satellite regions. 

A WDL file can be run locally using Cromwell, which is an open-source Workflow Management System for bioinformatics. The latest releases of Cromwell are available [here](https://github.com/broadinstitute/cromwell/releases) and the documentation is available [here](https://cromwell.readthedocs.io/en/stable/CommandLine/).

It is recommended to run the whole pipeline using [hmm_flagger_end_to_end_with_mapping.wdl](https://github.com/mobinasri/flagger/blob/v1.0.0/wdls/workflows/hmm_flagger_end_to_end_with_mapping.wdl).

Here is a list of input parameters for hmm_flagger_end_to_end_with_mapping.wdl (The parameters marked as **"(Mandatory)"** are mandatory to be defined in the input json):


| Parameter | Description | Type | Default | 
| --- | --- | --- | ------------ |
|sampleName| Sample name; for example 'HG002'| String | **No Default (Mandatory)** |
|suffixForFlagger| Suffix string that contains information about this analysis; for example 'hifi_winnowmap_flagger_for_hprc'| String | **No Default (Mandatory)** |
|suffixForMapping| Suffix string that contains information about this alignment. It will be appended to the name of the final alignment file. For example 'hifi_winnowmap_v2.03_hprc_y2'| String | **No Default (Mandatory)** |
|hap1AssemblyFasta| Path to uncompressed or gzip-compressed fasta file of the 1st haplotype.| File | **No Default (Mandatory)** |
|hap2AssemblyFasta| Path to uncompressed or gzip-compressed fasta file of the 2nd haplotype.| File | **No Default (Mandatory)** |
|readfiles| An array of read files. Their format can be either fastq, fq, fastq.gz, fq.gz, bam or cram. For cram format referenceFastaForReadExtraction should also be passed.| Array[File] | **No Default (Mandatory)** |
|preset| Paremeter preset should be selected based on aligner and sequencing platform. Common presets are map-pb/map-hifi/map-ont for minimap2, map-pb/map-ont for winnowmap and hifi-haploid/hifi-haploid-complete/hifi-diploid/ont-haploid-complete for veritymap| String | **No Default (Mandatory)** |
|aligner| Name of the aligner. It can be either minimap2, winnowmap or veritymap.| String | winnowmap |
|kmerSize| The kmer size for using minimap2 or winnowmap. With winnowmap kmer size should be 15 and with minimap2 kmer size should be 17 and 19 for using the presets map-ont and map-hifi/map-pb respectively.| Int | 15 |
|alignerOptions| Aligner options. It can be something like '--eqx --cs -Y -L -y' for minimap2/winnowmap. Note that if assembly is diploid and aligner is either minimap2 or winnowmap '-I8g' is necessary. If the reads contain modification tags and these tags are supposed to be present in the final alignment file, alignerOptions should contain '-y' and the aligner should be either minimap2 or winnowmap. If running secphase is enabled it is recommended to add '-p0.5' to alignerOptions; it will keep more secondary alignments so secphase will have more candidates per read. For veritymap '--careful' can be used but not recommended for whole-genome assembly since it increases the runtime dramatically.| String | --eqx -Y -L -y |
|readExtractionOptions| The options to be used while converting bam to fastq with samtools fastq. If the reads contain epigentic modification tags it is necessary to use '-TMm,Ml,MM,ML'| String | -TMm,Ml,MM,ML |
|referenceFastaForReadExtraction| If reads are in CRAM format then the related reference should be passed to this parameter for read extraction.| File | No Default (Optional) |
|enableAddingMDTag| If true it will call samtools calmd to add MD tags to the final bam file.| Boolean | true |
|splitNumber| The number of chunks which the input reads should be equally split into. Note that enableSplittingReadsEqually should be set to true if user wants to split reads into equally sized chunks.| Int | 16 |
|enableSplittingReadsEqually| If true it will merge all reads together and then split them into multiple chunks of roughly equal size. Each chunk will then be aligned via a separate task. This feature is useful for running alignment on cloud/slurm systems where there are  multiple nodes available with enough computing power and having alignments tasks distributed among small nodes is more efficient or cheaper than running a single alignment task in a large node. If the  whole workflow is being on a single node it is not recommened to use this feature since merging and splitting reads takes its own time. | Boolean | false |
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
|projectionReferenceFasta| The given bed files are in the coordinates of this reference. A reference should be passed if enableProjectingBedsFromRef2Asm is true. | File | No Default (Optional) |
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

### Using CHM13 annotation files (Optional)

A set of annotations files in the coordinates of chm13v2.0 are prepared beforehand and they can be used as inputs to the workflow. These bed files are useful when there is no denovo annotation for the assemblies and we like to project annotations from chm13v2.0 to the assembly coordinates for two main purposes:

1. Detecting regions with coverage biases: Satellite repeat arrays might have coverage biases so before running HMM-Flagger it is best to pass those annotations to bam2cov program (with `--runBiasDetection` enabled). bam2cov will inspect those regions and put this information inside the created coverage file if any coverage-biased region was detected.
2. Stratifying final results with the projected annotations.

Using these bed files are optional and the related parameters can be left undefined (being absent from the input json). In this case the workflow should still work properly. However users should be aware that in this case potential coverage biases in HSat arrays may mislead HMM-Flagger.

|Parameter| Value|
|:--------|:-----|
|enableProjectingBedsFromRef2Asm|true|
|projectionReferenceFastaGz| https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz|
|potentialBiasesBedArray| ['https://raw.githubusercontent.com/mobinasri/flagger/v0.4.0/misc/potential_biases/chm13v2.0_bsat.bed','https://raw.githubusercontent.com/mobinasri/flagger/v0.4.0/misc/potential_biases/chm13v2.0_hsat1A.bed','https://raw.githubusercontent.com/mobinasri/flagger/v0.4.0/misc/potential_biases/chm13v2.0_hsat1B.bed','https://raw.githubusercontent.com/mobinasri/flagger/v0.4.0/misc/potential_biases/chm13v2.0_hsat2.bed','https://raw.githubusercontent.com/mobinasri/flagger/v0.4.0/misc/potential_biases/chm13v2.0_hsat3.bed','https://raw.githubusercontent.com/mobinasri/flagger/v0.4.0/misc/potential_biases/chm13v2.0_hor.bed'] |
|censat_bed|https://raw.githubusercontent.com/mobinasri/flagger/v0.4.0/misc/stratifications/censat/chm13v2.0_no_ct.bed|
|sd_bed|https://raw.githubusercontent.com/mobinasri/flagger/v0.4.0/misc/stratifications/sd/chm13v2.0_SD.all.bed|
|sex_bed|https://raw.githubusercontent.com/mobinasri/flagger/v0.4.0/misc/stratifications/sex/chm13v2.0_sex.bed|
|additional_bed_array| ['https://raw.githubusercontent.com/mobinasri/flagger/v0.4.0/misc/stratifications/censat/chm13v2.0_no_ct.bed','https://raw.githubusercontent.com/mobinasri/flagger/v0.4.0/misc/stratifications/censat/chm13v2.0_only_ct.bed','https://raw.githubusercontent.com/mobinasri/flagger/v0.4.0/misc/stratifications/censat/chm13v2.0_bsat.bed','https://raw.githubusercontent.com/mobinasri/flagger/v0.4.0/misc/stratifications/censat/chm13v2.0_gsat.bed','https://raw.githubusercontent.com/mobinasri/flagger/v0.4.0/misc/stratifications/censat/chm13v2.0_hor.bed','https://raw.githubusercontent.com/mobinasri/flagger/v0.4.0/misc/stratifications/censat/chm13v2.0_hsat1A.bed','https://raw.githubusercontent.com/mobinasri/flagger/v0.4.0/misc/stratifications/censat/chm13v2.0_hsat1B.bed','https://raw.githubusercontent.com/mobinasri/flagger/v0.4.0/misc/stratifications/censat/chm13v2.0_hsat2.bed','https://raw.githubusercontent.com/mobinasri/flagger/v0.4.0/misc/stratifications/censat/chm13v2.0_hsat3.bed','https://raw.githubusercontent.com/mobinasri/flagger/v0.4.0/misc/stratifications/censat/chm13v2.0_mon.bed','https://raw.githubusercontent.com/mobinasri/flagger/v0.4.0/misc/stratifications/sd/chm13v2.0_SD.g99.bed','https://raw.githubusercontent.com/mobinasri/flagger/v0.4.0/misc/stratifications/sd/chm13v2.0_SD.g98_le99.bed','https://raw.githubusercontent.com/mobinasri/flagger/v0.4.0/misc/stratifications/sd/chm13v2.0_SD.g90_le98.bed','https://raw.githubusercontent.com/mobinasri/flagger/v0.4.0/misc/stratifications/sd/chm13v2.0_SD.le90.bed','https://raw.githubusercontent.com/mobinasri/flagger/v0.4.0/misc/stratifications/sd/chm13v2.0_SD.all.bed','https://raw.githubusercontent.com/mobinasri/flagger/v0.4.0/misc/stratifications/repeat_masker/chm13v2.0_RM_4.1.2p1_le6_STR.bed','https://raw.githubusercontent.com/mobinasri/flagger/v0.4.0/misc/stratifications/repeat_masker/chm13v2.0_RM_4.1.2p1_ge7_VNTR.bed'] |
|additional_name_array|['sat_no_ct','only_ct','bsat','gsat','hor','hsat1A','hsat1B','hsat2','hsat3','mon','sd_g99','sd_g98_le99','sd_g90_le98','sd_le90','sd_all','STR','VNTR']|

All files with git and s3 urls are publicly accessible so if you are running the WDL on Terra/AnVIL platforms you can use the same urls. They are also available in the directories `misc/annotations` and `misc/biased_regions` of this repository for those who want to run locally. 

### Other major workflows

[flagger_end_to_end_with_mapping.wdl](https://github.com/mobinasri/flagger/blob/v0.4.0/wdls/workflows/flagger_end_to_end_with_mapping.wdl) is running two major workflows; [long_read_aligner_scattered.wdl](https://github.com/mobinasri/flagger/blob/v0.4.0/wdls/workflows/long_read_aligner_scattered.wdl) and [flagger_end_to_end.wdl](https://github.com/mobinasri/flagger/blob/v0.4.0/wdls/workflows/flagger_end_to_end.wdl). The first one runs read mapping and the second one runs the end-to-end flagger pipeline (including annotation projection, bias detection, secphase and flagger). 

Users can run each of them separately. For example if there is a read alignment file available beforehand users can run only [flagger_end_to_end.wdl](https://github.com/mobinasri/flagger/blob/v0.4.0/wdls/workflows/flagger_end_to_end.wdl). The parameters for each of these workflows is a subset of the parameters listed above for [flagger_end_to_end_with_mapping.wdl](https://github.com/mobinasri/flagger/blob/v0.4.0/wdls/workflows/flagger_end_to_end_with_mapping.wdl) therefore this table can still be used as a reference for either long_read_aligner_scattered.wdl or flagger_end_to_end.wdl.


### Running WDLs locally using Cromwell

Below are the main commands for running flagger_end_to_end_with_mapping.wdl locally using Cromwell.
```
wget https://github.com/broadinstitute/cromwell/releases/download/85/cromwell-85.jar
wget https://github.com/broadinstitute/cromwell/releases/download/85/womtool-85.jar

# Get version 0.4.0 of Flagger
wget https://github.com/mobinasri/flagger/archive/refs/tags/v0.4.0.zip

unzip v0.4.0.zip


# make a directory for saving outputs and json files
mkdir workdir 

cd workdir


java -jar ../womtool-85.jar inputs ../flagger-0.4.0/wdls/workflows/flagger_end_to_end_with_mapping.wdl > inputs.json

```

After modifying `inputs.json`, setting mandatory parameters: `sampleName`, `suffixForFlagger`, `suffixForMapping`, `hap1AssemblyFasta`, `hap2AssemblyFasta`, `readfiles`, `preset` and removing unspecified parameters (they will be set to default values), users can run the command below:


```
## run flagger workflow
java -jar ../cromwell-85.jar run ../flagger-0.4.0/wdls/workflows/flagger_end_to_end_with_mapping.wdl -i inputs.json -m outputs.json
```

The paths to output files will be saved in `outputs.json`. The instructions for running any other WDL is similar.

### Running WDLs on Slurm using Toil
Instructions for running WDLs on Slurm are provided [here](https://github.com/mobinasri/flagger/tree/v1.0.0/test_wdls/toil_on_slurm) , which includes some test data sets for each of the workflows:
- [long_read_aligner_scattered.wdl](https://github.com/mobinasri/flagger/tree/v1.0.0/test_wdls/toil_on_slurm#running-long_read_aligner_scatteredwdl-on-test-datasets)
- [hmm_flagger_end_to_end.wdl](https://github.com/mobinasri/flagger/tree/v1.0.0/test_wdls/toil_on_slurm#running-flagger_end_to_endwdl-on-test-datasets)
- [hmm_flagger_end_to_end_with_mapping.wdl](https://github.com/mobinasri/flagger/tree/v1.0.0/test_wdls/toil_on_slurm#running-flagger_end_to_end_with_mappingwdl-on-test-datasets)


### Dockstore links

All WDLs are uploaded to Dockstore for easier import into platforms like Terra or AnVIL.

- [Dockstore link for long_read_aligner_scattered.wdl](https://dockstore.org/workflows/github.com/mobinasri/flagger/LongReadAlignerScattered:v1.0.0?tab=info)
- [Dockstore link for hmm_flagger_end_to_end.wdl](https://dockstore.org/workflows/github.com/mobinasri/flagger/HMMFlaggerEndToEnd:v1.0.0?tab=info)
- [Dockstore link for hmm_flagger_end_to_end_with_mapping.wdl](https://dockstore.org/workflows/github.com/mobinasri/flagger/HMMFlaggerEndToEndWithMapping:v1.0.0?tab=info)

### Components

|Component| Status| Color |Description|
|:--------|:-----|:-----|:----------|
|Err  |**Erroneous** |Red| This block has low read coverage. If it is located in the middle of a contig it could be either a misjoin or a region that needs polishing|
|Dup  |**Duplicated** |Orange| This block is potentially a false duplication of another block. It should mainly include low-MAPQ alignments with half of the expected coverage. Probably one of the copies has to be polished or removed to fix this issue|
|Hap  | **Haploid** |Green| This block is correctly assembled and has the expected read coverage |
|Col |**Collapsed** |Purple| Two or more highly similar haplotypes are collapsed into this block |

Each of these components has their own color when they are shown in the IGV or the UCSC Genome Browser.

## Publications
Liao, Wen-Wei, Mobin Asri, Jana Ebler, Daniel Doerr, Marina Haukness, Glenn Hickey, Shuangjia Lu et al. "[A draft human pangenome reference.](https://www.nature.com/articles/s41586-023-05896-x)" Nature 617, no. 7960 (2023): 312-324.
