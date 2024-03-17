# Running workflows on a Slurm system using test datasets

## Workflows Overview

This directory contains test files and tables for the workflows listed below:

1. **long_read_aligner_scattered.wdl**: This workflow is designed to align long reads to any reference or assembly which can be either diploid or haploid. Currently this workflow supports using minimap2, winnowmap and veritymap (veritymap is experimental and not tested reliably yet). Detailed descriptions of input parameters are included in the related WDL file.

2. **flagger_end_to_end.wdl**: This workflow is designed to run Flagger for evaluating a diploid assembly. Its primary inputs consist of one fasta file per assembled haplotype and one bam file containing the long read alignments to the diploid assembly.


These wdls can be found in `wdls/workflows/` including detailed descriptions of their input parameters.

Below are the guidelines for executing the workflows on the provided test datasets. These instructions can be useful for users 
who are either running these workflows for the first time or aiming to confirm their system's capability to handle the provided WDL scripts. 
They can also serve as sanity checks for future developements.

## Prerequisite softwares

### Toil 
Toil is an open-source pure-Python workflow engine. Toil is able to take workflows that are written in WDL format and execute them on systems 
with Slurm workload manager.
You can install Toil on your system using this command. More information about installing Toil is avaliable [here](https://toil.readthedocs.io/en/latest/gettingStarted/install.html#installation).
```
pip install toil[all]
```

## Steps for executing workflows

### long_read_aligner_scattered.wdl

#### 1. Cloning Flagger repository
```
git clone -b dev-0.3.0 https://github.com/mobinasri/flagger

# Go to the flagger directory
cd flagger
FLAGGER_DIR=${PWD}

# Go to the related directory for testing long_read_aligner_scattered 
cd test_wdls/toil_on_slurm/test_long_read_aligner_scattered
WORKING_DIR=${PWD}
```

#### 2. Downloading datasets
For long_read_aligner_scattered.wdl there is one set of test files (test_1). Its related data table (data_table_test_1_template.csv) 
constitutes of 7 different rows. Each row reprents a distinct combination of input parameters. For instance, they may include variations in the format of reads (e.g., bam or fastq.gz) or the specification of different mappers (e.g., minimap2 or winnowmap) .

```
## Download test_1.tar.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/e093fd72-e31a-11ee-b020-27964ee37032--flagger_test_files/flagger_v0.4.0/test_files/test_long_read_aligner_scattered/test_1.tar.gz

## Extract test_1 files
tar -xvzf test_1.tar.gz
```

```
# List files
tree test_1

test_1
├── bam_files
│   └── HG002.trio_hifiasm_0.19.5.DC_1.2.diploid.DC_1.2_40x.winnowmap_2.03.chr15_only.subsample_0.1.bam
├── fasta_files
│   ├── HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.chr15_only.fa
│   ├── HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.chr15_only.fa.fai
│   └── HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.chr15_only.fa.gz
└── fastq_files
    ├── HG002.trio_hifiasm_0.19.5.DC_1.2.diploid.DC_1.2_40x.winnowmap_2.03.chr15_only.subsample_0.1.part_1.fq.gz
    └── HG002.trio_hifiasm_0.19.5.DC_1.2.diploid.DC_1.2_40x.winnowmap_2.03.chr15_only.subsample_0.1.part_2.fq.gz
```
Description the files in test_1:
- `fasta_files` folder contains a gz-compressed fasta file that includes the contigs assembled by hifiasm_0.19.5 and subsetted to only those attributed to chromosome 15.
- `bam_files` folder contains a bam file with approximately 4x HiFi reads aligned to the chr15 contigs.
- `fastq_files` folder contains two gz-compressed fastq files including HiFi reads from chr15 with approximately 4x coverage altogether.

#### 3. Creating input json files

- Replace WORK_DIR with `$PWD`. The template data table (data_table_test_1_template.csv) does not contain actual file paths so modify the template paths to contain actual paths.
```
sed 's|WORK_DIR|'${PWD}'|g' data_table_test_1_template.csv > data_table_test_1.csv 
```

Then we need to create an input json file for Toil. We are using [launch_from_table.py](https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/main/hpc/launch_from_table.py) 
for this aim. More information about this script can be found in
[one of the HPRC githubs](https://github.com/human-pangenomics/hprc_intermediate_assembly/tree/main/hpc).

```
## Make sure you are in the working directory. Check step 1 for setting ${WORKING_DIR} if it's not set already
cd ${WORKING_DIR}

## Get the script for creating input json files.
wget https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/1f61ff0043442d8350a282ef3533def588bee8dc/hpc/launch_from_table.py

## Save WDL path and name in environment variables
WDL_PATH=${FLAGGER_DIR}/wdls/workflows/long_read_aligner_scattered.wdl
WDL_FILENAME=$(basename ${WDL_PATH})
WDL_NAME=${WDL_FILENAME%%.wdl}

## Make a folder for saving files related to run e.g. input and output jsons

mkdir -p run_test_1_toil_slurm
cd run_test_1_toil_slurm

## Make a directory for saving input json files
mkdir -p ${WDL_NAME}_input_jsons
cd ${WDL_NAME}_input_jsons

## Make input json files
python3  ${WORKING_DIR}/launch_from_table.py \
            --data_table ${WORKING_DIR}/data_table_test_1.csv \
            --field_mapping ${WORKING_DIR}/input_mapping_test_1.csv \
            --workflow_name ${WDL_NAME}
```
The output log should be
```
Creating json for HG002_hifiasm_chr15_only_test_secphase_and_md_tag
Creating json for HG002_hifiasm_chr15_only_test_secphase
Creating json for HG002_hifiasm_chr15_only_test_md_tag
Creating json for HG002_hifiasm_chr15_only_test
Creating json for HG002_hifiasm_chr15_only_test_secphase_and_md_tag_fasta_gz
Creating json for HG002_hifiasm_chr15_only_test_fasta_gz
Creating json for HG002_hifiasm_chr15_only_test_reads_bam
```

It has created 7 different input json files each of each contain a different combination of parameters. For example we take a look at the content of one json file for example
```
cat HG002_hifiasm_chr15_only_test_secphase_and_md_tag_long_read_aligner_scattered.json 
{
  "longReadAlignmentScattered.preset": "map-hifi",
  "longReadAlignmentScattered.suffix": "HiFi_minimap2",
  "longReadAlignmentScattered.alignerOptions": "-I8g --eqx -Y -L --cs",
  "longReadAlignmentScattered.readExtractionOptions": " ",
  "longReadAlignmentScattered.alignment.dockerImage": "mobinasri/long_read_aligner:v0.4.0",
  "longReadAlignmentScattered.enableSplittingReadsEqually": false,
  "longReadAlignmentScattered.splitNumber": 16,
  "longReadAlignmentScattered.alignment.threadCount": 8,
  "longReadAlignmentScattered.aligner": "minimap2",
  "longReadAlignmentScattered.kmerSize": 19,
  "longReadAlignmentScattered.sampleName": "HG002_hifiasm_chr15_only_test_secphase_and_md_tag",
  "longReadAlignmentScattered.assemblyFasta": "/private/groups/patenlab/masri/apps/test_flagger_tutorial/flagger/test_wdls/toil_on_slurm/test_long_read_aligner_scattered/test_1/fasta_files/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.chr15_only.fa",
  "longReadAlignmentScattered.alignment.memSize": 32,
  "longReadAlignmentScattered.readFiles": [
    "/private/groups/patenlab/masri/apps/test_flagger_tutorial/flagger/test_wdls/toil_on_slurm/test_long_read_aligner_scattered/test_1/fastq_files/HG002.trio_hifiasm_0.19.5.DC_1.2.diploid.DC_1.2_40x.winnowmap_2.03.chr15_only.subsample_0.1.part_1.fq.gz",
    "/private/groups/patenlab/masri/apps/test_flagger_tutorial/flagger/test_wdls/toil_on_slurm/test_long_read_aligner_scattered/test_1/fastq_files/HG002.trio_hifiasm_0.19.5.DC_1.2.diploid.DC_1.2_40x.winnowmap_2.03.chr15_only.subsample_0.1.part_2.fq.gz"
  ],
  "longReadAlignmentScattered.enableAddingMDTag": true,
  "longReadAlignmentScattered.enableRunningSecphase": true
```

In this json we are using `minimap2` aligner with the parameter preset of `map-hifi` and kmer size of 19 for aligning hifi reads. `longReadAlignmentScattered.readFiles` points to two read files in the `fq.gz` format. Both `longReadAlignmentScattered.enableAddingMDTag` and `longReadAlignmentScattered.enableRunningSecphase` are true which means that the pipeline will add MD tags to the final bam file and also run secphase for correcting potentially wrong alignments. For more information about other parameters take a look at the WDL file.

#### 4. Executing workflow using input json files
For running this WDL on Slurm we are using [a bash script](https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/2e5155690ec365e906dc82e72be39014dc38de27/hpc/toil_sbatch_single_machine.sh) that can execute an array of jobs by taking the data table csv file. For each row in the csv file Toil will create a separate job after acquiring the speficied cpu `--cpus-per-task` and memory `--mem`.
```
## Make sure you are in the working directory. Check step 1 for setting ${WORKING_DIR} if it's not set already
cd ${WORKING_DIR}

## Get the bash script for running WDLs on Slurm using Toil
wget https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/2e5155690ec365e906dc82e72be39014dc38de27/hpc/toil_sbatch_single_machine.sh

## Create folders for saving output json files
mkdir -p run_test_1_toil_slurm/${WDL_NAME}_output_jsons
mkdir -p run_test_1_toil_slurm/${WDL_NAME}_logs

## Set environment variables for sbatch
USERNAME="your_user_name"
EMAIL="your@email"
TIME_LIMIT="5:00:00"
# Partition should be modifed based on the available partitions on the server
PARTITION="medium"

## Go to the execution directory
cd run_test_1_toil_slurm

## Run jobs arrays
## --array=1-7%7 will make 7 jobs; one per json file (ordered by rows in csv file)
## It can be modified based on the desired jobs
## For example --array=4-5%2 will run the jobs related to the 4th and 5th rows in the csv file
## or --array=7-7%1 will run only the job related to the 7th row
sbatch      --job-name=${WDL_NAME}_${USERNAME} \
            --cpus-per-task=32 \
            --mem=64G \
            --mail-user=${EMAIL} \
            --output=${WDL_NAME}_logs/${WDL_NAME}_%A_%a.log \
            --array=7-7%1  \
            --time=${TIME_LIMIT} \
            --partition=${PARTITION} \
            ${WORKING_DIR}/toil_sbatch_single_machine.sh \
            --wdl ${WDL_PATH} \
            --sample_csv  ${WORKING_DIR}/data_table_test_1.csv \
            --input_json_path ${WORKING_DIR}/run_test_1_toil_slurm/${WDL_NAME}_input_jsons/\${SAMPLE_ID}_${WDL_NAME}.json
```

#### flagger_end_to_end.wdl
For flagger_end_to_end.wdl there are two set of test files (test_1 and test_2). Its related data table (data_table_test_1_template.csv) contain 7 
different rows. Each row contain different combinations of input parameters (for example reads can be in bam or fastq.gz format or 
different mappers may be specified)
```
cd test_long_read_aligner_scattered
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/e093fd72-e31a-11ee-b020-27964ee37032--flagger_test_files/flagger_v0.4.0/test_files/test_long_read_aligner_scattered/test_1.tar.gz
```


