# Running workflows on a Slurm system using test data sets

## Workflows Overview

In this directory test files and tables are provided for the workflows listed below:

### 1. long_read_aligner_scattered.wdl
This workflow is for aligning long reads to any reference or assembly which can be either diploid or haploid. Currently this workflow
supports using minimap2, winnowmap and veritymap (veritymap is experimental and not tested reliably yet). Detailed descriptions of 
input parameters are included in the related WDL file.

### 2. flagger_end_to_end.wdl
This workflow is for running Flagger for evaluating a diploid assembly. Its main inputs are one fasta file per assembled haplotype and the long 
read alignments to the diploid assembly.


These wdls can be found in `wdls/workflows/` containing detailed descriptions of their input parameters.

Below are the guidelines for executing the workflows on the provided test datasets. These instructions can be useful for users 
who are either running these workflows for the first time or aiming to confirm their system's capability to handle the provided WDL scripts. 
They can also serve as sanity checks for future developements.

## Prerequisite Steps

### Installing Toil 
Toil is an open-source pure-Python workflow engine. Toil is able to take workflows that are written in WDL format and execute them on systems 
with Slurm workload manager.
You can install Toil on your system using this command. More information about installing Toil is avaliable [here](https://toil.readthedocs.io/en/latest/gettingStarted/install.html#installation).
```
pip install toil[all]
```

## Steps for executing workflows

### long_read_aligner_scattered.wdl

#### 1. Clone Flagger repo
```

git clone -b dev-0.3.0 https://github.com/mobinasri/flagger

# Go to the flagger directory
cd flagger
FLAGGER_DIR=${PWD}

# Go to the related directory for testing long_read_aligner_scattered 
cd test_wdls/toil_on_slurm/test_long_read_aligner_scattered
WORKING_DIR=${PWD}

```

#### 2. Download datasets
For long_read_aligner_scattered.wdl there is one set of test files (test_1). Its related data table (data_table_test_1_template.csv) 
constitutes of 7 different rows. Each row contains a different combination of input parameters (for example reads can be in bam or fastq.gz 
format or different mappers may be specified)
```
# Download test_1.tar.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/e093fd72-e31a-11ee-b020-27964ee37032--flagger_test_files/flagger_v0.4.0/test_files/test_long_read_aligner_scattered/test_1.tar.gz

# Extract test_1 files
tar -xvzf test_1.tar.gz
```

```
# List the files
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
List of the files in test_1:
- `HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.chr15_only.fa.gz` is a gz-compressed fasta file that contains the contigs assembled by hifiasm_0.19.5 and attributed to chromosome 15.
- `HG002.trio_hifiasm_0.19.5.DC_1.2.diploid.DC_1.2_40x.winnowmap_2.03.chr15_only.subsample_0.1.bam` is a bam file that contains 
~4x HiFi reads aligned to the chr15 contigs.
- `HG002.trio_hifiasm_0.19.5.DC_1.2.diploid.DC_1.2_40x.winnowmap_2.03.chr15_only.subsample_0.1.part_1.fq.gz` and 
`HG002.trio_hifiasm_0.19.5.DC_1.2.diploid.DC_1.2_40x.winnowmap_2.03.chr15_only.subsample_0.1.part_2.fq.gz` are two gz-compressed
fastq files containing HiFi reads from chr15 with ~4x coverage altogether.

#### 3. Make input json files

The template data table does not contain actual file paths so first we need to modify the template paths to contain actual paths.
```
# Replace WORK_DIR with $PWD
sed 's|WORK_DIR|'${PWD}'|g' data_table_test_1_template.csv
```

Then we need to create an input json file for Toil. We are using [launch_from_table.py](https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/main/hpc/launch_from_table.py) 
for this aim. More information about this script can be found in
[one of the HPRC githubs](https://github.com/human-pangenomics/hprc_intermediate_assembly/tree/main/hpc).

```
# Get the script
wget https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/1f61ff0043442d8350a282ef3533def588bee8dc/hpc/launch_from_table.py

## Set environment variables
# Get WDL path and name
WDL_PATH=${FLAGGER_DIR}/wdls/workflows/long_read_aligner_scattered.wdl
WDL_FILENAME=$(basename ${WDL_PATH)})
WDL_NAME=${WDL_FILENAME}
INPUT_MAPPING_CSV=${PWD}/input_mapping_test_1.csv
INPUT_DATA_TABLE_CSV=${PWD}/data_table_test_1.csv
WORKING_DIR=${PWD}/run_test_1_toil_slurm

# make directories for saving input and output json files
mkdir -p ${WORKING_DIR}/${WDL_NAME}_input_jsons
cd ${WORKING_DIR}/${WDL_NAME}_input_jsons

python3  ../launch_from_table.py \
            --data_table data_table_test_1.csv \
            --field_mapping input_mapping_test_1.csv \
            --workflow_name ${WDL_NAME}
```

#### 4. Execute workflow with input json files

```
cd ${WORKING_DIR_1}
cd ${WDL_NAME}_output_jsons
mkdir -p ${WDL_NAME}_logs
sbatch      --job-name=${WDL_NAME}_${USERNAME} \
            --cpus-per-task=32 \
            --mem=64G \
            --mail-user=${EMAIL} \
            --output=${WDL_NAME}_logs/${WDL_NAME}_%A_%a.log \
            --array=7-7%1  \
            --time=${TIME_LIMIT} \
            --partition=medium \
            ${LAUNCH_WORKFLOW_JOB_ARRAY_BASH} \
            --wdl ${WDL_PATH} \
            --sample_csv  ${INPUT_DATA_TABLE_CSV_1} \
            --input_json_path ${WORKING_DIR_1}/${WDL_NAME}_input_jsons/\${SAMPLE_ID}_${WDL_NAME}.json
```

#### flagger_end_to_end.wdl
For flagger_end_to_end.wdl there are two set of test files (test_1 and test_2). Its related data table (data_table_test_1_template.csv) contain 7 
different rows. Each row contain different combinations of input parameters (for example reads can be in bam or fastq.gz format or 
different mappers may be specified)
```
cd test_long_read_aligner_scattered
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/e093fd72-e31a-11ee-b020-27964ee37032--flagger_test_files/flagger_v0.4.0/test_files/test_long_read_aligner_scattered/test_1.tar.gz
```


