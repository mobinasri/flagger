version 1.0

workflow scatterGuppyGPU {
    input {
        Array[File] fast5_tar_files
        String sample_name
        String guppy_version
        Int desired_size_GB
        
    }

    parameter_meta {
        fast5_tar_files: "Input fast5 tar, files must be in .tar format"
        sample_name: "Name of sample, used for output file names"
        guppy_version: "Guppy version, used for output file names"
        desired_size_GB: "Choose size to split input tar file by. With a 300GB fast5_tar_file and 30GB desired_size_GB, the fast5_tar_file will be split in 10 pieces."
    }


    # scatter in case multiple tar files are given
    scatter (fast5_tar in fast5_tar_files) {
        call splitFast5s {
            input:
                file_to_split_tar = fast5_tar,
                desired_size_GB = desired_size_GB

        }

        # call guppyGPU on each of the smaller "split" tar files
        scatter (split_fast5_tar in splitFast5s.split_fast5_tar) {
            call guppyGPU {
                input:
                    fast5_tar_file = split_fast5_tar
            }
        }

        call concatenateBam {
            input:
                files = flatten(guppyGPU.pass_bam),
                sample_name = sample_name,
                guppy_version = guppy_version
        }

        call concatenateFastq {
            input:
                files = flatten(guppyGPU.pass_fastq),
                sample_name = sample_name,
                guppy_version = guppy_version
        }

        call concatenateSummary {
            input:
                files = guppyGPU.summary,
                sample_name = sample_name,
                guppy_version = guppy_version
        }

    }

    output {
        Array[File] bams = concatenateBam.concatenatedBam
        Array[File] fastqs = concatenateFastq.concatenatedFastq
        Array[File] summaries = concatenateSummary.concatenatedSummary
    }
    
    meta {
        author: "Jimin Park"
        email: "jpark621@ucsc.edu"
        description: "Calls guppy_basecaller with GPUs. Takes in fast5 tar file and outputs unaligned bam with methylation calls, fastq and sequencing summary text file."
    }
}


task splitFast5s {
    input {
        File file_to_split_tar
        Int desired_size_GB

        String dockerImage = "jiminpark/guppy-wdl:latest" 

        Int preempts = 3
        Int memSizeGB = 8
        Int extraDisk = 5
        Int threadCount = 2
    }

    Int file_size = ceil(size(file_to_split_tar, "GB"))
    Int diskSizeGB = 3 * file_size + extraDisk

    command <<<
        ## Extract tar file to 
        mkdir tmp
        
        # place all extracted files into directory tmp
        tar xvf "~{file_to_split_tar}" --directory tmp
        rm ~{file_to_split_tar}

        # move files into folder until exceeds desired_size_GB
        # then tar contents of folder
        OUTPUT_IDX=0
        OUTPUT_DIR=fast5_tar_$OUTPUT_IDX
        mkdir $OUTPUT_DIR
        for FILE in `find tmp/ -name "*.fast5"`
        do
            size=$(du -s -BG $OUTPUT_DIR | sed 's/G.*//')
            if (( $size > ~{desired_size_GB} ))
            then
                tar -cvf fast5_tarball_$OUTPUT_IDX.tar $OUTPUT_DIR/*
                rm -r $OUTPUT_DIR
                OUTPUT_IDX=$(($OUTPUT_IDX + 1))
                OUTPUT_DIR=fast5_tar_$OUTPUT_IDX
                mkdir $OUTPUT_DIR
            fi
            mv $FILE $OUTPUT_DIR
        done
        tar -cvf fast5_tarball_$OUTPUT_IDX.tar $OUTPUT_DIR/*
        rm -r $OUTPUT_DIR

    >>>

    output {
        Array[File] split_fast5_tar = glob("*tar")
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible : preempts
    }
}


task guppyGPU {
    
    input {
    
        File fast5_tar_file

        String CONFIG_FILE = "dna_r9.4.1_450bps_modbases_5mc_cg_sup_prom.cfg"
        Int READ_BATCH_SIZE = 250000
        Int q = 250000

        String dockerImage = "jiminpark/guppy-wdl:latest" 

        String? additionalArgs

        Int preempts = 3
        Int memSizeGB = 64
        Int threadCount = 12
        Int extraDisk = 5
        Int gpuCount = 1
        Int maxRetries = 4
        String gpuType = "nvidia-tesla-v100"
        String nvidiaDriverVersion = "418.87.00"
        String zones = "us-west1-b"
    }

    # calculate needed disk size
    Int file_size = ceil(size(fast5_tar_file, "GB"))
    Int diskSizeGB = 3 * file_size + extraDisk


    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        ## Extract tar file to 
        mkdir input
        
        # place all extracted files into directory input
        tar xvf "~{fast5_tar_file}" --directory input

        mkdir output

        # check if length of "additionalArgs" is zero

        if [[ "~{additionalArgs}" == "" ]]
        then
            ADDITIONAL_ARGS=""
        else
            ADDITIONAL_ARGS="~{additionalArgs}"
        fi


        guppy_basecaller \
            -i input/ \
            -s output/ \
            -c /opt/ont/guppy/data/"~{CONFIG_FILE}" \
            --bam_out \
            -x cuda:all:100% \
            -r \
            --read_batch_size "~{READ_BATCH_SIZE}" \
            -q "~{q}" \
            ${ADDITIONAL_ARGS}

    >>>

    output {
        Array[File] pass_bam = glob("output/pass/*.bam")
        Array[File] pass_fastq = glob("output/pass/*.fastq")
        File summary = glob("output/sequencing_summary.txt")[0]
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        gpuCount: gpuCount
        gpuType: gpuType
        maxRetries : maxRetries
        nvidiaDriverVersion: nvidiaDriverVersion
        preemptible : preempts
        docker: dockerImage
        zones: zones
    }
}


task concatenateBam {
    input {
        Array[File] files
        
        String sample_name
        String guppy_version

        String dockerImage = "tpesout/megalodon:latest"

        Int preempts = 3
        Int memSizeGB = 8
        Int threadCount = 3
        Int diskSizeGB = 500
    }
    
    command {
        samtools merge -o "${sample_name}_${guppy_version}.bam" ${sep=" " files}
    }

    output {
        File concatenatedBam = "${sample_name}_${guppy_version}.bam"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible : preempts
    }
}


task concatenateFastq {
    input {
        Array[File] files
        String sample_name
        String guppy_version

        String dockerImage = "tpesout/megalodon:latest"

        # runtime
        Int preempts = 3
        Int memSizeGB = 8
        Int threadCount = 3
        Int diskSizeGB = 500
    }
    
    command {
        cat ${sep=" " files} | gzip -c > "${sample_name}_${guppy_version}.fastq.gz"
    }

    output {
        File concatenatedFastq = "${sample_name}_${guppy_version}.fastq.gz"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible : preempts
    }
}


task concatenateSummary {
    input {
        Array[File] files
        String sample_name
        String guppy_version

        String dockerImage = "tpesout/megalodon:latest"

        # runtime
        Int preempts = 3
        Int memSizeGB = 8
        Int threadCount = 3
        Int diskSizeGB = 50
    }
    
    command {
        cat ${sep=" " files} > "tmp.txt"
        # remove duplicate headers
        awk 'NR==1 || !/^filename/' "tmp.txt" > "${sample_name}_${guppy_version}_sequencing_summary.txt"
    }

    output {
        File concatenatedSummary = "${sample_name}_${guppy_version}_sequencing_summary.txt"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible : preempts
    }
}
