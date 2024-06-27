version 1.0

workflow chopperFilterWF {
    input {
        Array[File] readFiles
        Int minReadLength = 10000
        Int minReadQual   = 20
    }

    scatter (readFile in readFiles){
        call chopperFilter {
            input:
                readFastqGZ=readFile,
                minReadLength=minReadLength,
                minReadQual=minReadQual,
                memSizeGB=24,
                threadCount=6,
                diskSizeGB=ceil(3 * size(readFile, "GB")) + 64
        }
    }

    output {
        Array[File] filteredReads = chopperFilter.filteredFastqGZ 
    }
}


# This task does not work properly if sequences have methylation tags
task chopperFilter {
    input{
        File readFastqGZ
        Int minReadLength
        Int minReadQual

        # runtime configurations
        Int memSizeGB=24
        Int threadCount=6
        Int diskSizeGB=512
        Int preemptible=1
        String dockerImage="quay.io/biocontainers/chopper:0.2.0--hd03093a_0"
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        mkdir filtered

        FILENAME=$(basename -- "~{readFastqGZ}")
        PREFIX="${FILENAME%.*}"


        gunzip -c ~{readFastqGZ} \
            | chopper \
                -q ~{minReadQual} \
                -l ~{minReadLength} \
                | gzip \
                > filtered/${PREFIX}_gt~{minReadLength}_gt~{minReadQual}.fastq.gz
    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible: preemptible
    }

    output {
        File filteredFastqGZ = glob("filtered/*.fastq.gz")[0]
    }
}

