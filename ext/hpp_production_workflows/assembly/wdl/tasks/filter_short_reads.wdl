version 1.0

import "../../../QC/wdl/tasks/extract_reads.wdl" as extractReads_t


workflow FilterShortReads {
    input {
        Array[File] readFiles
        Int minReadLength
        File? referenceFasta
    }

    scatter (readFile in readFiles){
        call extractReads_t.extractReads as extractReads {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=ceil(3 * size(readFile, "GB")) + 64,
                dockerImage="tpesout/hpp_base:latest"
        }
        call filterShortReads{
            input:
                readFastq = extractReads.extractedRead,
                diskSizeGB = ceil(3 * size(extractReads.extractedRead, "GB")) + 64,
                minReadLength = minReadLength
        }
    }

    output {
        Array[File] longReadFastqs = filterShortReads.longReadFastq 
    }
}


# This task does not work properly if sequences have methylation tags
task filterShortReads {
    input{
        File readFastq
        Int minReadLength
        # runtime configurations
        Int memSizeGB=8
        Int threadCount=4
        Int diskSizeGB=512
        Int preemptible=1
        String dockerImage="tpesout/hpp_base:latest"
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        mkdir data
        cd data
        FILENAME=$(basename -- "~{readFastq}")
        PREFIX="${FILENAME%.*}"
        # filter reads shorter than minReadLength
        awk 'NR%4==1{a=$0} NR%4==2{b=$0} NR%4==3{c=$0} NR%4==0&&length(b)>~{minReadLength}{print a"\n"b"\n"c"\n"$0;}' ~{readFastq} > ${PREFIX}.long.fastq
        OUTPUTSIZE=`du -s -BG *.long.fastq | sed 's/G.*//'`
        echo $OUTPUTSIZE > outputsize
    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible: preemptible
    }

    output {
        File longReadFastq = glob("data/*.long.fastq")[0]
        Int fileSizeGB = read_int("data/outputsize")
    }
}

