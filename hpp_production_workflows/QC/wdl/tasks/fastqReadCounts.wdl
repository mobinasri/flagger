version 1.0

workflow fastqReadCounts {

    call countFastqReads 

    output {
        File totalReadsFile = countFastqReads.totalReadsFile
    }
}



task countFastqReads {

    input {
        Array[File] inputFastq

        Int memSizeGB = 4
        Int diskSizeGB = 128
        String dockerImage = "biocontainers/samtools:v1.9-4-deb_cv1"
    }

    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace

        READ_COUNT=0

        for fq in ~{sep=' ' inputFastq}
        do
              FILE_COUNT=$(zcat "${fq}" | wc -l )/4
              READ_COUNT=$(( $READ_COUNT + $FILE_COUNT ))
        done

        echo $READ_COUNT > total_reads.txt
    >>>

    output {

        File totalReadsFile  = "total_reads.txt"
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}