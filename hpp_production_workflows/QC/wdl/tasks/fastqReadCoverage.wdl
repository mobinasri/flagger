version 1.0

workflow fastqReadCoverage {

    call sumFastqReads 

    output {
        File coverageFile = sumFastqReads.coverageFile
    }
}



task sumFastqReads {

    input {
        Array[File] inputFastq

        Int memSizeGB = 4
        Int diskSizeGB = 128

        ## Must give container with gawk (not mawk which uses 4-byte integers)
        String dockerImage = "registry.access.redhat.com/ubi7/ubi:latest"
    }

    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace

        BASES=0

        for fq in ~{sep=' ' inputFastq}
        do
              FILE_BASES=$(zcat "${fq}" | awk '{if (NR % 2 == 0 && NR % 4 != 0) s+=length($0)} END {printf "%d" ,s}')
              BASES=$(( $BASES + $FILE_BASES ))
        done

        echo $BASES > total_coverage.txt
    >>>

    output {

        File coverageFile  = "total_coverage.txt"
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}