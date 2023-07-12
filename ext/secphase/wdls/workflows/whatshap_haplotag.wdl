version 1.0

workflow runHaplotag{
    call haplotag
    output{
        File haplotaggedBam = haplotag.haplotaggedBam
        File haplotaggedBai = haplotag.haplotaggedBai
    }
}

task haplotag{
    input{
        File bam
        File bai
        File vcfGz
        File refFastaGz
        # runtime configurations
        Int memSize=8
        Int threadCount=4
        Int diskSize=512
        String dockerImage="quay.io/masri2019/hpp_coverage:latest"
        Int preemptible=2
    }
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
        
        FILENAME=`basename ~{bam}`
        PREFIX="${FILENAME%.*.*}"

        ln -s ~{bam} ${PREFIX}.bam
        ln -s ~{bai} ${PREFIX}.bai
        tabix ~{vcfGz}

        gunzip -c ~{refFastaGz} > ref.fa
        samtools faidx ref.fa

        mkdir output
        whatshap haplotag --reference ref.fa -o output/${PREFIX}.haplotagged.bam ~{vcfGz} ${PREFIX}.bam --ignore-read-groups --skip-missing-contigs
        samtools index output/${PREFIX}.haplotagged.bam 

    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output{
        File haplotaggedBam = glob("output/*.bam")[0]
        File haplotaggedBai = glob("output/*.bai")[0]
    }
}

