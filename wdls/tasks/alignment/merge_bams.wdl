version 1.0

workflow MergeBamFiles{
    call merge
    output{
        File mergedBam = merge.mergedBam
        File mergedBai = merge.mergedBai
    }
}

task merge{
    input{
        Array[File] sortedBamFiles
        String sampleName
        String sampleSuffix
        # runtime configurations
        Int memSize=32
        Int threadCount=8
        Int diskSize=256
        String dockerImage="mobinasri/bio_base:latest"
        Int preemptible=2
        String zones
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
        
        samtools merge -@ ~{threadCount} ~{sampleName}.~{sampleSuffix}.bam ~{sep=" " sortedBamFiles}
        samtools index -@ ~{threadCount} ~{sampleName}.~{sampleSuffix}.bam
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones: zones
    }
    output {
        File mergedBam = "~{sampleName}.~{sampleSuffix}.bam"
        File mergedBai = "~{sampleName}.~{sampleSuffix}.bam.bai"
    }
}

