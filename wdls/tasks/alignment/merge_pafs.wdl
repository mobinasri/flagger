version 1.0

workflow MergePafFiles{
    call merge
    output{
        File mergedPaf = merge.mergedPaf
    }
}

task merge{
    input{
        Array[File] pafFiles
        String sampleName
        String sampleSuffix
        # runtime configurations
        Int memSize=8
        Int threadCount=2
        Int diskSize=256
        String dockerImage="mobinasri/bio_base:v0.1"
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
        
        cat ~{sep=" " pafFiles} | sort -k6,8 -n > ~{sampleName}.~{sampleSuffix}.paf
    >>>
    runtime {
        docker: dockerImage
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones: zones
    }
    output {
        File mergedPaf = "~{sampleName}.~{sampleSuffix}.paf"
    }
}

