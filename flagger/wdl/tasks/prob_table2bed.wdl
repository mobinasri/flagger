version 1.0 

workflow runProbTable2Bed{
    call table2bed
    output{
        File unreliableBed = table2bed.unreliableBed
    }
}

task table2bed {
    input {
        File hifiTable
        File ontTable
        File hifiCovGz
        File ontCovGz
        Float threshold=0.5
        String sampleName
        String sampleSuffix
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=64
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

        gunzip -c ~{hifiCovGz} > ~{sampleName}.~{sampleSuffix}.hifi.cov
        gunzip -c ~{ontCovGz} > ~{sampleName}.~{sampleSuffix}.ont.cov
        find_blocks_from_tables -t ~{hifiTable} -T ~{ontTable} -c ~{sampleName}.~{sampleSuffix}.hifi.cov -C ~{sampleName}.~{sampleSuffix}.ont.cov  -d ~{threshold} -o ~{sampleName}.~{sampleSuffix}.model.unreliable.unmerged.bed
        bedtools merge -i ~{sampleName}.~{sampleSuffix}.model.unreliable.unmerged.bed > ~{sampleName}.~{sampleSuffix}.model.unreliable.bed
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File unreliableBed = "~{sampleName}.~{sampleSuffix}.model.unreliable.bed"
    }
}

