version 1.0 

workflow runSamtoolsView{
    call view
    output{
        File outputBamFile = view.outputBamFile
        File outputBaiFile = view.outputBaiFile
    }
}
task view {
    input {
        File bamFile
        File baiFile
        String prefix
        String options
        String location
        # runtime configurations
        Int memSize=8
        Int threadCount=8
        Int diskSize=2 * ceil(size(bamFile, "GB")) + 64
        String dockerImage="mobinasri/bio_base:v0.1"
        Int preemptible=2
        String zones="us-west2-a"
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
        
        BAM_FILENAME=$(basename ~{bamFile})
        BAM_PREFIX=${BAM_FILENAME%.bam}

        ln ~{bamFile} ${BAM_PREFIX}.bam
        ln ~{baiFile} ${BAM_PREFIX}.bam.bai

        mkdir output
        samtools view -@8 -h -b ~{options} ${BAM_PREFIX}.bam ~{location} > output/${BAM_PREFIX}.~{prefix}.bam
        samtools index output/${BAM_PREFIX}.~{prefix}.bam
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones : zones
    }
    output {
        File outputBamFile = glob("output/*.bam")[0]
        File outputBaiFile = glob("output/*.bam.bai")[0]
    }
}

