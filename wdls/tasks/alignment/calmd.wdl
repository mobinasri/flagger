version 1.0 

workflow runCalMD{
    call calmd
    output{
        File outputBamFile = calmd.outputBamFile
    }
}
task calmd {
    input {
        File bamFile
        File assemblyFastaGz
        # runtime configurations
        Int memSize=16
        Int threadCount=8
        Int diskSize=500
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
       
        ASM_FILENAME=$(basename ~{assemblyFastaGz})
        ASM_PREFIX=${ASM_FILENAME%.fa.gz}
        gunzip -c  ~{assemblyFastaGz} > ${ASM_PREFIX}.fa
        samtools calmd -@8 -b ~{bamFile} ${ASM_PREFIX}.fa > ${BAM_PREFIX}.bam
        samtools index ${BAM_PREFIX}.bam
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
        File outputBamFile = glob("*.bam")[0]
        File outputBaiFile = glob("*.bam.bai")[0]
    }
}

