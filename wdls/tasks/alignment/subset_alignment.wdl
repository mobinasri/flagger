version 1.0 

workflow runSubsetAlignment{
    call subsetAlignment
    output{
        File outputBam = subsetAlignment.outputBam
        File outputBai = subsetAlignment.outputBai
    }
}
task subsetAlignment {
    input {
        File inputBam
        File inputBai
        String region
        String suffix
        # runtime configurations
        Int memSize=16
        Int threadCount=4
        Int diskSize=500
        String dockerImage="mobinasri/bio_base:v0.2"
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
        
        BAM_FILENAME=$(basename ~{inputBam})
        BAM_PREFIX=${BAM_FILENAME%.bam}
       
        ln ~{inputBam} ${BAM_PREFIX}.bam
        ln ~{inputBai} ${BAM_PREFIX}.bam.bai
        #samtools index ${BAM_PREFIX}.bam

        mkdir output
        samtools view -hb ${BAM_PREFIX}.bam ~{region} > output/${BAM_PREFIX}.~{suffix}.bam 
        java -jar /home/apps/jvarkit/dist/jvarkit.jar biostar84452 -o output/${BAM_PREFIX}.~{suffix}.no_clip_seq.bam --samoutputformat BAM output/${BAM_PREFIX}.~{suffix}.bam
        samtools index output/${BAM_PREFIX}.~{suffix}.no_clip_seq.bam
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones : zones
        cpuPlatform: "Intel Cascade Lake"
    }
    output {
        File outputBam = glob("output/*.no_clip_seq.bam")[0]
        File outputBai = glob("output/*.bai")[0]
    }
}

