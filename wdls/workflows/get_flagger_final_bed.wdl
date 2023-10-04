version 1.0 

workflow runGetFinalBed{
    call getFinalBed
    output {
        File finalBed = getFinalBed.finalBed
        File simplifiedFinalBed = getFinalBed.simplifiedFinalBed
    }
}

task getFinalBed {
    input {
        File correctedBedsTarGz
        File altRemovedBedsTarGz
        String sampleName
        String suffix
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=32
        String dockerImage="mobinasri/flagger:v0.3.1"
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
       
        mkdir output
        bash /home/scripts/combine_alt_removed_beds.sh \
            -a ~{correctedBedsTarGz} \
            -b ~{altRemovedBedsTarGz} \
            -m /home/scripts/colors.txt \
            -t ~{sampleName}.~{suffix} \
            -o output/~{sampleName}.~{suffix}.flagger_final.bed
   
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File finalBed = glob("output/*.flagger_final.bed")[0]
        File simplifiedFinalBed = glob("output/*.flagger_final.simplified.bed")[0]
    }
}

