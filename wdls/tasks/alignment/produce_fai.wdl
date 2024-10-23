version 1.0 

workflow runProduceFai{
    call produceFai
    output {
        File fai = produceFai.fai
    }
}
task produceFai {
    input {
        File fasta
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=ceil(size(fasta, "GB")) + 32
        String dockerImage="mobinasri/bio_base:v0.4.0"
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

        FA_PREFIX=$(echo $(basename ~{fasta}) | sed -e 's/\.fa$//' -e 's/\.fa.gz$//' -e 's/\.fasta$//' -e 's/\.fasta.gz$//')

        if [[ "~{fasta}" == *.gz ]]; then
            gunzip -c ~{fasta} > ${FA_PREFIX}.fa
        else
            cp ~{fasta} ${FA_PREFIX}.fa
        fi
        samtools faidx ${FA_PREFIX}.fa

    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File fai = glob("*.fai")[0]
    }
}

