version 1.0

workflow runPhasedGFAs2Fasta{
    call phasedGFAs2Fasta
    output {
        File paternalFastaGz = phasedGFAs2Fasta.outputPaternalFastaGz
        File maternalFastaGz = phasedGFAs2Fasta.outputMaternalFastaGz
    }
}


task phasedGFAs2Fasta {
    input{
        File paternalGfa
        File maternalGfa
        String patSuffix = "pat"
        String matSuffix = "mat"
        String childID
        # runtime configurations
        Int memSizeGB=32
        Int threadCount=8
        Int diskSizeGB=128
        Int preemptible=1
        String dockerImage="quay.io/masri2019/hpp_hifiasm:latest"
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

        # Convert contig GFA to FASTA with https://github.com/lh3/gfatools (this can be done with awk, too). Run in parallel:
        gfatools gfa2fa ~{paternalGfa} | pigz -p~{threadCount} > ~{childID}.~{patSuffix}.fa.gz &
        gfatools gfa2fa ~{maternalGfa} | pigz -p~{threadCount} > ~{childID}.~{matSuffix}.fa.gz &
        wait
    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible: preemptible
    }

    output {
        File outputPaternalFastaGz = "~{childID}.~{patSuffix}.fa.gz"
        File outputMaternalFastaGz = "~{childID}.~{matSuffix}.fa.gz"
    }
}

