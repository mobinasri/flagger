version 1.0

workflow runFindBlocks{
    call findBlocks
    output {
        File bedsTarGz = findBlocks.bedsTarGz
    }
}


task findBlocks {
    input {
        File coverageGz
        File table
        String suffix="whole_genome"
        # runtime configurations
        Int memSize=8
        Int threadCount=4
        Int diskSize=128
        String dockerImage="mobinasri/flagger:dev-v0.1"
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

        FILENAME=$(basename ~{coverageGz})
        PREFIX=${FILENAME%.cov.gz}

        gunzip -c ~{coverageGz} > $PREFIX.cov
        if [ -z ~{suffix} ]; then
            mkdir blocks
            find_blocks_from_table -c $PREFIX.cov -t ~{table} -p blocks/${PREFIX}
            for bed in $(ls ~{suffix})
            do
                sort -k1,1 -k2,2n blocks/${bed} > tmp.bed
                mv tmp.bed blocks/${bed}
                rm -rf tmp.bed
            done
            tar -cf ${PREFIX}.beds.tar blocks
            gzip ${PREFIX}.beds.tar
        else
            mkdir ~{suffix}
            find_blocks_from_table -c $PREFIX.cov -t ~{table} -p ~{suffix}/${PREFIX}.~{suffix}
            for bed in $(ls ~{suffix}) 
            do
                sort -k1,1 -k2,2n ~{suffix}/${bed} > tmp.bed
                mv tmp.bed ~{suffix}/${bed}
                rm -rf tmp.bed 
            done
            tar -cf ${PREFIX}.beds.~{suffix}.tar ~{suffix}
            gzip ${PREFIX}.beds.~{suffix}.tar
        fi
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File bedsTarGz = glob("*.tar.gz")[0]
    }
}
