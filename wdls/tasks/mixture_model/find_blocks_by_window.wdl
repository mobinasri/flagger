version 1.0

workflow runFindBlocksByWindow{
    call findBlocksByWindow
    output {
        File windowBedsTarGz = findBlocksByWindow.windowBedsTarGz
    }
}


task findBlocksByWindow {
    input {
        File windowCovsTarGz
        File windowProbTablesTarGz
        Int minWindowSize=5000000
        File windowsText
        String suffix="window_based"
        # runtime configurations
        Int memSize=8
        Int threadCount=8
        Int diskSize=64
        String dockerImage="mobinasri/flagger:v0.1"
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

        mkdir covs
        tar --strip-components 1 -xvzf ~{windowCovsTarGz} --directory covs
        FILENAME=$(basename ~{windowCovsTarGz})
        export PREFIX=${FILENAME%.covs.tar.gz}
        
        mkdir tables
        tar --strip-components 1 -xvzf ~{windowProbTablesTarGz} --directory tables
        
        mkdir tmp
        cat ~{windowsText} | awk '~{minWindowSize} <= $3' | xargs -n3 -P ~{threadCount} sh -c 'find_blocks_from_table -c covs/"${PREFIX}".$0_$1_$2.cov -t tables/"${PREFIX}".$0_$1_$2.table -p tmp/"${PREFIX}".$0_$1_$2'
        mkdir ~{suffix}
        cat tmp/*.error.bed | bedtools sort -i - | bedtools merge -i - > ~{suffix}/${PREFIX}.~{suffix}.error.bed
        cat tmp/*.duplicated.bed | bedtools sort -i - | bedtools merge -i - > ~{suffix}/${PREFIX}.~{suffix}.duplicated.bed
        cat tmp/*.haploid.bed | bedtools sort -i - | bedtools merge -i - > ~{suffix}/${PREFIX}.~{suffix}.haploid.bed
        cat tmp/*.collapsed.bed | bedtools sort -i - | bedtools merge -i - > ~{suffix}/${PREFIX}.~{suffix}.collapsed.bed
        tar -cf ${PREFIX}.beds.~{suffix}.tar ~{suffix}
        gzip ${PREFIX}.beds.~{suffix}.tar
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File windowBedsTarGz = glob("*.${suffix}.tar.gz")[0]
    }
}
