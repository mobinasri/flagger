version 1.0

workflow runCov2Wig{
    call cov2wig
    output{
        File wig = cov2wig.wig
    }
}

task cov2bigwig{
    input{
        File coverage
        Int windowLen=1000
        Int minHighMapqCov
        Int minMapq
        String trackName
        File hap1Fai
        File hap2Fai
        # runtime configurations
        Int memSize=32
        Int threadCount=8
        Int diskSize=ceil(size(coverage, "GB"))  + 64
        String dockerImage="mobinasri/flagger:v1.2.0"
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

        COV_FILE_PATH="~{coverage}"
        PREFIX=$(echo $(basename ${COV_FILE_PATH%.gz}) | sed -e 's/\.cov$//' -e 's/\.bed$//')

        cat ~{hap1Fai} ~{hap2Fai} > asm_dip.fa.fai

        mkdir output
        # create bedgraph
        coverage_format_converter -w ~{windowLen}  -i ~{coverage} -f asm_dip.fa.fai -n ~{trackName} -t ~{threadCount} -o output/${PREFIX}.bedgraph

        # convert bedgraph to bigwig with UCSC tool
        cat asm_dip.fa.fai | cut -f1-2 > asm.chrom.sizes

        for BEDGRAPH in $(ls output/*.bedgraph)
        do
            bedGraphToBigWig ${BEDGRAPH} asm.chrom.sizes ${BEDGRAPH%%.bedgraph}.bigwig 
        done

        # saving mappable regions in a BED file
        cat output/${PREFIX}.high_mapq.bedgraph | \
            grep -v "track type" | \
            awk '$4 >= ~{minHighMapqCov} {print $1"\t"$2"\t"$3}' | \
            bedtools sort -i - | \
            bedtools merge -i - > output/${PREFIX}.mapq_ge_~{minMapq}.cov_ge_~{minHighMapqCov}.mappable.bed


        cut -f1 ~{hap1Fai} | \
            grep -F -f - output/${PREFIX}.mapq_ge_~{minMapq}.cov_ge_~{minHighMapqCov}.mappable.bed > output/${PREFIX}.mapq_ge_~{minMapq}.cov_ge_~{minHighMapqCov}.mappable.hap1.bed
        cut -f1 ~{hap2Fai} | \
            grep -F -f - output/${PREFIX}.mapq_ge_~{minMapq}.cov_ge_~{minHighMapqCov}.mappable.bed > output/${PREFIX}.mapq_ge_~{minMapq}.cov_ge_~{minHighMapqCov}.mappable.hap2.bed


    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output{
        Array[File] bigwigArray = glob("output/*.bigwig")
        File mappableHap1Bed = glob("output/*.mappable.hap1.bed")[0]
        File mappableHap2Bed = glob("output/*.mappable.hap2.bed")[0]
    }
}



task cov2wig{
    input{
        File covGz
        Int segmentSize=1024
        Int threshold=250
        String trackName
        File fai
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=256
        String dockerImage="mobinasri/flagger:v0.4.0"
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
        
        FILENAME=`basename ~{covGz}`
        PREFIX="${FILENAME%.cov.gz}"

        gunzip -c ~{covGz} > ${PREFIX}.cov
        mkdir output
        cov2wig -i ${PREFIX}.cov -s ~{segmentSize} -t ~{threshold} -f ~{fai} -o output/${PREFIX}.wig -n ~{trackName}
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output{
        File wig = glob("output/*.wig")[0]
    }
}

