version 1.0 

workflow runSubsetCoverage{
    call subsetCoverage
    output {
        File subsetCoverageGz = subsetCoverage.outputCoverageGz
    }
}


task subsetCoverage {
    input {
        File coverageGz
        File blocksBed
        String suffix
        # runtime configurations
        Int memSize=8
        Int threadCount=4
        Int diskSize=128
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

        FILENAME=$(basename ~{coverageGz})
        PREFIX=${FILENAME%.cov.gz}

        zcat ~{coverageGz} | \
            awk '{if(substr($1,1,1) == ">") {contig=substr($1,2,40); len_contig=$2} else {print contig"\t"$1-1"\t"$2"\t"$3"\t"len_contig}}' | \
            bedtools intersect -a - -b ~{blocksBed} | \
            awk '{if(contig != $1){contig=$1; print ">"contig" "$5}; print $2+1"\t"$3"\t"$4}' | pigz -p4 > ${PREFIX}.~{suffix}.cov.gz
    >>>
    runtime {
        docker: dockerImage
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File outputCoverageGz = glob("*.${suffix}.cov.gz")[0]
    }
}
