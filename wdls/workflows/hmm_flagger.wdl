version 1.0 

workflow runHmmFlagger{
    call hmmFlagger
    output {
        File bed = hmmFlagger.bed
    }
}

task hmmFlagger {
    input {
        File mergedCovGz
        Float coverage
        Int regions
        String regionFactors
        Int chunkLen = 20000000
        Int iterations = 5
        Int windowLen = 1000
        Float minColScore = 1.8
        Int minColLen = 40000
        Float maxDupScore = 0.4
        Int minDupLen = 40000
        Float maxHighMapqRatio = 0.2
        # runtime configurations
        Int memSize=32
        Int threadCount=8
        Int diskSize=32
        String dockerImage="quay.io/masri2019/hpp_coverage:latest"
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
        FILENAME=$(basename ~{mergedCovGz})
        PREFIX=${FILENAME%%.cov.gz}

        gunzip -c ~{mergedCovGz} > ${PREFIX}.cov

        # index coverage file
        index_cov -c ${PREFIX}.cov -l ~{chunkLen} 

        # run hmm_flagger
        hmm_flagger --inputCov ${PREFIX}.cov \
                    --threads ~{threadCount} \
                    --chunkLen ~{chunkLen} \
                    --iterations ~{iterations} \
                    --windowLen ~{windowLen} \
                    --coverage ~{coverage} \
                    --regions  ~{regions} \
                    --trackName ${PREFIX} \
                    --outputDir output \
                    --minColScore ~{minColScore} \
                    --minColLen ~{minColLen} \
                    --maxDupScore ~{maxDupScore} \
                    --minDupLen ~{minDupLen} \
                    --regionFactors ~{regionFactors} \
                    --maxHighMapqRatio ~{maxHighMapqRatio}

    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File bed = glob("output/*.bed")[0]
    }
}

