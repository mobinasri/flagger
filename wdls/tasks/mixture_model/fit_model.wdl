version 1.0 

workflow runFitModel{
    call fitModel
    output {
       File probabilityTable = fitModel.probabilityTable 
    }
}
task fitModel {
    input {
        File counts
        Float cov=20.0
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=32
        String dockerImage="mobinasri/flagger:v0.3.2"
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

        FILENAME=$(basename ~{counts})
        PREFIX=${FILENAME%.counts}
        python3 ${FIT_GMM_PY} --counts ~{counts} --output ${PREFIX}.table --cov ~{cov}
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File probabilityTable = glob("*.table")[0]
    }
}

