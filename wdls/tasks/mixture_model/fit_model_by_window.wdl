version 1.0 

workflow runFitModelByWindow{
    call fitModelByWindow
    output {
       File windowProbTablesTarGz = fitModelByWindow.windowProbTablesTarGz
    }
}

task fitModelByWindow {
    input {
        File windowsText
        File countsTarGz
        Float cov=20
        Int minContigSize=5000000
        # runtime configurations
        Int memSize=16
        Int threadCount=16
        Int diskSize=64
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

        mkdir counts tables
        tar --strip-components 1 -xvzf ~{countsTarGz} --directory counts
        
        
        FILENAME="$(basename ~{countsTarGz})"
        export PREFIX=${FILENAME%.counts.tar.gz}
        
        cat ~{windowsText} | awk '~{minContigSize} <= $3' | xargs -n 3 -P ~{threadCount} sh -c 'python3 ${FIT_GMM_PY}  --cov ~{cov} --counts counts/"${PREFIX}".$0_$1_$2.counts --output tables/"${PREFIX}".$0_$1_$2.table'
        tar -cf ${PREFIX}.tables.tar tables
        gzip ${PREFIX}.tables.tar
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File windowProbTablesTarGz = glob("*.tables.tar.gz")[0]
    }
}

