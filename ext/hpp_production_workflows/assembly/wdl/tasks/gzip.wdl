version 1.0 

workflow runGzip{
    call gzip
}
task gzip {
    input {
        File fileInput
        # runtime configurations
        Int memSize=32
        Int threadCount=8
        Int diskSize=256
        String dockerImage="quay.io/masri2019/hpp_hifiasm:latest"
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
        
        FILENAME=$(basename "~{fileInput}")
        ln ~{fileInput} ${FILENAME}
        pigz -p~{threadCount} -c ${FILENAME} > ${FILENAME}.gz
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File fileGz = glob("*.gz")[0]
    }
}

