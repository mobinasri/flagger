version 1.0

workflow runCov2Bin{
    call cov2bin
    output{
        File bin = cov2bin.bin
    }
}

task cov2bin{
    input{
        File coverage
        Int windowLen=4000
        Int chunkLen=20000000
        File fai
        # runtime configurations
        Int memSize=32
        Int threadCount=8
        Int diskSize=ceil(size(coverage, "GB"))  + 64
        String dockerImage="mobinasri/flagger:v1.1.0"
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

        FILENAME=$(basename ${COV_FILE_PATH})
        ln -s ${COV_FILE_PATH} ${FILENAME}

        mkdir output
        # create bedgraph
        coverage_format_converter -c ~{chunkLen} -w ~{windowLen} -i ${FILENAME} -f ~{fai} -t ~{threadCount} -o output/${PREFIX}.bin

    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File bin = glob("output/*.bin")[0]
    }
}


