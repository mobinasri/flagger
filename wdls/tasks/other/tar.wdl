version 1.1 

workflow runTarGz{
    call tarGz
}
task tarGz {
    input {
        String tarGzName
        Array[File] files
        # runtime configurations
        Int memSize=16
        Int threadCount=8
        Int diskSize=128
        String dockerImage="mobinasri/bio_base:v0.4.0"
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
        
        mkdir ~{tarGzName}
        cp ~{sep=" " files} ~{tarGzName}
        tar -cf ~{tarGzName}.tar ~{tarGzName}
        gzip ~{tarGzName}.tar
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File fileTarGz = "~{tarGzName}.tar.gz"
    }
}

