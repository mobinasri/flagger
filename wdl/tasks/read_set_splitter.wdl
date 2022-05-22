version 1.0 

workflow runReadSetSplitter{
    call readSetSplitter
}

task readSetSplitter {
    input {
        Array[File] readFastqs
        Int splitNumber
        # runtime configurations
        Int memSize=16
        Int threadCount=8
        Int diskSize=1024
        String dockerImage="tpesout/hpp_base:latest"
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
        
        readsPerFile=$(cat ~{sep=" " readFastqs} | wc -l | awk '{print int($1/(4 * ~{splitNumber})) + 1}')
        linesPerFile=$(( readsPerFile * 4 ))     
        cat ~{sep=" " readFastqs} | split -l $linesPerFile  --additional-suffix ".split.fq"
        for i in *.split.fq;do du -s -BG $i | sed 's/G.*//';done > outputsize.txt
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        Array[File] splitReadFastqs = glob("*.split.fq")
        Array[Float] readSizes = read_lines("outputsize.txt")
    }
}

