version 1.0

workflow runWig2Tdf{
    call wig2tdf
    output{
        File tdf = wig2tdf.tdf
    }
}

task wig2tdf{
    input{
        File wig
        Int z=7
        File assemblyFastaGz
        # runtime configurations
        Int memSize=8
        Int threadCount=4
        Int diskSize=256
        String dockerImage="mobinasri/flagger:v0.3.0"
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
        
        FILENAME=`basename ~{wig}`
        PREFIX="${FILENAME%.wig}"

        gunzip -c ~{assemblyFastaGz} > asm.fa
        mkdir output
        igvtools toTDF -z ~{z} ~{wig} output/${PREFIX}.tdf asm.fa
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output{
        File tdf = glob("output/*.tdf")[0]
    }
}

