version 1.0

workflow runGenerateWig{
    call generateWig
    output{
        File wig = generateWig.wig
    }
}

task generateWig{
    input{
        String sampleName
        String sampleSuffix
        String platform
        File bamFile
        Int minMAPQ
        String extraOptions = "-z 5 -w 25 -e 250"
        File assemblyFastaGz
        # runtime configurations
        Int memSize=16
        Int threadCount=4
        Int diskSize=512
        String dockerImage="mobinasri/flagger:v0.2"
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
        
        # Extract assembly and index
        FILENAME=`basename ~{assemblyFastaGz}`
        PREFIX="${FILENAME%.*.*}"
        gunzip -c ~{assemblyFastaGz} > $PREFIX.fa
        samtools faidx $PREFIX.fa

        mkdir output
        samtools view -h -F256 -F4 -q ~{minMAPQ} ~{bamFile} > output/pri.bam
        igvtools count ~{extraOptions} output/pri.bam output/~{sampleName}.~{sampleSuffix}.~{platform}.wig $PREFIX.fa
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

