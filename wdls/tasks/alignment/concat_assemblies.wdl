version 1.0 

workflow runConcatAssemblies{
    call concatAssemblies
    output {
        File diploidAssembly = concatAssemblies.diploidAssembly
    }
}

task concatAssemblies {
    input {
        File assemblyHap1FastaGz
        File assemblyHap2FastaGz
        String sampleName
        String sampleSuffix
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=32
        String dockerImage="mobinasri/bio_base:v0.1"
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
        
        zcat ~{assemblyHap1FastaGz} ~{assemblyHap2FastaGz} | pigz -p2 >  ~{sampleName}.~{sampleSuffix}.dip.fa.gz
        
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File diploidAssembly = glob("*.dip.fa.gz")[0]
    }
}
