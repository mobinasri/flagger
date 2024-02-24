version 1.0 

workflow runPdfGenerator{
    call pdfGenerator
    output {
        File pdf = pdfGenerator.pdf
    }
}
task pdfGenerator{
    input {
        File windowProbTablesTarGz
        File genomeProbTable
        Boolean sortPdfPagesByHaplotype=false
        String hap1ContigPattern = "h1tg"
        String hap2ContigPattern = "h2tg"
        # runtime configurations
        Int memSize=16
        Int threadCount=8
        Int diskSize=128
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
        
        mkdir tables
        tar --strip-components 1 -xvzf ~{windowProbTablesTarGz} --directory tables
        FILENAME=$(basename ~{windowProbTablesTarGz})
        PREFIX=${FILENAME%.tables.tar.gz}
        python3 ${PDF_GENERATOR_PY} \
            --table  ~{genomeProbTable} \
            --dir tables \
            --pdf ${PREFIX}.cov_dist.pdf \
            --hap1Pattern "~{hap1ContigPattern}" \
            --hap2Pattern "~{hap2ContigPattern}" \
            ~{true="--diploid" false="" sortPdfPagesByHaplotype}
        
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File pdf = glob("*.pdf")[0]
    }
}

