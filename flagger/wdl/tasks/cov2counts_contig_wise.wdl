version 1.0 

workflow runCov2CountsContigWise{
    call cov2countsContigWise
    output {
        File contigCountsTarGz = cov2countsContigWise.contigCountsTarGz
        File contigCovsTarGz = cov2countsContigWise.contigCovsTarGz
        File windowsText = cov2countsContigWise.windowsText
    }
}

task cov2countsContigWise {
    input {
        File coverageGz
        File fai
        Int windowSize
        # runtime configurations
        Int memSize=8
        Int threadCount=8
        Int diskSize=64
        String dockerImage="quay.io/masri2019/hpp_coverage:latest"
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

        FILENAME=$(basename ~{fai})
        PREFIX_FAI=${FILENAME%.fa.fai}

        FILENAME=$(basename ~{coverageGz})
        PREFIX_COV=${FILENAME%.cov.gz}
        
        gunzip -c ~{coverageGz} > ${PREFIX_COV}.cov
        mkdir covs counts
        # Make a separate cov file for each contig
        split_contigs_cov_v2 -c ${PREFIX_COV}.cov -f ~{fai} -p covs/${PREFIX_COV} -s ~{windowSize} > ${PREFIX_FAI}.windows.txt
        # Count each window-specific cov file
        for c in $(ls covs);do cov2counts -i covs/$c -o counts/${c/.cov/.counts}; echo $c" finished";done

        # Compress Counts files
        tar -cf ${PREFIX_COV}.counts.tar counts
        gzip ${PREFIX_COV}.counts.tar

        # Compress Cov files
        tar -cf ${PREFIX_COV}.covs.tar covs
        gzip ${PREFIX_COV}.covs.tar
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File contigCountsTarGz = glob("*.counts.tar.gz")[0]
        File contigCovsTarGz = glob("*.covs.tar.gz")[0]
        File windowsText = glob("*.windows.txt")[0]
    }
}

