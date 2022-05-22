version 1.0

workflow runBamCoverage{
    call bamCoverage
    output{
        File counts = bamCoverage.counts
        File coverageGz = bamCoverage.coverageGz
        Float coverageMeanFloat = bamCoverage.coverageMean
        Float coverageSdFloat = bamCoverage.coverageSD
    }
}

task bamCoverage{
    input{
        File bam
        Int minMAPQ = 0
        Int maxMAPQ = 60
        String extraOptions = ""
        File assemblyFastaGz
        # runtime configurations
        Int memSize=16
        Int threadCount=4
        Int diskSize=256
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
        
        # Extract assembly and index
        ASM_FILENAME=`basename ~{assemblyFastaGz}`
        ASM_PREFIX="${ASM_FILENAME%.*.*}"
        gunzip -c ~{assemblyFastaGz} > ${ASM_PREFIX}.fa
        samtools faidx ${ASM_PREFIX}.fa

        BAM_FILENAME=$(basename ~{bam})
        BAM_PREFIX=${BAM_FILENAME%.bam}

        if [ ~{maxMAPQ} -lt 60 ]; then
            correct_bam -x~{maxMAPQ} -a 0 -m 0 -p -n~{threadCount} -i ~{bam} -o ${BAM_PREFIX}.bam
        else
            ln ~{bam} ${BAM_PREFIX}.bam 
        fi
        
        # create be file to only output the regions in the given assembly
        cat ${ASM_PREFIX}.fa.fai | awk '{print $1"\t"0"\t"$2}' | bedtools sort -i - > ${ASM_PREFIX}.fa.fai.bed
        samtools depth -b ${ASM_PREFIX}.fa.fai.bed -aa -Q ~{minMAPQ} ~{extraOptions} ${BAM_PREFIX}.bam  > ${BAM_PREFIX}.depth

        # Convert the output of samtools depth into a compressed format
        depth2cov -d ${BAM_PREFIX}.depth -f ${ASM_PREFIX}.fa.fai -o ${BAM_PREFIX}.cov
        # Convert cov to counts
        cov2counts -i ${BAM_PREFIX}.cov -o ${BAM_PREFIX}.counts
        # Calculate mean and standard deviation
        python3 ${CALC_MEAN_SD_PY} --countsInput ${BAM_PREFIX}.counts --meanOutput cov_mean.txt --sdOutput cov_sd.txt
        gzip ${BAM_PREFIX}.cov
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output{
        File coverageGz = glob("*.cov.gz")[0]
        File counts = glob("*.counts")[0]
        Float coverageMean = read_float("cov_mean.txt") 
        Float coverageSD = read_float("cov_sd.txt")
    }
}

