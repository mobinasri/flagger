version 1.1

workflow runBamCoverage{
    call bamCoverageFast
    output{
        File counts = bamCoverageFast.counts
        File coverageGz = bamCoverageFast.coverageGz
        Float coverageModeFloat = bamCoverageFast.coverageMode
        Float coverageSdFloat = bamCoverageFast.coverageSD
    }
}

task bamCoverageFast{
    input{
        File bam
        File bai
        Int minMAPQ
        String output_format = "only_total" # it can be "only_high_mapq"
        File assemblyFastaGz
        # runtime configurations
        Int memSize=32
        Int threadCount=16
        Int diskSize=ceil(size(bam, "GB"))  + 512
        String dockerImage="mobinasri/flagger:v0.4.0"
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
        FA_FILENAME=`basename ~{assemblyFastaGz}`
        FA_PREFIX="${FA_FILENAME%.*.*}"
        gunzip -c ~{assemblyFastaGz} > ${FA_PREFIX}.fa
        samtools faidx ${FA_PREFIX}.fa

        cat ${FA_PREFIX}.fa.fai | awk '{print $1"\t0\t"$2}' | bedtools sort -i - > asm_wg.bed

        # make a json file pointing to the asm_wg.bed
        echo "{" > bed_file.json
        echo \"asm_wg\" : \"asm_wg.bed\" >> bed_file.json
        echo "}" >> bed_file.json

        BAM_FILENAME=`basename ~{bam}`
        BAM_PREFIX="${BAM_FILENAME%.*}"
        ln -s ~{bam} ${BAM_PREFIX}.bam
        ln -s ~{bai} ${BAM_PREFIX}.bam.bai

      
        COV_PREFIX=${BAM_PREFIX}.mapq~{minMAPQ}.~{output_format}
        # convert bam to cov
        # -f only_total is for printing the total coverage only
        bam2cov -i ${BAM_PREFIX}.bam -O "c" -m ~{minMAPQ} -o ${COV_PREFIX}.cov -j bed_file.json -t~{threadCount} -f ~{output_format}

        # Convert cov to counts
        cov2counts -i ${COV_PREFIX}.cov -o ${COV_PREFIX}.counts
        # Calculate mod and standard deviation
        python3 ${CALC_MODE_SD_PY} --countsInput ${COV_PREFIX}.counts --minCoverage 5 --modeOutput cov_mode.txt --sdOutput cov_sd.txt
        gzip ${COV_PREFIX}.cov
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
        Float coverageMode = read_float("cov_mode.txt")
        Float coverageSD = read_float("cov_sd.txt")
    }
}



task bamCoverage{
    input{
        File bam
        Int minMAPQ
        File assemblyFastaGz
        # runtime configurations
        Int memSize=16
        Int threadCount=4
        Int diskSize=ceil(size(bam, "GB"))  + 512
        String dockerImage="mobinasri/flagger:v0.4.0"
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
        FA_FILENAME=`basename ~{assemblyFastaGz}`
        FA_PREFIX="${FA_FILENAME%.*.*}"
        gunzip -c ~{assemblyFastaGz} > ${FA_PREFIX}.fa
        samtools faidx ${FA_PREFIX}.fa

        BAM_FILENAME=`basename ~{bam}`
        BAM_PREFIX="${BAM_FILENAME%.*}"
        samtools depth -aa -Q ~{minMAPQ} ~{bam}  > ${BAM_PREFIX}.depth
        
        COV_PREFIX=${BAM_PREFIX}.mapq~{minMAPQ}
        # Convert the output of samtools depth into a compressed format
        depth2cov -d ${BAM_PREFIX}.depth -f ${FA_PREFIX}.fa.fai -o ${COV_PREFIX}.cov
        # Convert cov to counts
        cov2counts -i ${COV_PREFIX}.cov -o ${COV_PREFIX}.counts
        # Calculate mean and standard deviation
        python3 ${CALC_MODE_SD_PY} --countsInput ${COV_PREFIX}.counts --minCoverage 5 --modeOutput cov_mode.txt --sdOutput cov_sd.txt
        gzip ${COV_PREFIX}.cov
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
        Float coverageMode = read_float("cov_mode.txt") 
        Float coverageSD = read_float("cov_sd.txt")
    }
}

