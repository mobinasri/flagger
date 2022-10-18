version 1.0

workflow runPepperMarginDeepVariant{
    input {
        File assemblyFastaGz
        File bam
        File bamIndex
        ## Model types can be:
        ##     "ont_r9_guppy5_sup"
        ##     "ont_r10_q20"
        ##     "hifi"
        String pmdvModelType
        Int minMAPQ = 0
        String includeSupplementary="False"
        Boolean flagRemoveMultiplePrimary = true
    }
    if (flagRemoveMultiplePrimary) { 
        call removeMultiplePrimary{
            input:
                bam = bam,
                diskSize= 2 * ceil(size(bam, "GB")) + 64
        }
    }

    File bamForCalling =  select_first([removeMultiplePrimary.correctedBam, bam])
    File baiForCalling =  select_first([removeMultiplePrimary.correctedBai, bamIndex])
 
    call pmdv{
        input:
            modelType = pmdvModelType,
            assemblyFastaGz = assemblyFastaGz,
            bam = bamForCalling,
            bamIndex = baiForCalling,
            includeSupplementary = includeSupplementary,
            minMAPQ = minMAPQ,
            diskSize= 4 * ceil(size(bam, "GB")) + 256
    }
    output{
        File vcfGz = pmdv.vcfGz
        File intermediateTar = pmdv.intermediateTar
    }
}

task removeMultiplePrimary{
    input{
        File bam
        # runtime configurations
        Int memSize=8
        Int threadCount=8
        Int diskSize=512
        String dockerImage="mobinasri/flagger:dev-0.1"
        Int preemptible=3
        String zones="us-west2-a"
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

        BAM_NAME=$(basename ~{bam})
        BAM_PREFIX=${BAM_NAME%%.bam}

        mkdir output
        # Remove reads with multiple primary alignment (it caused an issue for margin before)
        samtools view -F 0x904 ~{bam} | cut -f 1 | sort | uniq -c | awk '$1 > 1' | cut -f2 > read_ids_multiple_primary.txt
        correct_bam -p -i ~{bam} -m 0 -a 0 -e read_ids_multiple_primary.txt -o output/${BAM_PREFIX}.bam -n ~{threadCount}
        samtools index output/${BAM_PREFIX}.bam
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones: zones
    }
    output {
        File correctedBam = glob("output/*.bam")[0]
        File correctedBai = glob("output/*.bai")[0]
    }

}

task pmdv{
    input{
        File bam
        File? bamIndex
        File assemblyFastaGz
        Int minMAPQ
        String includeSupplementary="False"
        String modelType
        # runtime configurations
        Int memSize=256
        Int threadCount=64
        Int diskSize=512
        String dockerImage="kishwars/pepper_deepvariant:r0.7"
        Int preemptible=3
        String zones="us-west2-a"
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
        
        BAM_NAME=$(basename ~{bam})
        BAM_PREFIX=${BAM_NAME%%.bam}
        ln ~{bam} ${BAM_PREFIX}.bam

        if [ -n "~{bamIndex}" ]
        then
            ln ~{bamIndex} ${BAM_PREFIX}.bam.bai
        else
            samtools index ${BAM_PREFIX}.bam
        fi

        ## unzip the fasta file and produce its index
        gunzip -c ~{assemblyFastaGz} > asm.fa
        samtools faidx asm.fa

        MORE_OPTIONS="--dv_min_mapping_quality ~{minMAPQ} --pepper_min_mapping_quality ~{minMAPQ}"
        if [ ~{includeSupplementary} == "True" ]; then
            MORE_OPTIONS="${MORE_OPTIONS} --pepper_include_supplementary"
        fi
       
        mkdir output

        ## call pepper-margin-deepvariant 
        run_pepper_margin_deepvariant call_variant \
        -b  ${BAM_PREFIX}.bam \
        -f  asm.fa \
        -o  output \
        -t ~{threadCount} \
        --~{modelType} ${MORE_OPTIONS}

        mv output/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz output/${BAM_PREFIX}.vcf.gz
        cd output
        tar -cf ${BAM_PREFIX}.intermediate_files.tar intermediate_files

    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones: zones
    }
    output {
        File vcfGz = glob("output/*.vcf.gz")[0]
        File intermediateTar = glob("output/*.tar")[0]
    }
}

