version 1.0

workflow runDeepVariant{
    input {
        File assemblyFastaGz
        File bam
        File? bamIndex
        String deepVariantModelType
        Int minMAPQ = 0
        String includeSecondary="False"
        String includeSupplementary="True"
        String openvinoString="use_openvino=false" # for deepvariant_1.4 it should be false
    }
    call deepVariant{
        input:
            modelType = deepVariantModelType,
            assemblyFastaGz = assemblyFastaGz,
            bam = bam,
            bamIndex = bamIndex,
            includeSecondary = includeSecondary,
            includeSupplementary = includeSupplementary,
            minMAPQ = minMAPQ,
            threadCount=64,
            memSize=256,
            diskSize= 2 * ceil(size(bam, "GB")) + 64,
            openvinoString=openvinoString
    }
    output{
        File vcfGz = deepVariant.vcfGz
    }
}


task deepVariant{
    input{
        File bam
        File? bamIndex
        File assemblyFastaGz
        File? bed
        Int minMAPQ
        String includeSecondary="False"
        String includeSupplementary="False"
        String modelType
        String openvinoString="openvino=false"
        # runtime configurations
        Int memSize=32
        Int threadCount=16
        Int diskSize=64
        String dockerImage="google/deepvariant:latest"
        Int preemptible=2
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
        
        ## hard link the bam file to the working directory and produce its index file
        BAM_NAME=$(basename ~{bam})
        BAM_PREFIX=${BAM_NAME%%.bam}
        ln -f ~{bam} > ${BAM_PREFIX}.bam
        if [ -n "~{bamIndex}" ]; then
            ln -f ~{bamIndex} > ${BAM_PREFIX}.bam.bai
        else       
            samtools index -@~{threadCount} ${BAM_PREFIX}.bam
        fi

        ## unzip the fasta file and produce its index
        gunzip -c ~{assemblyFastaGz} > asm.fa
        samtools faidx asm.fa

        MAKE_EXAMPLES_EXTRA_ARGS="min_mapping_quality=~{minMAPQ}"
        if [ ~{includeSecondary} == "True"]; then
            MAKE_EXAMPLES_EXTRA_ARGS="${MAKE_EXAMPLES_EXTRA_ARGS},keep_secondary_alignments=true"
        fi
        if [ ~{includeSupplementary} == "True" ]; then
            MAKE_EXAMPLES_EXTRA_ARGS="${MAKE_EXAMPLES_EXTRA_ARGS},keep_supplementary_alignments=true"
        fi

        MORE_OPTIONS=""
        if [ -n "~{bed}" ]; then
            MORE_OPTIONS="--regions=~{bed}"
        fi

        ## call deepvariant 
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=~{modelType} \
        --ref=asm.fa \
        --reads=${BAM_PREFIX}.bam \
        --output_vcf=${BAM_PREFIX}.vcf \
        --make_examples_extra_args=${MAKE_EXAMPLES_EXTRA_ARGS} \
        --call_variants_extra_args=~{openvinoString} \
        --num_shards=$(nproc) \
        --dry_run=false ${MORE_OPTIONS}

        gzip -c ${BAM_PREFIX}.vcf > ${BAM_PREFIX}.vcf.gz 
        
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
        File vcfGz = glob("*.vcf.gz")[0]
    }
}

