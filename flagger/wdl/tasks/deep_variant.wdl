version 1.0

import "../../../asset/wdl/tasks/variant_calling.wdl" as var_t

workflow runVariantCalling{
    input {
        File assemblyFastaGz
        File bam
        File bamIndex
        String deepVariantModelType
        Int minMAPQ = 0
        Int numberOfCallerNodes=16
        Int nodeThreadCount=16
        String includeSecondary="False"
        String includeSupplementary="False"
    }
    call splitBamContigWise {
        input:
            assemblyFastaGz = assemblyFastaGz,
            bam = bam,
            bamIndex = bamIndex,
            splitNumber = numberOfCallerNodes,
            threadCount = numberOfCallerNodes,
            diskSize = 2 * ceil(size(bam, "GB")) + 64
    }
    scatter (part in zip(splitBamContigWise.splitBams, splitBamContigWise.splitBeds)) {
        #call increaseMapq{
        #    input:
        #       bam = part.left
        #}
        call deepVariant{
            input:
                modelType = deepVariantModelType,
                assemblyFastaGz = assemblyFastaGz,
                bam = part.left,
                bed = part.right,
                includeSecondary = includeSecondary,
                includeSupplementary = includeSupplementary,
                minMAPQ = minMAPQ,
                diskSize = 2 * ceil(size(part.left, "GB")) + 64
        }
    }
    call var_t.mergeVcf{
        input:
            vcfGzFiles = deepVariant.vcfGz,
            outputName = basename("${bam}", ".bam")
    }
    output{
        File vcfGz = mergeVcf.vcfGz 
    }
}

task splitBamContigWise{
    input{
        File assemblyFastaGz
        File bam
        File bamIndex
        Int splitNumber
        Int memSize=32
        Int threadCount
        Int diskSize=512
        String dockerImage="quay.io/masri2019/hpp_bcftools:latest"
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

        ## unzip fasta file and produce its index file        
        ASSEMBLY_NAME=$(basename ~{assemblyFastaGz})
        ASSEMBLY_PREFIX=${ASSEMBLY_NAME%%.fa.gz}
        gunzip -c ~{assemblyFastaGz} > ${ASSEMBLY_PREFIX}.fa
        samtools faidx ${ASSEMBLY_PREFIX}.fa

        ## hard link the bam and bai files to the working directory
        BAM_NAME=$(basename ~{bam})
        BAM_PREFIX=${BAM_NAME%%.bam}
        ln -f ~{bam} > ${BAM_PREFIX}.bam
        ln -f ~{bamIndex} > ${BAM_PREFIX}.bam.bai

        ## make a bed file that covers the whole assembly
        cat ${ASSEMBLY_PREFIX}.fa.fai | awk '{print $1"\t"0"\t"$2}' > ${ASSEMBLY_PREFIX}.bed

        ## split the bed file of the whole assembly into multiple bed files
        mkdir split_beds split_bams
        python3 ${SPLIT_BED_CONTIG_WISE_PY} --bed ${ASSEMBLY_PREFIX}.bed --n ~{splitNumber} --dir split_beds --prefix ${ASSEMBLY_PREFIX}
        
        ## make a separate bam for each bed file
        n=$(ls split_beds/ | wc -l)
        seq 1 ${n} | xargs -I {} -n 1 -P ~{threadCount} sh -c "samtools view -h -b -L split_beds/${ASSEMBLY_PREFIX}_{}.bed ${BAM_PREFIX}.bam > split_bams/${BAM_PREFIX}_{}.bam"
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
        Array[File] splitBams = glob("split_bams/*.bam")
        Array[File] splitBeds = glob("split_beds/*.bed")
    }
}


task increaseMapq{
    input{
        File bam
        Int threshold=20
        Int memSize=4
        Int threadCount=2
        Int diskSize=64
        String dockerImage="quay.io/masri2019/hpp_coverage:latest"
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

        ## hard link the bam and bai files to the working directory
        BAM_NAME=$(basename ~{bam})
        BAM_PREFIX=${BAM_NAME%%.bam}

        mkdir output
        increase_mapq -i ~{bam} -o output/${BAM_PREFIX}.increased_mapq.bam -t ~{threshold}
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
        File outputBam = glob("output/*.bam")[0]
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
        --call_variants_extra_args="use_openvino=true" \
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

