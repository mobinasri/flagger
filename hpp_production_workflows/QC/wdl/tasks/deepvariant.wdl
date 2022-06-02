version 1.0

workflow runDeepVariant {

    call DeepVariant

    output{
        File deepVariantVCF    = DeepVariant.vcfOut
        File deepVariantVCFTBI = DeepVariant.vcfIdxOut
        File deepVariantHTML   = DeepVariant.visualReport
    }
}

task DeepVariant{
    input{
        File inputReads
        File inputReadsIdx
        File assembly
        File assemblyIndex
        String sample

        String modelType = "HYBRID_PACBIO_ILLUMINA"
        String? callRegions

        Int memSizeGB = 128
        Int threadCount = 64
        Int diskSizeGB = 128
        String dockerImage = "google/deepvariant@sha256:440074e9cf854e20e3e05eab0a4fbbc32652c1f0d71b2aacbbd35da47f84faae" # 1.2.0


    }

    parameter_meta {
        inputReads: "Reads aligned to assembly. must be in BAM or CRAM format."
        inputReadsIdx: "Index file for BAM/CRAM."
        assembly: "Assembly (reference) to that reads are aligned to. If compressed, must be compressed w/ bgzip."
        assemblyIndex: "Index (fai) file for assembly"
        modelType: "modelType type to use in variant calling. Should be one of: WGS,WES,PACBIO,HYBRID_PACBIO_ILLUMINA"
    }


    String outputVCF    = "~{sample}_deepvariant.vcf.gz"

    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        ## Soft link fasta and index so they are in the same directory
        REF=$(basename ~{assembly})
        REF_IDX=$(basename ~{assemblyIndex}) 

        ln -s ~{assembly} ./$REF
        ln -s ~{assemblyIndex} ./$REF_IDX

        
        ## Soft link reads and index so they are in the same directory
        READS=$(basename ~{inputReads})
        READS_IDX=$(basename ~{inputReadsIdx}) 

        ln -s ~{inputReads} ./$READS
        ln -s ~{inputReadsIdx} ./$READS_IDX

        
        ## Pass regions argument if callRegions is set, if not just pass empty string
        if [ -z "~{callRegions}" ]
        then
            REGIONS_TOKEN=""
        else
            REGIONS_TOKEN="--regions ~{callRegions}"
        fi


        ## Run DeepVariant
        /opt/deepvariant/bin/run_deepvariant \
            --model_type ~{modelType} \
            --ref ${REF} \
            --reads ${READS} \
            --output_vcf ~{outputVCF} \
            --num_shards ~{threadCount} \
            ${REGIONS_TOKEN}

    >>>
    output{
        File vcfOut        = outputVCF
        File vcfIdxOut     = "~{outputVCF}.tbi"
        File visualReport  = glob("*visual_report.html")[0]
    }

    runtime{
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}