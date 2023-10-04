version 1.0

workflow runpepperMarginDeepVariant {

    call pepperMarginDeepVariant

    call bcftoolsFilter {
        input:
            inputVCF = pepperMarginDeepVariant.vcfOut
    }
    output {
        File pepperMarginDeepVariantVCF     = pepperMarginDeepVariant.vcfOut
        File pepperMarginDeepVariantVCFIdx  = pepperMarginDeepVariant.vcfIdxOut
        File pepperMarginDeepVariantHTML    = pepperMarginDeepVariant.visualReport

        File filtPepperMarginDeepVariantVCF = bcftoolsFilter.vcfOut
    }
}


task pepperMarginDeepVariant {
    input {
        File inputReads
        File inputReadsIdx
        File assembly
        File assemblyIndex
        String sample

        String readTypeFlag = "ont_r9_guppy5_sup"
        String? extraArgs

        Int memSize = 128
        Int threadCount = 64
        Int diskSizeGB = 128
        String dockerImage = "kishwars/pepper_deepvariant@sha256:70908591ad67e8567a6e4551119b2cfc33d957ad39701c8af51b36b516214645" # r0.8
    
    }

    parameter_meta {
        inputReads: "Reads aligned to assembly. must be in BAM format."
        inputReadsIdx: "Index file for BAM."
        assembly: "Assembly (reference) to that reads are aligned to."
        assemblyIndex: "Index (fai) file for assembly"
        sample: "Sample name. Will be used in output VCF file."
        readTypeFlag: "Read type flag to pass to pepperMarginDeepVariant. See documentation for allowable values."
    }


    String outputPrefix = "~{sample}_PEPPER_DeepVariant"

    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        ## Soft link fasta and index so they are in the same directory
        REF=$(basename ~{assembly})
        REF_IDX=$(basename ~{assemblyIndex}) 

        ln -s ~{assembly} ./$REF
        ln -s ~{assemblyIndex} ./$REF_IDX

        
        ## Soft link fasta and index so they are in the same directory
        READS=$(basename ~{inputReads})
        READS_IDX=$(basename ~{inputReadsIdx}) 

        ln -s ~{inputReads} ./$READS
        ln -s ~{inputReadsIdx} ./$READS_IDX

        ## Pass optional arguments if extraArgs is set, if not just pass empty string
        if [ -z "~{extraArgs}" ]
        then
            EXTRA_ARGS=""
        else
            EXTRA_ARGS="~{extraArgs}"
        fi


        ## Run pepperMarginDeepVariant
        run_pepper_margin_deepvariant call_variant \
            -b ${READS} \
            -f ${REF} \
            -p "~{outputPrefix}" \
            -s "~{sample}" \
            -o "pepper_deepvariant_output" \
            -t "~{threadCount}" \
            --"~{readTypeFlag}" \
            ${EXTRA_ARGS}

    >>>

    output {
        File vcfOut        = glob("pepper_deepvariant_output/~{outputPrefix}*")[0]
        File vcfIdxOut     = glob("pepper_deepvariant_output/~{outputPrefix}*.tbi")[0]
        File visualReport  = glob("pepper_deepvariant_output/~{outputPrefix}*visual_report.html")[0]
    }

    runtime {
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 2
    }
}


task bcftoolsFilter {
    input {
        File inputVCF

        String excludeExpr   = "'FORMAT/VAF<=0.5 | FORMAT/GQ<=30'"
        String applyFilters  = "PASS"
        String exludeTypes   = ""

        Int memSize = 8
        Int threadCount = 4
        Int diskSizeGB = 50
        String dockerImage = "kishwars/pepper_deepvariant@sha256:70908591ad67e8567a6e4551119b2cfc33d957ad39701c8af51b36b516214645" # r0.8
    
    }

    parameter_meta {
        inputVCF: "VCF to filter. Must be gzipped."
    }

    String inputPrefix = basename(inputVCF, ".vcf.gz")
    String outputFile  = "~{inputPrefix}.filtered.vcf.gz"

    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace


        ## Pass optional argument if applyFilters is set, if not just pass empty string
        if [ -z "~{applyFilters}" ]
        then
            APPLY_FILTERS_TOKEN=""
        else
            APPLY_FILTERS_TOKEN="--apply-filters ~{applyFilters}"
        fi

        ## Pass optional argument if exludeTypes is set, if not just pass empty string
        if [ -z "~{exludeTypes}" ]
        then
            EXCLUDE_TYPES_TOKEN=""
        else
            EXCLUDE_TYPES_TOKEN="--exclude-types ~{exludeTypes}"
        fi
        

        ## Call bcftools
        bcftools view \
            $APPLY_FILTERS_TOKEN \
            $EXCLUDE_TYPES_TOKEN \
            -e ~{excludeExpr} \
            -Oz \
            ~{inputVCF} \
            > ~{outputFile}
    >>>

    output {
        File vcfOut = outputFile
    }

    runtime {
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 2
    }
}