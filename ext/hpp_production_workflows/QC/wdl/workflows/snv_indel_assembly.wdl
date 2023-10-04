version 1.0

import "../tasks/deepvariant.wdl" as runDeepVariant
import "../tasks/pepperMarginDeepVariant.wdl" as runPepperMarginDeepVariant
import "../tasks/hapDotPy.wdl" as runHappy


workflow snv_indel_assembly {

    meta {
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
        description: "WDL implementation of the small variant calling section of the [T2T Polishing Case Study](https://github.com/arangrhie/T2T-Polish/blob/master/doc/T2T_polishing_case_study.md)."
    }

    input {
        File assembly
        File assemblyIndex
        String sample
    }

    ## Call DeepVariant (recommended input of HiFi + Ilmn data)
    call runDeepVariant.DeepVariant as deepVariant_t {
        input:
            assembly      = assembly,
            assemblyIndex = assemblyIndex,
            sample        = sample
    }

    ## Call PepperMarginDeepVariant (ONT data)
    call runPepperMarginDeepVariant.pepperMarginDeepVariant as pmdv_t {
        input:
            assembly      = assembly,
            assemblyIndex = assemblyIndex,
            sample        = sample
    }    

    ## Filter DeepVariant calls
    call runPepperMarginDeepVariant.bcftoolsFilter as filt_1 {
        input:
            inputVCF      = deepVariant_t.vcfOut,
            excludeExpr   = "'FORMAT/VAF<=0.5 | FORMAT/GQ<=30'",
            applyFilters  = "PASS"
    }

    ## Filter Pepper Margin DeepVariant calls
    call runPepperMarginDeepVariant.bcftoolsFilter as filt_2 {
        input:
            inputVCF      = pmdv_t.vcfOut,
            excludeExpr   = "'FORMAT/VAF<=0.5 | FORMAT/GQ<=30'",
            applyFilters  = "PASS",
            exludeTypes   = "indels"
    }

    ## Compare filtered callsets (DeepVariant & PMDV) to see where they agree
    call runHappy.hapDotPy as happy_t {
        input:    
            truthVCF      = filt_1.vcfOut,
            queryVCF      = filt_2.vcfOut,
            assembly      = assembly,
            assemblyIndex = assemblyIndex,
            sample        = sample        
    }

    ## Output the union of the two filtered callsets (DeepVariant & PMDV)
    call mergeVCF {
        input:
            VCF1     = filt_1.vcfOut,
            VCF2     = filt_2.vcfOut,
            happyVCF = happy_t.vcfOut,
            sample   = sample
    }

    ## create text file w/ stats of final output from the last step
    call createVCFStats {
        input:
            inputVCF = mergeVCF.vcfOut
    }


    output{
        File deepVariantVCF    = deepVariant_t.vcfOut
        File deepVariantVCFTBI = deepVariant_t.vcfIdxOut
        File deepVariantHTML   = deepVariant_t.visualReport

        File pepperMarginDeepVariantVCF     = pmdv_t.vcfOut
        File pepperMarginDeepVariantVCFIdx  = pmdv_t.vcfIdxOut
        File pepperMarginDeepVariantHTML    = pmdv_t.visualReport

        File filtDeepVarVCF = filt_1.vcfOut
        File filtPmdvVCF    = filt_2.vcfOut

        File hapDotPyVCF    = happy_t.vcfOut
        File hapDotPyVCFIdx = happy_t.vcfIdxOut
        File hapDotPyTar    = happy_t.happyTar

        File mergedVCF      = mergeVCF.vcfOut
        File mergedVCFIdx   = createVCFStats.vcfTBIOut
        File mergedVCFStats = createVCFStats.statsFileOut
    }
}

task mergeVCF {
    input{
        File VCF1
        File VCF2
        File happyVCF
        String sample

        Int memSize = 8
        Int threadCount = 2
        Int diskSizeGB = 128
        String dockerImage = "kishwars/t2t_polishing@sha256:418486a1e88c48555ad4f7158c0a9923762182e7c9cd883342ffe0a161d89de6" # v0.1


    }

    String outputVCF = "~{sample}_MERGED_SMALL_VARIANTS.vcf.gz"

    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        ln -s /opt/vcf_merge_t2t.py .

        python3 vcf_merge_t2t.py \
            -v1 ~{VCF1} \
            -v2 ~{VCF2} \
            -hv ~{happyVCF} \
            -o ~{outputVCF}

    >>>
    output{
        File vcfOut = outputVCF
    }

    runtime{
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}

task createVCFStats {
    input{
        File inputVCF

        Int memSize = 8
        Int threadCount = 2
        Int diskSizeGB = 128
        String dockerImage = "kishwars/pepper_deepvariant@sha256:70908591ad67e8567a6e4551119b2cfc33d957ad39701c8af51b36b516214645" # r0.8


    }

    String filePrefix   = basename(inputVCF, ".vcf.gz")
    String statsFile    = "~{filePrefix}.stats.txt"

    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        
        ## copy file so index is in local directory (easier for delocalization)
        INPUT_VCF="~{filePrefix}.vcf.gz"
        cp ~{inputVCF} $INPUT_VCF
        
        ## Index VCF (not related to stats, just nice to have)
        tabix -p vcf $INPUT_VCF

        ## call bcftools stats
        bcftools stats $INPUT_VCF > ~{statsFile}
    >>>
    output{
        File vcfTBIOut    = glob("~{filePrefix}*.tbi")[0]
        File statsFileOut = statsFile
    }

    runtime{
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}