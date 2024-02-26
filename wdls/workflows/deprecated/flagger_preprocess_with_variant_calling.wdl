version 1.0

import "../../ext/secphase/wdls/workflows/correct_bam.wdl" as correct_bam_t
import "../tasks/alignment/filter_alt_reads.wdl" as filter_alt_reads_t
import "../tasks/coverage/bam_coverage.wdl" as bam_coverage_t
import "../tasks/variant_calling/deep_variant_scattered.wdl" as dv_scat_t
import "../tasks/variant_calling/pepper_margin_deep_variant_scattered.wdl" as pmdv_scat_t

workflow runFlaggerPreprocess{
    input {
        File bam
        File assemblyFastaGz
        File? phasingLogText
        Int minReadLength = 5000
        Int minAlignmentLength = 5000
        Float maxDivergence # For HiFi 0.02, For ONT guppy 6 -> 0.09
        String deepVariantModelType = "PACBIO"
        String pepperModelType = "ont_r9_guppy5_sup"
        String variantCaller  # For HiFi dv, for ONT pmdv
        String pmdvDockerImage = "kishwars/pepper_deepvariant:r0.8"
        String dvDockerImage = "google/deepvariant:1.4.0"
        Float vafCutoff = 0.3 # for filterAltReads
        Int qCutoff = 10 # for filterAltReads
        String moreOptions = "-m1000 -r0.4" # for filterAltReads
        Int variantCallingMemory = 48
    }
    
    ## Correct the bam file by swapping pri/sec tags for the wrongly phased reads
    call correct_bam_t.correctBam {
        input:
            bam = bam,
            phasingLogText = phasingLogText,
            suffix = "corrected",
            options = "--primaryOnly --minReadLen ${minReadLength} --minAlignment ${minAlignmentLength} --maxDiv ${maxDivergence}",
            flagRemoveSupplementary = true,
            flagRemoveMultiplePrimary = true
    }
    
    ## If the user selected deepvariant as the variant caller
    if ("${variantCaller}" == "dv") { 
        ## Call variants to be used for finding the reads with alternative alleles
        call dv_scat_t.runDeepVariantScattered as dv {
            input:
                deepVariantModelType = deepVariantModelType,
                assemblyFastaGz = assemblyFastaGz,
                bam = correctBam.correctedBam,
                bamIndex = correctBam.correctedBamIndex,
                minMAPQ = 0,
                includeSecondary="False",
                includeSupplementary="False",
                dockerImage = dvDockerImage,
                variantCallingMemory = variantCallingMemory
        }
    }

    ## If the user selected pepper-margin-deepvariant as the variant caller 
    if ("${variantCaller}" == "pmdv") {
        ## Call variants to be used for finding the reads with alternative alleles
        call pmdv_scat_t.runPepperMarginDeepVariantScattered as pmdv {
            input:
                pmdvModelType = pepperModelType,
                assemblyFastaGz = assemblyFastaGz,
                bam = correctBam.correctedBam,
                bamIndex = correctBam.correctedBamIndex,
                minMAPQ = 0,
                includeSupplementary="False",
                flagRemoveMultiplePrimary = false,
                dockerImage = pmdvDockerImage,
                variantCallingMemory = variantCallingMemory
        }
    }
   
    File vcfGz = select_first([dv.vcfGz, pmdv.vcfGz])     

    ## Filter the reads with alternative alleles
    call filter_alt_reads_t.filterAltReads {
        input:
            vcf = vcfGz,
            bam = correctBam.correctedBam,
            moreOptions = moreOptions,
            vafCutoff = vafCutoff,
            qCutoff = qCutoff
    }
    
    ## Calculate coverage for the corrected bam file (without filtering)
    call bam_coverage_t.bamCoverageFast as bam2cov_corrected{
        input:
            bam = correctBam.correctedBam,
            minMAPQ = 0,
            assemblyFastaGz = assemblyFastaGz
    }
    
    ## Calculate coverage for the corrected bam file in which the 
    ## reads with alternative alleles are removed
    call bam_coverage_t.bamCoverageFast as bam2cov_altRemoved{
        input:
            bam = filterAltReads.filteredBam,
            minMAPQ = 0,
            assemblyFastaGz = assemblyFastaGz
    }

    ## Calculate coverage of reads with high mapqs (20<) for the 
    ## corrected bam file (without filtering)
    ##
    ## This coverage will be used for checking the false duplications
    ## in the 2nd phase of the FLAGGER pipeline
    call bam_coverage_t.bamCoverageFast as bam2cov_corrected_highMapq{
        input:
            bam = correctBam.correctedBam,
            minMAPQ = 20,
            assemblyFastaGz = assemblyFastaGz
    }

    ## Calculate coverage of reads with high mapqs (20<) for the 
    ## corrected bam file in which the
    ## reads with alternative alleles are removed
    ##
    ## This coverage will be used for checking the false duplications
    ## in the 2nd phase of the FLAGGER pipeline
    call bam_coverage_t.bamCoverageFast as bam2cov_altRemoved_highMapq{
        input:
            bam = filterAltReads.filteredBam,
            minMAPQ = 20,
            assemblyFastaGz = assemblyFastaGz
    }
    output {
        File outputVcfGz = vcfGz
        File altBam = filterAltReads.altBam
        File altBai = filterAltReads.altBamIndex
        Float meanCorrectedCoverageFloat = bam2cov_corrected.coverageMean
        Float meanAltRemovedCoverageFloat = bam2cov_altRemoved.coverageMean
        File correctedCovGz = bam2cov_corrected.coverageGz
        File altRemovedCovGz = bam2cov_altRemoved.coverageGz
        File correctedHighMapqCovGz = bam2cov_corrected_highMapq.coverageGz
        File altRemovedHighMapqCovGz = bam2cov_altRemoved_highMapq.coverageGz
        File excludedReadIdsText = correctBam.excludedReadIdsText
    }
}
