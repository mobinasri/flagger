version 1.0

import "../../ext/secphase/wdls/workflows/correct_bam.wdl" as correct_bam_t
import "../tasks/coverage/bam_coverage.wdl" as bam_coverage_t

workflow runFlaggerPreprocess{
    input {
        File bam
        File assemblyFastaGz
        File? phasingLogText
        Int minReadLength = 5000
        Int minAlignmentLength = 5000
        Float maxDivergence # For HiFi 0.02, For ONT guppy 6 -> 0.09
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
    

    
    ## Calculate coverage for the corrected bam file (without filtering)
    call bam_coverage_t.bamCoverageFast as bam2cov_corrected{
        input:
            bam = correctBam.correctedBam,
            bai = correctBam.correctedBamIndex,
            output_format = "only_total",
            minMAPQ = 0,
            assemblyFastaGz = assemblyFastaGz
    }
    

    ## Calculate coverage of reads with high mapqs (> 20) for the 
    ## corrected bam file (without filtering)
    ##
    ## This coverage will be used for checking the false duplications
    ## in the 2nd phase of the FLAGGER pipeline
    call bam_coverage_t.bamCoverageFast as bam2cov_corrected_highMapq{
        input:
            bam = correctBam.correctedBam,
            bai = correctBam.correctedBamIndex,
            output_format = "only_high_mapq",
            minMAPQ = 20,
            assemblyFastaGz = assemblyFastaGz
    }

    output {
        Float meanCorrectedCoverageFloat = bam2cov_corrected.coverageMean
        File correctedCovGz = bam2cov_corrected.coverageGz
        File correctedHighMapqCovGz = bam2cov_corrected_highMapq.coverageGz
        File excludedReadIdsText = correctBam.excludedReadIdsText
    }
}
