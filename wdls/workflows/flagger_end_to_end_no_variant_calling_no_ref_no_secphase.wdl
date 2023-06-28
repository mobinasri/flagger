version 1.0

import "flagger.wdl" as flagger_t
import "flagger_preprocess_no_variant_calling.wdl" as preprocess_t
import "../../ext/secphase/wdls/workflows/secphase.wdl" as secphase_t

workflow FlaggerEndToEndNoVariantCallingNoRef{
    input{
        File assemblyFastaGz
        File readAlignmentBam
        String secphaseDockerImage = "mobinasri/secphase:v0.3.0"
        String secphaseOptions
        String secphaseVersion = "v0.3.0"
        Float maxReadDivergence
        String sampleName
        String suffix
        File fai
    }
    call preprocess_t.runFlaggerPreprocess as preprocess{
        input:
            bam = readAlignmentBam,
            assemblyFastaGz = assemblyFastaGz,
            maxDivergence = maxReadDivergence
    }
    call flagger_t.runFlagger as flagger{
        input:
            coverageGz = preprocess.correctedCovGz,
            highMapqCoverageGz = preprocess.correctedHighMapqCovGz,
            fai = fai,
            sampleName = sampleName,
            suffix = suffix,
            covFloat = preprocess.meanCorrectedCoverageFloat
    }
    
    output {
        # flagger preprocess files
        Float meanCoverageFloat = preprocess.meanCorrectedCoverageFloat
        File covGz = preprocess.correctedCovGz
        File highMapqCovGz = preprocess.correctedHighMapqCovGz
        File excludedReadIdsText = preprocess.excludedReadIdsText

        # flagger outputs for all alignments
        File finalBed = flagger.finalBed
        File miscFilesTarGz = flagger.miscFilesTarGz
        File pdf = flagger.pdf
    }
}
