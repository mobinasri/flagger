version 1.0

import "flagger.wdl" as flagger_t
import "flagger_preprocess_no_variant_calling.wdl" as preprocess_t
import "../tasks/other/project_blocks_for_flagger.wdl" as project_t
import "../tasks/other/flagger_stats.wdl" as stats_t
import "../../ext/secphase/wdls/workflows/secphase.wdl" as secphase_t

workflow FlaggerEndToEndNoVariantCalling{
    input{
        File assemblyFastaGz
        File readAlignmentBam
        File hap1ToRefBam
        File hap2ToRefBam
        String secphaseDockerImage
        String secphaseOptions
        String secphaseVersion
        Float maxReadDivergence
        String sampleName
        String suffix
        String refName
        File fai
        Array[File] refBiasedBlocksBedArray
        Array[String] refBiasedRegionNameArray
        Array[Float] refBiasedRegionFactorArray
        File refSexBed
        File refSDBed
        File refCntrBed
        File refCntrCtBed
    }
    call secphase_t.runSecPhase as secphase{
        input:
            inputBam = readAlignmentBam,
            diploidAssemblyFastaGz = assemblyFastaGz,
            secphaseOptions = secphaseOptions,
            secphaseDockerImage = secphaseDockerImage,
            version = secphaseVersion
    }
    call preprocess_t.runFlaggerPreprocess as preprocess{
        input:
            bam = readAlignmentBam,
            assemblyFastaGz = assemblyFastaGz,
            phasingLogText  = secphase.outLog,
            maxDivergence = maxReadDivergence
    }
    call project_t.runProjectBlocksForFlagger as project{
        input:
            sampleName = sampleName,
            hap1AssemblyBam = hap1ToRefBam,
            hap2AssemblyBam = hap2ToRefBam,
            refSuffix = refName, 
            refBiasedBlocksBedArray = refBiasedBlocksBedArray,
            biasedBlocksNameStringArray = refBiasedRegionNameArray,
            refSexBed = refSexBed,
            refSDBed = refSDBed,
            refCntrBed = refCntrBed,
            refCntrCtBed = refCntrCtBed
    }
    call flagger_t.runFlagger as flagger{
        input:
            biasedRegionBedArray = project.projectionBiasedBedArray,
            biasedRegionNameArray = flatten([prefix("hap1_", refBiasedRegionNameArray), prefix("hap2_", refBiasedRegionNameArray)]),
            biasedRegionFactorArray = flatten([refBiasedRegionFactorArray, refBiasedRegionFactorArray]),
            coverageGz = preprocess.correctedCovGz,
            highMapqCoverageGz = preprocess.correctedHighMapqCovGz,
            fai = fai,
            sampleName = sampleName,
            suffix = suffix,
            covFloat = preprocess.meanCorrectedCoverageFloat
    }
    call stats_t.flaggerStats as stats{
        input:
            fastaGz = assemblyFastaGz,
            flaggerBed = flagger.finalBed,
            difficultBed_1 = project.projectionCntrBed,
            difficultString_1 = "Cntr",
            difficultBed_2 = project.projectionSDBed,
            difficultString_2 = "SD",
            sexBed = project.projectionSexBed,
            sample = sampleName,
            prefix = suffix
    }
    output {
        # flagger preprocess files
        Float meanCoverageFloat = preprocess.meanCorrectedCoverageFloat
        File covGz = preprocess.correctedCovGz
        File highMapqCovGz = preprocess.correctedHighMapqCovGz
        File excludedReadIdsText = preprocess.excludedReadIdsText

        # projection bed files
        Array[File] projectionBiasedBedArray = project.projectionBiasedBedArray
        File projectionSDBed = project.projectionSDBed
        File projectionSexBed = project.projectionSexBed
        File projectionCntrBed = project.projectionCntrBed

        # flagger outputs for all alignments
        File finalBed = flagger.finalBed
        File miscFilesTarGz = flagger.miscFilesTarGz
        File pdf = flagger.pdf


        # flagger statistics stratified by long contigs, centeromere, SD and sex
        # all alignments
        File statsTsv = stats.flaggerStatsTsv
        File statsPercOnlyTsv = stats.flaggerStatsPercOnlyTsv

    }
}
