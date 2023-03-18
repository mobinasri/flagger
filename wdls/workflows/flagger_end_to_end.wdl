version 1.0

import "flagger.wdl" as flagger_t
import "flagger_preprocess.wdl" as preprocess_t
import "../tasks/other/project_blocks_for_flagger.wdl" as project_t
import "../tasks/other/flagger_stats.wdl" as stats_t
import "../../ext/secphase/wdls/workflows/secphase.wdl" as secphase_t

workflow FlaggerEndToEnd{
    input{
        File assemblyFastaGz
        File readAlignmentBam
        File hap1ToRefBam
        File hap2ToRefBam
        String secphaseDockerImage
        String secphaseOptions
        Float maxReadDivergence
        String variantCaller 
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
            secphaseDockerImage = secphaseDockerImage
    }
    call preprocess_t.runFlaggerPreprocess as preprocess{
        input:
            bam = readAlignmentBam,
            assemblyFastaGz = assemblyFastaGz,
            phasingLogText  = secphase.outLog,
            maxDivergence = maxReadDivergence,
            variantCaller = variantCaller
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
    call flagger_t.runFlagger as flagger_alt_removed{
        input:
            biasedRegionBedArray = project.projectionBiasedBedArray,
            biasedRegionNameArray = flatten([prefix("hap1_", refBiasedRegionNameArray), prefix("hap2_", refBiasedRegionNameArray)]),
            biasedRegionFactorArray = flatten([refBiasedRegionFactorArray, refBiasedRegionFactorArray]),
            coverageGz = preprocess.altRemovedCovGz,
            highMapqCoverageGz = preprocess.altRemovedHighMapqCovGz,
            fai = fai,
            sampleName = sampleName,
            suffix = suffix + ".alt_removed",
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
    call stats_t.flaggerStats as stats_alt_removed{
        input:
            fastaGz = assemblyFastaGz,
            flaggerBed = flagger_alt_removed.finalBed,
            difficultBed_1 = project.projectionCntrBed,
            difficultString_1 = "Cntr",
            difficultBed_2 = project.projectionSDBed,
            difficultString_2 = "SD",
            sexBed = project.projectionSexBed,
            sample = sampleName,
            prefix = suffix + ".alt_removed"
    }
    output {
        # flagger preprocess files
        File outputVcfGz = preprocess.outputVcfGz
        File altBam = preprocess.altBam
        File altBai = preprocess.altBai
        Float meanCoverageFloat = preprocess.meanCorrectedCoverageFloat
        Float altRemovedMeanCoverageFloat = preprocess.meanAltRemovedCoverageFloat
        File covGz = preprocess.correctedCovGz
        File altRemovedCovGz = preprocess.altRemovedCovGz
        File highMapqCovGz = preprocess.correctedHighMapqCovGz
        File altRemovedHighMapqCovGz = preprocess.altRemovedHighMapqCovGz
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

        # flagger outputs after removing alt alignments
        File altRemovedMiscFilesTarGz = flagger_alt_removed.miscFilesTarGz
        File altRemovedPdf = flagger_alt_removed.pdf
        File altRemovedFinalBed = flagger_alt_removed.finalBed

        # flagger statistics stratified by long contigs, centeromere, SD and sex
        # all alignments
        File statsTsv = stats.flaggerStatsTsv
        File statsPercOnlyTsv = stats.flaggerStatsPercOnlyTsv

        # after removing alt alignments
        File altRemovedStatsTsv = stats_alt_removed.flaggerStatsTsv
        File altRemovedStatsPercOnlyTsv = stats_alt_removed.flaggerStatsPercOnlyTsv
    }
}
