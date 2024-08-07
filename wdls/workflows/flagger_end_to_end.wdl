version 1.0

import "../tasks/other/misc.wdl" as misc_t
import "flagger.wdl" as flagger_t
import "flagger_preprocess.wdl" as preprocess_t
import "../tasks/other/project_blocks_for_flagger.wdl" as project_t
import "../tasks/other/flagger_stats.wdl" as stats_t
import "../../ext/secphase/wdls/workflows/secphase.wdl" as secphase_t
import "../tasks/alignment/asm2asm_aligner.wdl" as asm2asm_t 
import "../tasks/coverage/bias_detector.wdl" as bias_t
import "../tasks/alignment/produce_fai.wdl" as fai_t
import "../tasks/coverage/cov2wig.wdl" as cov2wig_t

workflow FlaggerEndToEnd{
    meta {
        author: "Mobin Asri"
        email: "masri@ucsc.edu"
        description: "Running Flagger for evaluating a diploid assembly using long read alignments. More information at https://github.com/mobinasri/flagger"
    }
    parameter_meta {
        sampleName: "Sample name for example 'HG002'"
        suffix: "Suffix string that contains information about this analysis for example 'hifi_winnowmap_flagger_for_hprc'"
        hap1AssemblyFasta: "Path to uncompressed or gzip-compressed fasta file of the 1st haplotype."
        hap2AssemblyFasta: "Path to uncompressed or gzip-compressed fasta file of the 2nd haplotype."
        readAlignmentBam: "Path to sorted read alignment bam."
        readAlignmentBai: "Path to bam index for read alignment"
        maxReadDivergence: "Alignments with gap-compressed ratio higher than this will be filtered in the pre-process step. (Default: 0.1)"
        potentialBiasesBedArray: "Array of bed files each of which contains regions with potential coverage bias for example one bed file can contain HSat2 regions in haplotype 1. (Default: [])"
        sexBed: "Optional bed file containing regions assigned to X/Y chromosomes. (can be either in ref or asm coordinates)"
        SDBed: "Optional Bed file containing Segmental Duplications. (can be either in ref or asm coordinates)"
        cntrBed: "Optional Bed file containing peri/centromeric satellites (ASat, HSat, bSat, gSat) without 'ct' blocks."
        cntrCtBed: "Optional Bed file containing centromere transition 'ct' blocks."
        additionalStratificationBedArray: "Array of additional stratification bed files for final stats tsv file. (Default: [])"
        additionalStratificationNameArray: "Array of names for the stratifications provided in the argument additionalStratificationBedArray. (Default: [])"
        enableProjectingBedsFromRef2Asm: "If True it means that the given bed files are in ref coors (e.g. chm13v2) and they have to be projected to asm coors. (Default: false)"
        projectionReferenceFasta: "The given bed files are in the coordinates of this reference. A reference should be passed if enableProjectingBedsFromRef2Asm is true. (Default: '')"
        enableRunningSecphase : "If True it will run secphase in the marker mode using the parameters starting with 'secphase' otherwise skip it. (Default: false)"
        secphaseDockerImage: "Docker image for running Secphase (Default: mobinasri/secphase:v0.4.3)"
        secphaseOptions: "String containing secphase options (can be either --hifi or --ont). (Default --hifi)"
        secphaseVersion: "Secphase version. (Default: v0.4.3)"
        enableOutputtingWig: "If True it will make wig files from cov files and output them. wig files can be easily imported into IGV sessions (Default: true)"
        windowSize: "The size of the window flagger uses for finding coverage distrubutions (Default: 5000000)"
        sortPdfPagesByHaplotype: "Sort the coverage distributions plotted in the output pdf by haplotype (Default: false)"
        hap1ContigPattern: "The pattern that will be used for finding the names of the contigs belonging to haplotype1. It will be skipped if sortPdfPagesByHaplotype is false. (Default: hap1)"
        hap2ContigPattern: "The pattern that will be used for finding the names of the contigs belonging to haplotype2. It will be skipped if sortPdfPagesByHaplotype is false. (Default: hap2)"
    }
    input{
        String sampleName
        String suffix

        File hap1AssemblyFasta
        File hap2AssemblyFasta
        File readAlignmentBam
        File readAlignmentBai
        Float maxReadDivergence = 0.1

        Array[File] potentialBiasesBedArray = []
        File? sexBed
        File? SDBed
        File? cntrBed # censat annotation with no "ct"
        File? cntrCtBed
        Array[File] additionalStratificationBedArray=[]
        Array[String] additionalStratificationNameArray=[]
        
        Boolean enableProjectingBedsFromRef2Asm = true
        File projectionReferenceFasta = ""

        Boolean enableRunningSecphase = false
        String secphaseDockerImage = "mobinasri/secphase:v0.4.3--c99e0e9f3561192e127b2d913c932c3e68aa21bf"
        String secphaseOptions = "--hifi"
        String secphaseVersion = "v0.4.3"

        Boolean enableOutputtingWig = true
        Int windowSize=5000000
        Boolean sortPdfPagesByHaplotype=false
        String hap1ContigPattern="hap1"
        String hap2ContigPattern="hap2"
    }

    # Make en empty bed file
    # it will be used instead of any bed file
    # not given as an input since stats_t.flaggerStats 
    # does not handle optional bed files
    call misc_t.createFile as createEmptyBed{
        input:
            content = "",
            filename = "mock.bed"
    }

    File emptyBed = createEmptyBed.outFile

    # Create a diploid assembly 
    # from the given haplotypes
    call misc_t.createDipAsm {
        input:
            hap1AssemblyFasta = hap1AssemblyFasta,
            hap2AssemblyFasta = hap2AssemblyFasta,
            outputName = "${sampleName}.dip.asm"
    }

    # Index diploid assembly
    call fai_t.produceFai {
        input:
            assemblyGz = createDipAsm.diploidAssemblyFastaGz
    }

    # Get coordinates of canonical bases only (no "N" which may come from scaffolding)
    call misc_t.getCanonicalBasesBed as dipCanonical{
        input: 
            assemblyFasta = createDipAsm.diploidAssemblyFastaGz
    }

    # Run Secphase if it is enabled by user
    # Secphase is for fixing reads that were
    # not mapped to the correct haplotype
    # more info https://github.com/mobinasri/secphase
    if (enableRunningSecphase) {
        call secphase_t.runSecPhase as secphase{
            input:
                inputBam = readAlignmentBam,
                diploidAssemblyFastaGz = createDipAsm.diploidAssemblyFastaGz,
                secphaseOptions = secphaseOptions,
                secphaseDockerImage = secphaseDockerImage,
                version = secphaseVersion
        }
    }

    # Preprocess read alignment file by applying
    # the changes suggested by secphase and
    # then creating coverage files that can be
    # used by Flagger
    call preprocess_t.runFlaggerPreprocess as preprocess{
        input:
            bam = readAlignmentBam,
            assemblyFastaGz = createDipAsm.diploidAssemblyFastaGz,
            phasingLogText  = secphase.outLog,
            maxDivergence = maxReadDivergence,
            correctBamDockerImage = secphaseDockerImage
    }

    # Map each haplotype to the given reference
    # and then project bed files to the assembly 
    # coordinates. User can disable this part if 
    # the bed files are already in the 
    # coordinates of the assembly. Note that 
    # after projecting biased bed files it will 
    # produced one bed file per haplotype
    # but for sex, SD, cntr and cntrCt it will make
    # a combined bed file of the two haplotypes
    if (enableProjectingBedsFromRef2Asm) {
        call asm2asm_t.asm2asmAlignment as hap1ToRef {
            input :
                aligner="minimap2",
                preset="asm5",
                queryAssemblyFasta=hap1AssemblyFasta,
                refAssemblyFasta=projectionReferenceFasta
        }
        call asm2asm_t.asm2asmAlignment as hap2ToRef {
            input :
                aligner="minimap2",
                preset="asm5",
                queryAssemblyFasta=hap2AssemblyFasta,
                refAssemblyFasta=projectionReferenceFasta
        }
        call project_t.runProjectBlocksForFlagger as project{
            input:
                sampleName = sampleName,
                hap1AssemblyBam = hap1ToRef.sortedBamFile,
                hap2AssemblyBam = hap2ToRef.sortedBamFile, 
                refBiasedBlocksBedArray = potentialBiasesBedArray,
                additionalBedArray = additionalStratificationBedArray,
                refSexBed = select_first([sexBed, emptyBed]),
                refSDBed = select_first([SDBed, emptyBed]),
                refCntrBed = select_first([cntrBed, emptyBed]),
                refCntrCtBed = select_first([cntrCtBed, emptyBed])
        }
    }


    # Get bed files containing potentially biased 
    # regions in asm coordinates, which can be 
    # projections from reference
    Array[File] potentialBiasesBedArrayInAsmCoor = select_first([project.projectionBiasedBedArray, potentialBiasesBedArray])
    
    # Get bed files containing additional stratifications
    # regions in asm coordinates, which can be
    # projections from reference
    Array[File] additionalStratificationBedArrayInAsmCoor = select_first([project.projectionAdditionalBedArray, additionalStratificationBedArray])
     
    # Get bed files for SD, centromere and 
    # sex which can be projections from reference
    File SDBedInAsmCoor = select_first([project.projectionSDBed, SDBed, emptyBed])
    File cntrBedInAsmCoor = select_first([project.projectionCntrBed, cntrBed, emptyBed])
    File sexBedInAsmCoor = select_first([project.projectionSexBed, sexBed, emptyBed])
    
    # Run bias detection if there is at least one
    # bed file given as a potentially biased region
    if (length(potentialBiasesBedArray) > 0) {
        call bias_t.biasDetector {
            input:
                bedArray = potentialBiasesBedArrayInAsmCoor,
                inputBam = readAlignmentBam,
                inputBai = readAlignmentBai,
                assemblyFastaGz = createDipAsm.diploidAssemblyFastaGz,
        }
    }

    # Run flagger
    call flagger_t.runFlagger as flagger{
        input:
            biasedRegionBedArray = select_first([biasDetector.biasedRegionBedArray, []]),
            biasedRegionNameArray = select_first([biasDetector.biasedRegionNameArray, []]),
            biasedRegionFactorArray = select_first([biasDetector.biasedRegionFactorArray, []]),
            coverageGz = preprocess.correctedCovGz,
            highMapqCoverageGz = preprocess.correctedHighMapqCovGz,
            fai = produceFai.fai,
            sampleName = sampleName,
            suffix = suffix,
            covFloat = preprocess.modeCorrectedCoverageFloat,
            canonicalBasesDiploidBed = dipCanonical.canonicalBasesBed,
            windowSize = windowSize,
            sortPdfPagesByHaplotype = sortPdfPagesByHaplotype,
            hap1ContigPattern = hap1ContigPattern,
            hap2ContigPattern = hap2ContigPattern
    }

    # Get Flagger stats
    # Note that even if there is no sex, cntr and SD bed file
    # given by user it should work
    # but the numbers related to these annotations
    # will not be correct in the output tsv
    call stats_t.flaggerStats as stats{
        input:
            fastaGz = createDipAsm.diploidAssemblyFastaGz,
            flaggerBed = flagger.finalBed,
            difficultBed_1 = cntrBedInAsmCoor,
            difficultString_1 = "Cntr",
            difficultBed_2 = SDBedInAsmCoor,
            difficultString_2 = "SD",
            additionalBeds = additionalStratificationBedArrayInAsmCoor,
            additionalStrings = additionalStratificationNameArray,
            sexBed = sexBedInAsmCoor,
            minContigSize = 1000000,
            sample = sampleName,
            prefix = suffix
    }

    # make wig files from cov files
    # wig files can be easily imported into IGV sessions
    if (enableOutputtingWig) {
        call cov2wig_t.cov2wig as cov2wig{
            input:
                covGz = preprocess.correctedCovGz,
                segmentSize = 1024,
                threshold = 250,
                trackName = "${sampleName}_tot",
                fai = produceFai.fai
        }
        call cov2wig_t.cov2wig as cov2wigHighMapq{
            input:
                covGz = preprocess.correctedHighMapqCovGz,
                segmentSize = 1024,
                threshold = 250,
                trackName = "${sampleName}_high_mapq",
                fai = produceFai.fai
        }
    }

    output {
        # bed files and factors for biased regions
        Array[File] detectedBiasedRegionBedArray = select_first([biasDetector.biasedRegionBedArray, []])
        Array[Float] detectedBiasedRegionFactorArray = select_first([biasDetector.biasedRegionFactorArray, []])

        # get projected bed files if there is any
        File? projectionSexBed = project.projectionSexBed
        File? projectionSDBed = project.projectionSDBed
        File? projectionCntrBed = project.projectionCntrBed
        Array[File]? projectionAdditionalStratificationBedArray = project.projectionAdditionalBedArray

        # flagger preprocess files
        Float modeCoverageFloat = preprocess.modeCorrectedCoverageFloat
        File covGz = preprocess.correctedCovGz
        File highMapqCovGz = preprocess.correctedHighMapqCovGz
        File excludedReadIdsText = preprocess.excludedReadIdsText

        # wig files if user enabled outputting them
        File? wig = cov2wig.wig
        File? highMapqWig = cov2wigHighMapq.wig

        # secphase output
        File? secphaseOutputLog = secphase.outLog
        File? secphaseModifiedReadBlocksMarkersBed = secphase.modifiedReadBlocksMarkersBed
        File? secphaseMarkerBlocksBed = secphase.markerBlocksBed

        # flagger outputs for all alignments
        File finalBed = flagger.finalBed
        File finalBedNoHap = flagger.finalBedNoHap
        File miscFilesTarGz = flagger.miscFilesTarGz
        File pdf = flagger.pdf

        # flagger statistics stratified by long contigs, centeromere, SD and sex
        # all alignments
        File statsTsv = stats.flaggerStatsTsv
        File statsPercOnlyTsv = stats.flaggerStatsPercOnlyTsv
    }
}

