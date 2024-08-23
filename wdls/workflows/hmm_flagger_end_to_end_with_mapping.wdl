version 1.0

import "../tasks/other/misc.wdl" as misc_t
import "hmm_flagger_end_to_end.wdl" as flagger_t
import "long_read_aligner_scattered.wdl" as aligner_t



workflow FlaggerEndToEndWithMapping{
    meta {
        author: "Mobin Asri"
        email: "masri@ucsc.edu"
        description: "Mapping reads and running Flagger for evaluating a diploid assembly using long read alignments.More information at https://github.com/mobinasri/flagger"
    }
    parameter_meta {
        sampleName: "Sample name for example 'HG002'"
        suffixForFlagger: "Suffix string that contains information about this analysis for example 'hifi_winnowmap_flagger_for_hprc'"
        suffixForMapping: "Suffix string that contains information about this alignment. It will be appended to the name of the final alignment. For example 'hifi_winnowmap_v2.03_hprc_y2'"
        hap1AssemblyFasta: "Path to uncompressed or gzip-compressed fasta file of the 1st haplotype."
        hap2AssemblyFasta: "Path to uncompressed or gzip-compressed fasta file of the 2nd haplotype."
        readfiles: "An array of read files. Their format can be either fastq, fq, fastq.gz, fq.gz, bam or cram. For cram format referenceFastaForReadExtraction should also be passed."
        aligner: "Name of the aligner. It can be either minimap2, winnowmap or veritymap. (Default = winnowmap)"
        preset: "Paremeter preset should be selected based on aligner and sequencing platform. Common presets are map-pb/map-hifi/map-ont for minimap2, map-pb/map-ont for winnowmap and hifi-haploid/hifi-haploid-complete/hifi-diploid/ont-haploid-complete for veritymap"
        kmerSize: "The kmer size for using minimap2 or winnowmap. With winnowmap kmer size should be 15 and with minimap2 kmer size should be 17 and 19 for using the presets map-ont and map-hifi/map-pb respectively."
        alignerOptions: "Aligner options. It can be something like '--eqx --cs -Y -L -y' for minimap2/winnowmap. Note that if assembly is diploid and aligner is either minimap2 or winnowmap '-I8g' is necessary. If the reads contain modification tags and these tags are supposed to be present in the final alignment file, alignerOptions should contain '-y' and the aligner should be either minimap2 or winnowmap. If running secphase is enabled it is recommended to add '-p0.5' to alignerOptions; it will keep more secondary alignments so secphase will have more candidates per read. For veritymap '--careful' can be used but not recommended for whole-genome assembly since it increases the runtime dramatically."
        readExtractionOptions: "The options to be used while converting bam to fastq with samtools fastq. If the reads contain epigentic modification tags it is necessary to use '-TMm,Ml,MM,ML'"
        referenceFastaForReadExtraction: "If reads are in CRAM format then the related reference should be passed to this parameter for read extraction."
        enableAddingMDTag: "If true it will call samtools calmd to add MD tags to the final bam file. [Default = true]"
        splitNumber: "The number of chunks which the input reads should be equally split into. Note that enableSplittingReadsEqually should be set to true if user wants to split reads into equally sized chunks. [Default = 16]"
        enableSplittingReadsEqually: "If true it will merge all reads together and then split them into multiple chunks of roughly equal size. Each chunk will then be aligned via a separate task. This feature is useful for running alignment on cloud/slurm systems where there are  multiple nodes available with enough computing power and having alignments tasks distributed among small nodes is more efficient or cheaper than running a single alignment task in a large node. If the  whole workflow is being on a single node it is not recommened to use this feature since mergin and splitting reads takes its own time. [Default = false]"
        minReadLengthForMapping: "If it is greater than zero, a task will be executed for filtering reads shorter than this value before alignment. [Default = 0]"
        alignerThreadCount : "The number of threads for mapping in each alignment task [Default = 16]"
        alignerMemSize : "The size of the memory in Gb for mapping in each alignment task [Default = 48]"
        alignerDockerImage : "The mapping docker image [Default = 'mobinasri/long_read_aligner:v0.4.0']"
    }

    input{
        String sampleName
        String suffixForMapping
        String suffixForFlagger
        File hap1AssemblyFasta
        File hap2AssemblyFasta
        Array[File] readFiles
        String presetForMapping

        String aligner="winnowmap"
        Int kmerSize = 15
        String alignerOptions="--eqx -Y -L -y"
        String readExtractionOptions="-TMM,ML,Mm,Ml"
        File? referenceFastaForReadExtraction
        Boolean enableAddingMDTag=true
        Int splitNumber = 16
        Boolean enableSplittingReadsEqually=false

        Int minReadLengthForMapping = 0
        String? correctBamOptions = "--primaryOnly --minReadLen 5000 --minAlignment 5000 --maxDiv 0.1" 
        Int alignerThreadCount = 16
        Int alignerMemSize = 48
        String alignerDockerImage = "mobinasri/long_read_aligner:v0.4.0"
        Float downSamplingRateForFlagger = 1.0

        File? includeContigListText
        File? binSizeArrayTsv
        Int chunkLen = 20000000
        Int windowLen = 4000
        String labelNames = "Err,Dup,Hap,Col"
        String trackName = "hmm_flagger_v1.0"
        Int numberOfIterationsForFlagger = 100
        Float convergenceToleranceForFlagger = 0.001
        Float maxHighMapqRatio=0.25
        String? flaggerMoreOptions
        File? alphaTsv
        String modelType = "gaussian"
        Array[Int] flaggerMinimumBlockLenArray = []
        Int flaggerMemSize=32
        Int flaggerThreadCount=8
        String flaggerDockerImage="mobinasri/flagger:v1.0.0_alpha"

        File? sexBed
        File? SDBed
        File? cntrBed # censat annotation with no "ct"
        File? cntrCtBed

        File? sexBedToBeProjected
        File? SDBedToBeProjected
        File? cntrBedToBeProjected # censat annotation with no "ct"
        File? cntrCtBedToBeProjected

        Array[File] annotationsBedArray = []
        Array[File] annotationsBedArrayToBeProjected = []
        Array[File] biasAnnotationsBedArray = []
        Array[File] biasAnnotationsBedArrayToBeProjected = []

        File? projectionReferenceFasta

        Boolean enableRunningSecphase = false
        String secphaseDockerImage = "mobinasri/secphase:v0.4.3--c99e0e9f3561192e127b2d913c932c3e68aa21bf"
        String secphaseOptions = "--hifi"
        String secphaseVersion = "v0.4.3"

        Boolean enableOutputtingBigWig = true
        Boolean enableOutputtingBam = true
        File? truthBedForMisassemblies
        
    }

    # Create a diploid assembly
    # from the given haplotypes
    call misc_t.createDipAsm {
        input:
            hap1AssemblyFasta = hap1AssemblyFasta,
            hap2AssemblyFasta = hap2AssemblyFasta,
            outputName = "${sampleName}.dip.asm"
    }

    call aligner_t.longReadAlignmentScattered as alignmentTask{
        input:
            readFiles = readFiles,
            assemblyFasta = createDipAsm.diploidAssemblyFastaGz,
            aligner = aligner,
            preset = presetForMapping,
            kmerSize = kmerSize,
            alignerOptions = alignerOptions,
            readExtractionOptions = readExtractionOptions,
            sampleName = sampleName,
            suffix = suffixForMapping,
            referenceFastaForReadExtraction = referenceFastaForReadExtraction,
            enableAddingMDTag = enableAddingMDTag,
            splitNumber = splitNumber,
            enableSplittingReadsEqually = enableSplittingReadsEqually,

            # minReadLength will be used for filtering
            # short reads before mapping reads
            minReadLength = minReadLengthForMapping,

            alignerThreadCount = alignerThreadCount,
            alignerMemSize = alignerMemSize,
            alignerDockerImage = alignerDockerImage,

            enableRunningSecphase = enableRunningSecphase,
            secphaseOptions = secphaseOptions,
            secphaseDockerImage = secphaseDockerImage,
            secphaseVersion = secphaseVersion,
            correctBamOptions = correctBamOptions,
    }

    ## To filter alignments use correctBamOptions above to void running correctBam twice
    ## once in the alignment workflow and once in the HMM-Flagger workflow
    ## maxReadDivergence = 1.0, # set to 1.0 to deactivate running correctBam
    ## minReadLength = 0, # set to 0 to deactivate running correctBam
    ## minAlignmentLength = 0, # set to 0 to deactivate running correctBam
    call flagger_t.HMMFlaggerEndToEnd as hmmFlaggerTask{
        input:
            sampleName = sampleName,
            suffix = suffixForFlagger,

            hap1AssemblyFasta = hap1AssemblyFasta,
            hap2AssemblyFasta = hap2AssemblyFasta,
            readAlignmentBam = alignmentTask.bamFile,
            readAlignmentBai = alignmentTask.baiFile,
            maxReadDivergence = 1.0, # set to 1.0 to deactivate running correctBam
            minReadLength = 0, # set to 0 to deactivate running correctBam
            minAlignmentLength = 0, # set to 0 to deactivate running correctBam
            downSamplingRate = 1.0,
            includeContigListText = includeContigListText,
            binArrayTsv = binSizeArrayTsv,
            chunkLen = chunkLen,
            windowLen = windowLen,
            labelNames = labelNames,
            trackName = trackName,
            numberOfIterations = numberOfIterationsForFlagger,
            convergenceTolerance = convergenceToleranceForFlagger,
            maxHighMapqRatio = maxHighMapqRatio,
            flaggerMoreOptions = flaggerMoreOptions,
            alphaTsv = alphaTsv,
            modelType = modelType,
            flaggerMinimumBlockLenArray = flaggerMinimumBlockLenArray,
            flaggerMemSize = flaggerMemSize,
            flaggerThreadCount = flaggerThreadCount,
            flaggerDockerImage = flaggerDockerImage,

            sexBed = sexBed,
            SDBed = SDBed,
            cntrBed = cntrBed, # censat annotation with no "ct"
            cntrCtBed = cntrCtBed,

            sexBedToBeProjected = sexBedToBeProjected,
            SDBedToBeProjected = SDBedToBeProjected,
            cntrBedToBeProjected = cntrBedToBeProjected, # censat annotation with no "ct"
            cntrCtBedToBeProjected = cntrCtBedToBeProjected,

            annotationsBedArray = annotationsBedArray,
            annotationsBedArrayToBeProjected = annotationsBedArrayToBeProjected,
            biasAnnotationsBedArray = biasAnnotationsBedArray,
            biasAnnotationsBedArrayToBeProjected = biasAnnotationsBedArrayToBeProjected,

            projectionReferenceFasta = projectionReferenceFasta,

            enableRunningSecphase = false,
            enableOutputtingBigWig = enableOutputtingBigWig,
            truthBedForMisassemblies = truthBedForMisassemblies,
    }
    if (enableOutputtingBam){
        File readAlignmentBamToOutput = alignmentTask.bamFile
        File readAlignmentBaiToOutput = alignmentTask.baiFile
    }
    output {
        # HMM-Flagger coverage related files
        File coverageGz = hmmFlaggerTask.coverageGz
        File biasTableTsv = hmmFlaggerTask.biasTableTsv

        # HMM-Flagger output bed files
        File finalPredictionBed = hmmFlaggerTask.finalBed
        File intermediatePredictionBed = hmmFlaggerTask.predictionBed

        # HMM-Flagger summary files
        File predictionSummaryTsv = hmmFlaggerTask.predictionSummaryTsv
        File loglikelihoodTsv = hmmFlaggerTask.loglikelihoodTsv
        File miscFlaggerFilesTarGz = hmmFlaggerTask.miscFilesTarGz

        # Get projected bed files if there is any
        File? projectionSexBed = hmmFlaggerTask.projectionSexBed
        File? projectionSDBed = hmmFlaggerTask.projectionSDBed
        File? projectionCntrBed = hmmFlaggerTask.projectionCntrBed
        Array[File]? projectionAnnotationsBedArray = hmmFlaggerTask.projectionAnnotationsBedArray
        Array[File]? projectionBiasAnnotationsBedArray = hmmFlaggerTask.projectionBiasAnnotationsBedArray

        # bigwig files if user enabled outputting them
        Array[File]? bigwigArray = hmmFlaggerTask.bigwigArray

        # secphase output
        File? secphaseOutputLog = alignmentTask.secphaseOutputLog
        File? secphaseModifiedReadBlocksMarkersBed = alignmentTask.secphaseModifiedReadBlocksMarkersBed
        File? secphaseMarkerBlocksBed = alignmentTask.secphaseMarkerBlocksBed

    }
}
