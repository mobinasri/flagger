version 1.0

import "../tasks/other/misc.wdl" as misc_t
import "hmm_flagger_end_to_end.wdl" as flagger_t
import "long_read_aligner_scattered.wdl" as aligner_t



workflow HMMFlaggerEndToEndWithMapping{
    meta {
        author: "Mobin Asri"
        email: "masri@ucsc.edu"
        description: "Mapping reads and running Flagger for evaluating a diploid assembly using long read alignments.More information at https://github.com/mobinasri/flagger"
    }
    parameter_meta {
        sampleName: "(Required) Sample name for example 'HG002'"
        suffixForFlagger: "(Required) Suffix string that contains information about this analysis for example 'hifi_winnowmap_flagger_for_hprc'"
        suffixForMapping: "(Required) Suffix string that contains information about this alignment. It will be appended to the name of the final alignment. For example 'hifi_winnowmap_v2.03_hprc_y2'"
        hap1AssemblyFasta: "(Required) Path to uncompressed or gzip-compressed fasta file of the 1st haplotype."
        hap2AssemblyFasta: "(Required) Path to uncompressed or gzip-compressed fasta file of the 2nd haplotype."
        readFiles: "(Required) Array of read files. Their format can be either fastq, fq, fastq.gz, fq.gz, bam or cram. For cram format referenceFastaForReadExtraction should also be passed."
        presetForFlagger :  "(Required) HMM-Flagger preset (can be one of 'hifi', 'ont-r9', or 'ont-r10')"
        aligner: "Name of the aligner. It can be either minimap2, winnowmap or veritymap. (Default = minimap2)"
        presetForMapping: "(Required) Paremeter preset should be selected based on aligner and sequencing platform. Common presets are lr:hqae/map-ont for minimap2, map-pb/map-ont for winnowmap and hifi-haploid/hifi-haploid-complete/hifi-diploid/ont-haploid-complete for veritymap"
        alphaTsv : "(Optional) The dependency factors for adjusting emission parameters with previous emission. This parameter is a tsv file with 4 rows and 4 columns with no header line. All numbers should be between 0 and 1. (Default = will be determined based on presetForFlagger)"
        kmerSize: "The kmer size for using minimap2 or winnowmap. With winnowmap kmer size should be 15 and with minimap2 kmer size should be 15 and 25 for using the presets map-ont and lr:hqae respectively."
        alignerOptions: "Aligner options. It can be something like '--eqx --cs -Y -L -y' for minimap2/winnowmap. Note that if assembly is diploid and aligner is either minimap2 or winnowmap '-I8g' is necessary. If the reads contain modification tags and these tags are supposed to be present in the final alignment file, alignerOptions should contain '-y' and the aligner should be either minimap2 or winnowmap. If running secphase is enabled it is recommended to add '-p0.5' to alignerOptions; it will keep more secondary alignments so secphase will have more candidates per read. For veritymap '--careful' can be used but not recommended for whole-genome assembly since it increases the runtime dramatically."
        readExtractionOptions: "The options to be used while converting bam to fastq with samtools fastq. If the reads contain epigentic modification tags it is necessary to use '-TMm,Ml,MM,ML'"
        referenceFastaForReadExtraction: "If reads are in CRAM format then the related reference should be passed to this parameter for read extraction."
        enableAddingMDTag: "If true it will call samtools calmd to add MD tags to the final bam file. [Default = true]"
        splitNumber: "The number of chunks which the input reads should be equally split into. Note that enableSplittingReadsEqually should be set to true if user wants to split reads into equally sized chunks. [Default = 16]"
        enableSplittingReadsEqually: "If true it will merge all reads together and then split them into multiple chunks of roughly equal size. Each chunk will then be aligned via a separate task. This feature is useful for running alignment on cloud/slurm systems where there are  multiple nodes available with enough computing power and having alignments tasks distributed among small nodes is more efficient or cheaper than running a single alignment task in a large node. If the  whole workflow is being on a single node it is not recommened to use this feature since mergin and splitting reads takes its own time. [Default = false]"
        minReadLengthForMapping: "If it is greater than zero, a task will be executed for filtering reads shorter than this value before alignment. [Default = 0]"
        alignerThreadCount : "The number of threads for mapping in each alignment task [Default = 16]"
        alignerMemSize : "The size of the memory in Gb for mapping in each alignment task [Default = 48]"
        alignerDockerImage : "The mapping docker image [Default = 'mobinasri/long_read_aligner:v1.1.0']"
        correctBamOptions : "Options for the correct_bam program that can filters short/highly divergent alignments [ Default = '--primaryOnly --minReadLen 5000 --minAlignment 5000 --maxDiv 0.1' ]"
        enableRunningCorrectBam: "If true it will run correct_bam with correctBamOptions"
        downSamplingRateForFlagger: "Rate of downsampling (Default: 1.0 which means no down-sampling)"
        sexBed: "(Optional) bed file containing regions assigned to X/Y chromosomes. (in asm coordinates)"
        sexBedToBeProjected: "(Optional) bed file containing regions assigned to X/Y chromosomes. (in ref coordinates)"
        SDBed: "(Optional) Bed file containing Segmental Duplications. (in asm coordinates)"
        SDBedToBeProjected: "(Optional) Bed file containing Segmental Duplications. (in ref coordinates)"
        cntrBed: "(Optional) Bed file containing peri/centromeric satellites (ASat, HSat, bSat, gSat) without 'ct' blocks. (in asm coordinates)"
        enableDecomposingCntrBed: "If true it means that cntrBed contains different satellite families like ASat and HSat. This Bed file can be decomposed into multiple bed files; one for each family. This decomposition will benefit bias detection since different types of repeat arrays might have different coverage biases. If this option is true biasAnnotationsBedArray will be filled automatically with the decomposed censat bed files and users can exclude it from their input json. This option shouldn't be true if the names of the satellite families are not mentioned in the tracks of cntrBed. [list of patterns for decomposing (case insensitive): 'hsat1A', 'hsat1B', 'hsat2', 'hsat3', 'active_hor', 'bsat'] (Default: false)"
        cntrBedToBeProjected: "(Optional) Bed file containing peri/centromeric satellites (ASat, HSat, bSat, gSat) without 'ct' blocks. (in ref coordinates)"
        cntrCtBed: "(Optional) Bed file containing centromere transition 'ct' blocks. (in asm coordinates)"
        cntrCtBedToBeProjected: "(Optional) Bed file containing centromere transition 'ct' blocks. (in ref coordinates)"
        annotationsBedArray: "(Optional) list of annotations to be used for augmenting coverage files and stratifying final HMM-Flagger results. These annotations should be in the coordinates of the assembly under evaluation."
        annotationsBedArrayToBeProjected: "(Optional) list of annotations to be used for augmenting coverage files and stratifying final HMM-Flagger results. These annotations are not in the coordinates of the assembly and they need to be projected from the given projection reference to the assembly coordinates."
        biasAnnotationsBedArray: "(Optional) list similar to annotationsBedArray but these annotations potentially have read coverage biases like HSat2. (in the assembly coordinates)"
        biasAnnotationsBedArrayToBeProjected: "(Optional) list similar to annotationsBedArrayToBeProjected but these annotations potentially have read coverage biases like HSat2. (in the reference coordinates)"
        projectionReferenceFasta: "(Optional) If any of the parameters ending with 'ToBeProjected' is not empty a reference fasta should be passed for performing the necessary projections (Default: '')"
        enableRunningSecphase : "If true it will run secphase in the marker mode using the parameters starting with 'secphase' otherwise skip it. (Default: false)"
        secphaseDockerImage: "Docker image for running Secphase (Default: mobinasri/secphase:v0.4.4)"
        secphaseOptions: "String containing secphase options (can be either --hifi or --ont). (Default --hifi)"
        secphaseVersion: "Secphase version. (Default: v0.4.4)"
        includeContigListText : "(Optional) Create coverage file and run HMM-Flagger only on these contigs (listed in a text file with one contig name per line). (Default: all contigs)"
        binSizeArrayTsv : "(Optional)  A tsv file (tab-delimited) that contains bin arrays for stratifying results by event size. Bin intervals can have overlap. It should contain three columns. 1st column is the closed start of the bin and the 2nd column is the open end. The 3rd column has a name for each bin. (Default: all sizes in a single bin named ALL_SIZES)"
        chunkLen : "(Optional) The length of chunks for running HMM-Flagger. Each chunk will be processed in a separate thread before merging results together. (Default: 20000000)"
        windowLen : "(Optional) The length of windows for running HMM-Flagger. The coverage values will be averaged over the bases in each window and then the average value will be considered as an emission. (Default: will be determined based on presetForFlagger)"
        labelNames : "The names of the labels/states for reporting in the final summary tsv files (Default: 'Err,Dup,Hap,Col')"
        trackName : "The track name in the final BED file (Default: hmm_flagger)"
        numberOfIterationsForFlagger : "Number of EM iterations for estimating HMM parameters (Default:100)"
        convergenceToleranceForFlagger : "Convergence tolerance. The EM iteration will stop once the difference between all model parameter values in two consecutive iterations is less than this value. (Default = 0.001)"
        maxHighMapqRatio : "Maximum ratio of high mapq coverage for duplicated state (Default = 0.25)"
        minHighMapqRatio : "Minimum ratio of high mapq coverage for collapsed state (Default = 0.75)"
        flaggerMoreOptions : "(Optional) More options for HMM-Flagger provided in a single string (Default = '')"
        modelType : "(Optional) Model type can be either 'gaussian', 'negative_binomial', or 'trunc_exp_gaussian' (Default = will be determined based on presetForFlagger)"
        flaggerMinimumBlockLenArray : "Array of minimum lengths for converting short non-Hap blocks into Hap blocks. Given numbers should be related to the states Err, Dup and Col respectively. (Default: [0,0,0])"
        flaggerMemSize : "Memory size in GB for running HMM-Flagger (Default : 32)"
        flaggerThreadCount : "Number of threads for running HMM-Flagger (Default : 16)"
        flaggerDockerImage : "Docker image for HMM-Flagger (Default : mobinasri/flagger:v1.2.0)"
        enableOutputtingBigWig: "If true it will make bigwig files from cov files and output them. bigwig files can be easily imported into IGV sessions (Default: true)"
        enableCreatingConservativeBed: "If true it will map assembly contigs to themselves to create self-homology mappings and those mappings will be used for filtering HMM-Flagger calls. Among outputs there will be a conservative bed file and also its related summary tables. (Default: true)"
        enableOutputtingBam: "If true it will make bigwig files from cov files and output them. bigwig files can be easily imported into IGV sessions (Default: false)" 
        truthBedForMisassemblies : "(Optional) A BED file containing the coordinates and labels of the truth misassemblies. It can be useful when the misassemblies are simulated (e.g. with Falsifier) (Default: None)"
    }

    input{
        String sampleName
        String suffixForMapping
        String suffixForFlagger
        File hap1AssemblyFasta
        File hap2AssemblyFasta
        Array[File] readFiles
        String presetForFlagger
        String presetForMapping
        File? alphaTsv

        String aligner="minimap2"
        Int kmerSize = 15
        String alignerOptions="--eqx --cs -Y -L -y -I8g"
        String readExtractionOptions="-TMM,ML,Mm,Ml"
        File? referenceFastaForReadExtraction
        Boolean enableAddingMDTag=true
        Int splitNumber = 16
        Boolean enableSplittingReadsEqually=false

        Int minReadLengthForMapping = 0
        String correctBamOptions = "--primaryOnly --minReadLen 5000 --minAlignment 5000 --maxDiv 0.1" 
        Boolean enableRunningCorrectBam=false
        Int alignerThreadCount = 16
        Int alignerMemSize = 48
        String alignerDockerImage = "mobinasri/long_read_aligner:v1.1.0"
        Float downSamplingRateForFlagger = 1.0

        File? includeContigListText
        File? binSizeArrayTsv
        String? modelType
        Int? chunkLen
        Int? windowLen
        String labelNames = "Err,Dup,Hap,Col"
        String trackName = "hmm_flagger"
        Int numberOfIterationsForFlagger = 100
        Float convergenceToleranceForFlagger = 0.001
        Float maxHighMapqRatio=0.25
        Float minHighMapqRatio=0.75
        String? flaggerMoreOptions
        Array[Int] flaggerMinimumBlockLenArray = []
        Int flaggerMemSize=32
        Int flaggerThreadCount=16
        String flaggerDockerImage="mobinasri/flagger:v1.2.0"

        File? sexBed
        File? SDBed
        File? cntrBed # censat annotation with no "ct"
        File? cntrCtBed
        Boolean enableDecomposingCntrBed = false

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
        String secphaseDockerImage = "mobinasri/secphase:v0.4.4"
        String secphaseOptions = "--hifi"
        String secphaseVersion = "v0.4.4"

        Boolean enableOutputtingBigWig = true
        Boolean enableCreatingConservativeBed = true
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
            enableRunningCorrectBam = enableRunningCorrectBam
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
            presetForFlagger = presetForFlagger,

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
            minHighMapqRatio = minHighMapqRatio,
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
            enableDecomposingCntrBed = enableDecomposingCntrBed,

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
            enableCreatingConservativeBed = enableCreatingConservativeBed,
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
        File finalPredictionBedHap1 = hmmFlaggerTask.finalBedHap1
        File finalPredictionBedHap2 = hmmFlaggerTask.finalBedHap2
        File intermediatePredictionBed = hmmFlaggerTask.predictionBed

        # HMM-Flagger summary files
        File loglikelihoodTsv = hmmFlaggerTask.loglikelihoodTsv
        File miscFlaggerFilesTarGz = hmmFlaggerTask.miscFilesTarGz

        File? benchmarkingSummaryTsv = hmmFlaggerTask.benchmarkingSummaryTsv
        File? contiguitySummaryTsv = hmmFlaggerTask.contiguitySummaryTsv
        File fullStatsTsv = hmmFlaggerTask.fullStatsTsv

        # outputs for conservative calls (they will exist only if enableCreatingConservativeBed is true)
        File? finalPredictionBedConservative = hmmFlaggerTask.finalBedConservative
        File? finalPredictionBedConservativeHap1 = hmmFlaggerTask.finalBedConservativeHap1
        File? finalPredictionBedConservativeHap2 = hmmFlaggerTask.finalBedConservativeHap2
        File? intermediatePredictionBedConservative = hmmFlaggerTask.predictionBedConservative

        File? benchmarkingSummaryTsvConservative = hmmFlaggerTask.benchmarkingSummaryTsvConservative
        File? contiguitySummaryTsvConservative = hmmFlaggerTask.contiguitySummaryTsvConservative
        File? fullStatsTsvConservative = hmmFlaggerTask.fullStatsTsvConservative

 
        # Get projected bed files if there is any
        File? projectionSexBed = hmmFlaggerTask.projectionSexBed
        File? projectionSDBed = hmmFlaggerTask.projectionSDBed
        File? projectionCntrBed = hmmFlaggerTask.projectionCntrBed
        Array[File]? projectionAnnotationsBedArray = hmmFlaggerTask.projectionAnnotationsBedArray
        Array[File]? projectionBiasAnnotationsBedArray = hmmFlaggerTask.projectionBiasAnnotationsBedArray

        # bigwig files if user enabled outputting them
        Array[File]? bigwigArray = hmmFlaggerTask.bigwigArray
        File? mappableHap1Bed = hmmFlaggerTask.mappableHap1Bed
        File? mappableHap2Bed = hmmFlaggerTask.mappableHap2Bed
        
        # read alignment bam file if user enabled outputting it
        File? readAlignmentBam = readAlignmentBamToOutput
        File? readAlignmentBai = readAlignmentBaiToOutput

        # secphase output
        File? secphaseOutputLog = alignmentTask.secphaseOutputLog
        File? secphaseModifiedReadBlocksMarkersBed = alignmentTask.secphaseModifiedReadBlocksMarkersBed
        File? secphaseMarkerBlocksBed = alignmentTask.secphaseMarkerBlocksBed

    }
}
