version 1.0

import "../tasks/coverage/make_summary_table.wdl" as make_summary_table_t
import "../../ext/secphase/wdls/workflows/correct_bam.wdl" as correct_bam_t
import "../tasks/coverage/bam_coverage.wdl" as cov_t
import "../tasks/other/misc.wdl" as misc_t
import "../tasks/hmm_flagger/hmm_flagger.wdl" as hmm_flagger_t
import "../tasks/other/project_blocks_for_flagger.wdl" as project_t
import "../../ext/secphase/wdls/workflows/secphase.wdl" as secphase_t
import "../tasks/alignment/asm2asm_aligner.wdl" as asm2asm_t 
import "../tasks/alignment/produce_fai.wdl" as fai_t
import "../tasks/coverage/cov2wig.wdl" as cov2wig_t
import "../tasks/coverage/augment_coverage_by_labels.wdl" as augment_cov_t

workflow HMMFlaggerEndToEnd{
    meta {
        author: "Mobin Asri"
        email: "masri@ucsc.edu"
        description: "Running HMM-Flagger for evaluating a diploid assembly using long read alignments. More information at https://github.com/mobinasri/flagger"
    }
    parameter_meta {
        sampleName: "(Required) Sample name for example 'HG002'"
        suffix: "(Required) Suffix string that contains information about this analysis for example 'hifi_winnowmap_flagger_for_hprc'"
        hap1AssemblyFasta: "(Required) Path to uncompressed or gzip-compressed fasta file of the 1st haplotype."
        hap2AssemblyFasta: "(Required) Path to uncompressed or gzip-compressed fasta file of the 2nd haplotype."
        readAlignmentBam: "(Required) Path to sorted read alignment bam."
        readAlignmentBai: "(Required) Path to bam index for read alignment"
        presetForFlagger :  "(Required) HMM-Flagger preset (can be one of 'hifi', 'ont-r9', or 'ont-r10')"
        alphaTsv : "(Optional) The dependency factors for adjusting emission parameters with previous emission. This parameter is a tsv file with 4 rows and 4 columns with no header line. All numbers should be between 0 and 1. (Default = will be determined based on presetForFlagger)"
        maxReadDivergence: "Alignments with gap-compressed ratio higher than this will be filtered in the pre-process step. (Default: 0.1)"
        minReadLength: "Minimum read size (Default: 5000)"
        minAlignmentLength: "Minimum alignment size (Default: 5000)"
        downSamplingRate: "Rate of downsampling (Default: 1.0 which means no down-sampling)"
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
        enableOutputtingBigWig: "If true it will make bigwig files from cov files and output them. bigwig files can be easily imported into IGV sessions (Default: true)"
        enableCreatingConservativeBed: "If true it will map assembly contigs to themselves to create self-homology mappings and it will use them for filtering HMM-Flagger calls. Among outputs there will be a conservative bed file and also its related summary tables. (Default: true)"
        includeContigListText : "(Optional) Create coverage file and run HMM-Flagger only on these contigs (listed in a text file with one contig name per line). (Default: all contigs)"
        binArrayTsv : "(Optional)  A tsv file (tab-delimited) that contains bin arrays for stratifying results by event size. Bin intervals can have overlap. It should contain three columns. 1st column is the closed start of the bin and the 2nd column is the open end. The 3rd column has a name for each bin. (Default: all sizes in a single bin named ALL_SIZES)"
        chunkLen : "(Optional) The length of chunks for running HMM-Flagger. Each chunk will be processed in a separate thread before merging results together. (Default: 20000000)"
        windowLen : "(Optional) The length of windows for running HMM-Flagger. The coverage values will be averaged over the bases in each window and then the average value will be considered as an emission. (Default: will be determined based on presetForFlagger)"
        labelNames : "The names of the labels/states for reporting in the final summary tsv files (Default: 'Err,Dup,Hap,Col')"
        trackName : "The track name in the final BED file (Default: hmm_flagger)"
        numberOfIterations : "Number of EM iterations for estimating HMM parameters (Default:100)"
        convergenceTolerance : "Convergence tolerance. The EM iteration will stop once the difference between all model parameter values in two consecutive iterations is less than this value. (Default = 0.001)"
        maxHighMapqRatio : "Maximum ratio of high mapq coverage for duplicated state (Default = 0.25)"
        minHighMapqRatio : "Minimum ratio of high mapq coverage for collapsed state (Default = 0.75)"
        flaggerMoreOptions : "(Optional) More options for HMM-Flagger provided in a single string (Default = '')"
        modelType : "(Optional) Model type can be either 'gaussian', 'negative_binomial', or 'trunc_exp_gaussian' (Default = will be determined based on presetForFlagger)"
        flaggerMinimumBlockLenArray : "Array of minimum lengths for converting short non-Hap blocks into Hap blocks. Given numbers should be related to the states Err, Dup and Col respectively. (Default: [0,0,0])"
        flaggerMemSize : "Memory size in GB for running HMM-Flagger (Default : 32)"
        flaggerThreadCount : "Number of threads for running HMM-Flagger (Default : 16)"
        flaggerDockerImage : "Docker image for HMM-Flagger (Default : mobinasri/flagger:v1.2.0)"
        truthBedForMisassemblies : "(Optional) A BED file containing the coordinates and labels of the truth misassemblies. It can be useful when the misassemblies are simulated (e.g. with Falsifier) (Default: None)"
    }
    input{
        String sampleName
        String suffix

        File hap1AssemblyFasta
        File hap2AssemblyFasta
        File readAlignmentBam
        File readAlignmentBai
        File? alphaTsv
        Float maxReadDivergence = 0.1
        Int minReadLength = 5000
        Int minAlignmentLength = 5000
        Float downSamplingRate = 1.0
              
        File? includeContigListText
        File? binArrayTsv
        String labelNames = "Err,Dup,Hap,Col"
        String trackName = "hmm_flagger"
        Int numberOfIterations = 100
        Float convergenceTolerance = 0.001
        Float maxHighMapqRatio=0.25
        Float minHighMapqRatio=0.75
        String presetForFlagger
        String? flaggerMoreOptions
        String? modelType
        Int? chunkLen
        Int? windowLen
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

    # Index hap1 assembly
    call fai_t.produceFai as produceFaiHap1 {
        input:
            fasta = hap1AssemblyFasta
    }

    # Index hap2 assembly
    call fai_t.produceFai as produceFaiHap2 {
        input:
            fasta = hap2AssemblyFasta
    }

    # Index diploid assembly
    call fai_t.produceFai {
        input:
            fasta = createDipAsm.diploidAssemblyFastaGz
    }

    # if creating conservative bed is enabled by user
    # then map assembly to itself with -D parameter. 
    # This mapping will be used for filtering hmm-flagger 
    # predictions in a later task
    if (enableCreatingConservativeBed) {
        call asm2asm_t.asm2asmAlignment as selfHomologyMapping {
            input :
                aligner = "minimap2",
                preset = "asm5",
                options = "-D -I8g",
                queryAssemblyFasta = createDipAsm.diploidAssemblyFastaGz,
                refAssemblyFasta = createDipAsm.diploidAssemblyFastaGz
        }
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
                version = secphaseVersion,
        }
    }

    # Preprocess read alignment file by applying
    # the changes suggested by secphase and
    # filtering reads/alignments based on length
    # and divergence rate
    if (enableRunningSecphase || (0 < minReadLength) || (0 < minAlignmentLength) || (maxReadDivergence < 1.0)){
        call correct_bam_t.correctBam {
            input:
                bam = readAlignmentBam,
                phasingLogText = secphase.outLog,
                suffix = "corrected",
                options = "--primaryOnly --minReadLen ${minReadLength} --minAlignment ${minAlignmentLength} --maxDiv ${maxReadDivergence}",
                flagRemoveSupplementary = false,
                flagRemoveMultiplePrimary = false,
                dockerImage = secphaseDockerImage,
        }
    }


    Boolean needsAtLeastOneProjection = (0 < length(annotationsBedArrayToBeProjected)) || (0 < length(biasAnnotationsBedArrayToBeProjected) || defined(sexBedToBeProjected) || defined(SDBedToBeProjected) || defined(cntrBedToBeProjected))

    # Map each haplotype to the given reference
    # and then project bed files to the assembly 
    # coordinates. This projection will be performed
    # only for the annotations provided in the parameters
    # ending with ToBeProjected
    #
    #   1. annotationsBedArrayToBeProjected
    #   2. biasAnnotationsBedArrayToBeProjected
    #   3. sexBedToBeProjected
    #   4. SDBedToBeProjected
    #   5. cntrBedToBeProjected
    #   6. cntrCtBedToBeProjected (can be ignored since it is only for better projection of cntr blocks)
    #
    # if both arrays are empty then no projection will
    # be executed
    if (needsAtLeastOneProjection && defined(projectionReferenceFasta)) {
        call asm2asm_t.asm2asmAlignment as hap1ToRef {
            input :
                aligner="minimap2",
                preset="asm5",
                queryAssemblyFasta=hap1AssemblyFasta,
                refAssemblyFasta=select_first([projectionReferenceFasta]),
        }
        call asm2asm_t.asm2asmAlignment as hap2ToRef {
            input :
                aligner="minimap2",
                preset="asm5",
                queryAssemblyFasta=hap2AssemblyFasta,
                refAssemblyFasta=select_first([projectionReferenceFasta]),
        }
        call project_t.runProjectBlocksForFlagger as project{
            input:
                sampleName = sampleName,
                hap1AssemblyBam = hap1ToRef.sortedBamFile,
                hap2AssemblyBam = hap2ToRef.sortedBamFile, 
                refBiasedBlocksBedArray = biasAnnotationsBedArrayToBeProjected,
                additionalBedArray = annotationsBedArrayToBeProjected,
                refSexBed = sexBedToBeProjected,
                refSDBed = SDBedToBeProjected,
                refCntrBed = cntrBedToBeProjected,
                refCntrCtBed = cntrCtBedToBeProjected,
        }
    }

    if (enableDecomposingCntrBed && defined(cntrBed)){
        call decomposeCntrBed{
            input:
                cntrBed = select_first([cntrBed])
        }
    }

    # Get bed files containing potentially biased 
    # regions in asm coordinates, which can be 
    # projections from reference
    Array[File] potentialBiasesBedArrayInAsmCoor = flatten([select_first([decomposeCntrBed.cntrDecomposedBedFiles, []]), select_first([project.projectionBiasedBedArray, []]), select_first([biasAnnotationsBedArray, []])])
    
    # Get bed files containing additional stratifications
    # regions in asm coordinates, which can be
    # projections from reference
    Array[File] additionalStratificationBedArrayInAsmCoor = flatten([select_first([project.projectionAdditionalBedArray, []]), select_first([annotationsBedArray, []])])
     
    # Get bed files for SD, centromere and 
    # sex which can be projections from reference

    if (defined(project.projectionSDBed) || defined(SDBed)){
        File SDBedInAsmCoor = select_first([project.projectionSDBed, SDBed])
    }

    if (defined(project.projectionCntrBed) || defined(cntrBed)){
        File cntrBedInAsmCoor = select_first([decomposeCntrBed.cntrNoCtBed, project.projectionCntrBed, cntrBed])
    }

    if (defined(project.projectionSexBed) || defined(sexBed)){
        File sexBedInAsmCoor = select_first([project.projectionSexBed, sexBed])
    }

    call collectAnnotations{
        input:
            fasta = createDipAsm.diploidAssemblyFastaGz,
            minContigLen = 1000000,
            minContigLenName = "1Mb",
            includeContigListText = includeContigListText,
            biasAnnotationsBedArray = potentialBiasesBedArrayInAsmCoor,
            otherAnnotationsBedArray = additionalStratificationBedArrayInAsmCoor,
            difficultBed_1 = cntrBedInAsmCoor,
            difficultString_1 = "Cntr",
            difficultBed_2 = SDBedInAsmCoor,
            difficultString_2 = "SD",
            sexBed = sexBedInAsmCoor,
            dockerImage = flaggerDockerImage
    }

    # convert bam to cov and add annotation indices to the cov file
    call cov_t.bam2cov {
        input:
            bam = select_first([correctBam.correctedBam, readAlignmentBam]),
            bai = select_first([correctBam.correctedBamIndex, readAlignmentBai]),
            fasta = createDipAsm.diploidAssemblyFastaGz,
            suffix = "",
            mapqThreshold = 10,
            clipRatioThreshold = 0.1,
            downsampleRate = downSamplingRate,
            annotationBedArray = collectAnnotations.annotationBedArray,
            biasAnnotationNameArray = collectAnnotations.biasAnnotationNameArray,
            baselineAnnotationName = "WHOLE_GENOME_DEFAULT",
            includeContigListText = includeContigListText,
            runBiasDetection = (length(collectAnnotations.biasAnnotationNameArray) > 0),
            format = "all",
            memSize = 32,
            threadCount = 8,
            dockerImage = flaggerDockerImage,
    }


    call hmm_flagger_t.hmmFlagger {
        input:
            coverage = bam2cov.coverageGz,
            preset = presetForFlagger,
            suffix = suffix,
            binArrayTsv = binArrayTsv,
            labelNames = labelNames,
            trackName = trackName,
            numberOfIterations = numberOfIterations,
            convergenceTolerance = convergenceTolerance,
            maxHighMapqRatio = maxHighMapqRatio,
            minHighMapqRatio = minHighMapqRatio,
            moreOptions = flaggerMoreOptions,
            minimumBlockLenArray = flaggerMinimumBlockLenArray, 
            alphaTsv = alphaTsv,
            modelType = modelType,
            chunkLen = chunkLen,
            windowLen = windowLen,
            memSize = flaggerMemSize,
            threadCount = flaggerThreadCount,
            dockerImage = flaggerDockerImage,
    }


    # Get coordinates of canonical bases only (no "N" which may come from scaffolding)
    call misc_t.getCanonicalBasesBed as dipCanonical{
        input:
            assemblyFasta = createDipAsm.diploidAssemblyFastaGz,
            dockerImage = flaggerDockerImage
    }


    # make a truth bed file that contains
    # state indices instead of the names of the states
    if (defined(truthBedForMisassemblies)){
        call misc_t.getIndexLabeledBed as labelTruth{
            input:
                bed = select_first([truthBedForMisassemblies]),
                canonicalBasesDiploidBed = dipCanonical.canonicalBasesBed,
                addHapCoordinates = true
        }
    }

    # make a prediction bed file that contains
    # state indices instead of the names of the states
    call misc_t.getIndexLabeledBed as labelPrediction{
        input:
            bed = hmmFlagger.predictionBed
    }

    # Augment coverage file with prediction labels
    # it will add truth labels if available
    call augment_cov_t.augmentCoverageByLabels {
        input:
            coverage = bam2cov.coverageGz,
            fai = produceFai.fai,
            numberOfLabels = 4,
            truthBed = labelTruth.labeledBed,
            predictionBed = labelPrediction.labeledBed,
            includeContigListText = includeContigListText,
            suffix="augmented",
            dockerImage = flaggerDockerImage,
    }

    call make_summary_table_t.makeSummaryTable {
        input:
           coverage = augmentCoverageByLabels.augmentedCoverageGz,
           binArrayTsv = binArrayTsv,
           dockerImage = flaggerDockerImage
    }


    # if creating conservative bed was enabled by user
    # then 
    # - Filter hmm-flagger calls using self-homology mappings
    # - Make a coverage file augmented with conservative calls
    # - Make a summary table using conservative bed
    if (enableCreatingConservativeBed) {

        # filter calls with self-homology mappings
        call hmm_flagger_t.filterHmmFlaggerCalls as filterCalls {
            input:
                selfAsmMapBam = select_first([selfHomologyMapping.sortedBamFile]),
                flaggerBed = hmmFlagger.predictionBed
        }
        
        # make a prediction bed file that contains
        # state indices instead of the names of the states
        # for conservative calls
        call misc_t.getIndexLabeledBed as labelPredictionConservative{
            input:
                bed = filterCalls.conservativeBed
        }

        # Augment coverage file with prediction labels
        # it will add truth labels if available
        # for conservative calls
        call augment_cov_t.augmentCoverageByLabels as augmentCoverageByLabelsConservative{
            input:
                coverage = bam2cov.coverageGz,
                fai = produceFai.fai,
                numberOfLabels = 4,
                truthBed = labelTruth.labeledBed,
                predictionBed = labelPredictionConservative.labeledBed,
                includeContigListText = includeContigListText,
                suffix="augmented",
                dockerImage = flaggerDockerImage,
         }
         # make a summary table for conservative calls
         call make_summary_table_t.makeSummaryTable as makeSummaryTableConservative{
             input:
                 coverage = augmentCoverageByLabelsConservative.augmentedCoverageGz,
                 binArrayTsv = binArrayTsv,
                 dockerImage = flaggerDockerImage
         }
         # add labels for N bases and remove Hap labels from bed file
         # for conservative
         call getFinalBed as getFinalBedConservative {
            input:
                predictionBed = filterCalls.conservativeBed,
                canonicalBasesDiploidBed = dipCanonical.canonicalBasesBed,
                trackName = trackName,
                hap1Fai = produceFaiHap1.fai,
                hap2Fai = produceFaiHap2.fai,
                dockerImage = flaggerDockerImage
         }
    }


    # make bigwig files from cov files
    # bigwig files can be easily imported into IGV sessions
    if (enableOutputtingBigWig) {
        call cov2wig_t.cov2bigwig as cov2bigwig{
            input:
                coverage = bam2cov.coverageGz,
                windowLen = 1000,
                trackName = trackName,
                fai = produceFai.fai,
                dockerImage = flaggerDockerImage,
        }
    }

    # add labels for N bases and remove Hap labels from bed file
    call getFinalBed {
        input:
            predictionBed = hmmFlagger.predictionBed,
            canonicalBasesDiploidBed = dipCanonical.canonicalBasesBed, 
            trackName = trackName,
            hap1Fai = produceFaiHap1.fai,
            hap2Fai = produceFaiHap2.fai,
            dockerImage = flaggerDockerImage
    }

    if (defined(truthBedForMisassemblies)){
        File benchmarkingSummaryTsvOutput = makeSummaryTable.benchmarkingSummaryTsv
        File contiguitySummaryTsvOutput = makeSummaryTable.contiguitySummaryTsv
        if (enableCreatingConservativeBed) {
            File benchmarkingSummaryTsvOutputConservative = select_first([makeSummaryTableConservative.benchmarkingSummaryTsv])
            File contiguitySummaryTsvOutputConservative = select_first([makeSummaryTableConservative.contiguitySummaryTsv])
        }
    }

    output {
        File coverageGz = augmentCoverageByLabels.augmentedCoverageGz
        File biasTableTsv = bam2cov.biasTableTsv

        File? benchmarkingSummaryTsv = benchmarkingSummaryTsvOutput
        File? contiguitySummaryTsv = contiguitySummaryTsvOutput
        File fullStatsTsv = makeSummaryTable.fullStatsTsv

        File finalBed = getFinalBed.finalBed
        File finalBedHap1 = getFinalBed.finalBedHap1
        File finalBedHap2 = getFinalBed.finalBedHap2
        File predictionBed = hmmFlagger.predictionBed
        File loglikelihoodTsv = hmmFlagger.loglikelihoodTsv
        File miscFilesTarGz = hmmFlagger.outputTarGz

        # outputs for conservative calls (they will exist only if enableCreatingConservativeBed is true)
        File? benchmarkingSummaryTsvConservative = benchmarkingSummaryTsvOutputConservative
        File? contiguitySummaryTsvConservative = contiguitySummaryTsvOutputConservative
        File? fullStatsTsvConservative = makeSummaryTableConservative.fullStatsTsv
        
        File? finalBedConservative = getFinalBedConservative.finalBed
        File? finalBedConservativeHap1 = getFinalBedConservative.finalBedHap1
        File? finalBedConservativeHap2 = getFinalBedConservative.finalBedHap2
        File? predictionBedConservative = filterCalls.conservativeBed

        # get projected bed files if there is any
        File? projectionSexBed = project.projectionSexBed
        File? projectionSDBed = project.projectionSDBed
        File? projectionCntrBed = project.projectionCntrBed
        Array[File]? projectionAnnotationsBedArray = project.projectionAdditionalBedArray
        Array[File]? projectionBiasAnnotationsBedArray = project.projectionBiasedBedArray

        # bigwig files if user enabled outputting them
        Array[File]? bigwigArray = cov2bigwig.bigwigArray
        
        # secphase output
        File? secphaseOutputLog = secphase.outLog
        File? secphaseModifiedReadBlocksMarkersBed = secphase.modifiedReadBlocksMarkersBed
        File? secphaseMarkerBlocksBed = secphase.markerBlocksBed
    }
}

task decomposeCntrBed {
    input {
        File cntrBed
        Array[String] patterns = ["hsat1A", "hsat1B", "hsat2", "hsat3", "active_hor", "bsat"]
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=32
        String dockerImage="mobinasri/flagger:v1.2.0"
        Int preemptible=2
    }
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        mkdir output
        
        for pattern in ~{sep=" " patterns};
        do
            cat ~{cntrBed} | \
                grep -i ${pattern} | \
                awk '{print $1"\t"$2"\t"$3}' | \
                bedtools sort -i - | \
                bedtools merge -i - > output/censat_decomposed.${pattern}.bed
        done

        cat ~{cntrBed} | \
            awk '$4 != "ct"{print $1"\t"$2"\t"$3}' | \
            bedtools sort -i - | \
            bedtools merge -i - > output/censat.no_ct.bed
        
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        Array[File] cntrDecomposedBedFiles = glob("output/censat_decomposed.*")
        File cntrNoCtBed = glob("output/censat.no_ct.bed")[0]
    }
}



task getFinalBed {
    input {
        File predictionBed
        File canonicalBasesDiploidBed
        String trackName
        File hap1Fai
        File hap2Fai
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=32
        String dockerImage="mobinasri/flagger:v0.4.0"
        Int preemptible=2
    }
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        IN_BED_PATH="~{predictionBed}"
        PREFIX=$(basename ${IN_BED_PATH%.bed})

        mkdir output
        
        # make a bed file for Ns only
        # the color is black with an rgb of "0,0,0"
        # the label is "NNN"
        bedtools subtract \
            -a ~{predictionBed} \
            -b ~{canonicalBasesDiploidBed} | \
            cut -f1-3 | \
            bedtools merge -i - | \
            awk '{print $1"\t"$2"\t"$3"\tNNN\t0\t.\t"$2"\t"$3"\t0,0,0"}' > non_canonical.bed

        # make a BED with no Hap tracks
        cat ~{predictionBed} | \
            grep -v "Hap" | \
            grep -v "^track" | \
            bedtools intersect -a - -b ~{canonicalBasesDiploidBed}  > canonical.no_Hap.bed
        
        
        # add track name
        echo "track name=\"~{trackName}\" visibility=2 itemRgb=\"On\"" > output/${PREFIX}.canonical.no_Hap.bed

        # merge canonical and non-canonical tracks in the final bed
        cat non_canonical.bed canonical.no_Hap.bed | bedtools sort -i - >> output/${PREFIX}.canonical.no_Hap.bed

        # add track name for hap1
        echo "track name=\"~{trackName}_hap1\" visibility=2 itemRgb=\"On\"" > output/${PREFIX}.canonical.no_Hap.hap1.bed
        cut -f1 ~{hap1Fai} | grep -F -f - output/${PREFIX}.canonical.no_Hap.bed >> output/${PREFIX}.canonical.no_Hap.hap1.bed

        # add track name for hap2
        echo "track name=\"~{trackName}_hap2\" visibility=2 itemRgb=\"On\"" > output/${PREFIX}.canonical.no_Hap.hap2.bed
        cut -f1 ~{hap2Fai} | grep -F -f - output/${PREFIX}.canonical.no_Hap.bed >> output/${PREFIX}.canonical.no_Hap.hap2.bed

    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File finalBed = glob("output/*.no_Hap.bed")[0]
        File finalBedHap1 = glob("output/*.no_Hap.hap1.bed")[0]
        File finalBedHap2 = glob("output/*.no_Hap.hap2.bed")[0]
    }
}


task collectAnnotations{
    input{
        File fasta
        Int minContigLen = 1000000
        String minContigLenName = "1Mb"
        File? includeContigListText
        Array[File] biasAnnotationsBedArray = []
        Array[File] otherAnnotationsBedArray = []
        File? difficultBed_1
        File? difficultBed_2
        String? difficultString_1
        String? difficultString_2  
        File? sexBed
        # runtime configurations
        Int memSize=8
        Int threadCount=4
        Int diskSize=32
        String dockerImage="mobinasri/flagger:v1.2.0"
        Int preemptible=2
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace
        

        FA_PREFIX=$(echo $(basename ~{fasta}) | sed -e 's/\.fa$//' -e 's/\.fa.gz$//' -e 's/\.fasta$//' -e 's/\.fasta.gz$//')

        if [[ "~{fasta}" == *.gz ]]; then
            gunzip -c ~{fasta} > ${FA_PREFIX}.fa
        else
            cp ~{fasta} ${FA_PREFIX}.fa
        fi

        samtools faidx ${FA_PREFIX}.fa

        mkdir -p intermediate
        python3 /home/scripts/get_N_coords.py --inputFasta ${FA_PREFIX}.fa > intermediate/N_coords.bed

        if [ -n "~{includeContigListText}" ]
        then
            cat ${FA_PREFIX}.fa.fai | \
                awk '{print $1"\t0\t"$2}' | \
                bedtools sort -i - | \
                grep -F -f ~{includeContigListText} | \
                bedtools subtract -a - -b intermediate/N_coords.bed > whole_genome.bed
        else
            cat ${FA_PREFIX}.fa.fai | \
                awk '{print $1"\t0\t"$2}' | \
                bedtools sort -i - | \
                bedtools subtract -a - -b intermediate/N_coords.bed > whole_genome.bed
        fi

        cat whole_genome.bed | awk '($3-$2)>~{minContigLen}{print $0}' > whole_genome.gt_~{minContigLenName}.bed

        if [ -n "~{difficultBed_1}" ]
        then
            bedtools intersect -a ~{difficultBed_1} -b whole_genome.bed > ~{difficultString_1}.bed
            bedtools intersect -a ~{difficultString_1}.bed -b whole_genome.gt_~{minContigLenName}.bed > ~{difficultString_1}.gt_~{minContigLenName}.bed
            bedtools subtract -a whole_genome.bed -b ~{difficultString_1}.bed > Non_~{difficultString_1}.bed
        fi

        if [ -n "~{difficultBed_2}" ]
        then
            bedtools intersect -a ~{difficultBed_2} -b whole_genome.bed > ~{difficultString_2}.bed
            bedtools intersect -a ~{difficultString_2}.bed -b whole_genome.gt_~{minContigLenName}.bed > ~{difficultString_2}.gt_~{minContigLenName}.bed
            bedtools subtract -a whole_genome.bed -b ~{difficultString_2}.bed > Non_~{difficultString_2}.bed
        fi

        if [[ -n "~{difficultBed_2}" && -n "~{difficultBed_1}" ]]
        then
            cat ~{difficultString_1}.bed ~{difficultString_2}.bed | bedtools sort -i - | bedtools merge -i - > ~{difficultString_1}_or_~{difficultString_2}.bed   
            bedtools intersect -a Non_~{difficultString_1}.bed -b Non_~{difficultString_2}.bed > Non_~{difficultString_1}_and_Non_~{difficultString_2}.bed
            bedtools intersect -a Non_~{difficultString_1}_and_Non_~{difficultString_2}.bed -b whole_genome.gt_~{minContigLenName}.bed > Non_~{difficultString_1}_and_Non_~{difficultString_2}.gt_~{minContigLenName}.bed
        fi

        if [ -n "~{sexBed}" ]
        then
            bedtools intersect -a ~{sexBed} -b whole_genome.bed > sex.bed
            bedtools subtract -a whole_genome.bed -b sex.bed > autosome.bed
        fi

        if [[ -n "~{sexBed}" && -n "~{difficultBed_1}" ]]
        then
            bedtools intersect -a autosome.bed -b Non_~{difficultString_1}.bed > autosome_Non_~{difficultString_1}.bed
        fi

        if [[ -n "~{sexBed}" && -n "~{difficultBed_2}" ]]
        then
            bedtools intersect -a autosome.bed -b Non_~{difficultString_2}.bed > autosome_Non_~{difficultString_2}.bed
        fi

        if [[ -n "~{sexBed}" && -n "~{difficultBed_1}" && -n "~{difficultBed_2}" ]]
        then
            bedtools intersect -a autosome.bed -b Non_~{difficultString_1}_and_Non_~{difficultString_2}.bed > autosome_Non_~{difficultString_1}_and_Non_~{difficultString_2}.bed
            bedtools intersect -a autosome_Non_~{difficultString_1}_and_Non_~{difficultString_2}.bed \
                           -b whole_genome.gt_~{minContigLenName}.bed > autosome_Non_~{difficultString_1}_and_Non_~{difficultString_2}.gt_~{minContigLenName}.bed
        fi

        touch biased_annotation_names.txt

        if [ ~{length(biasAnnotationsBedArray)} -gt "0" ]
        then
            for BED_FILE in ~{sep=" " biasAnnotationsBedArray}
            do
                bedtools intersect -a ${BED_FILE} -b whole_genome.bed > $(basename ${BED_FILE})
                echo $(basename ${BED_FILE%%.bed}) >> biased_annotation_names.txt
            done
        fi

         
        if [ ~{length(otherAnnotationsBedArray)} -gt "0" ]
        then
            for BED_FILE in ~{sep=" " otherAnnotationsBedArray}
            do
                bedtools intersect -a ${BED_FILE} -b whole_genome.bed > $(basename ${BED_FILE})
            done
        fi

    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output{
        Array[File] annotationBedArray = glob("*.bed")
        Array[String] biasAnnotationNameArray = read_lines("biased_annotation_names.txt")
    }
}

