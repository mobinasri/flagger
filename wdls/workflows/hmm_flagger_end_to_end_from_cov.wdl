version 1.0

import "../tasks/coverage/make_summary_table.wdl" as make_summary_table_t
import "../tasks/other/misc.wdl" as misc_t
import "../tasks/hmm_flagger/hmm_flagger.wdl" as hmm_flagger_t
import "../tasks/alignment/asm2asm_aligner.wdl" as asm2asm_t
import "../tasks/alignment/produce_fai.wdl" as fai_t
import "../tasks/coverage/cov2wig.wdl" as cov2wig_t
import "../tasks/coverage/augment_coverage_by_labels.wdl" as augment_cov_t
import "hmm_flagger_end_to_end.wdl" as hmm_flagger_end_to_end_t

workflow HMMFlaggerEndToEndFromCov{
    meta {
        author: "Mobin Asri"
        email: "masri@ucsc.edu"
        description: "Running HMM-Flagger for evaluating a diploid assembly starting from an existing coverage file (cov.gz) instead of a read alignment bam. The coverage file is expected to be already annotated (and bias-detected) e.g. the 'coverageGz' output of HMMFlaggerEndToEnd. Any truth/prediction label columns in the input coverage file are removed before running HMM-Flagger. More information at https://github.com/mobinasri/flagger"
    }
    parameter_meta {
        sampleName: "(Required) Sample name for example 'HG002'"
        suffix: "(Required) Suffix string that contains information about this analysis for example 'hifi_winnowmap_flagger_for_hprc'"
        hap1AssemblyFasta: "(Required) Path to uncompressed or gzip-compressed fasta file of the 1st haplotype."
        hap2AssemblyFasta: "(Required) Path to uncompressed or gzip-compressed fasta file of the 2nd haplotype."
        coverage: "(Required) Path to a coverage file (cov or cov.gz) created for the diploid assembly (e.g. output of bam2cov or the 'coverageGz' output of HMMFlaggerEndToEnd). Annotation and region (bias) information in its header will be reused. Existing truth/prediction label columns will be removed."
        presetForFlagger :  "(Required) HMM-Flagger preset (can be one of 'hifi', 'ont-r9', or 'ont-r10')"
        alphaTsv : "(Optional) The dependency factors for adjusting emission parameters with previous emission. This parameter is a tsv file with 4 rows and 4 columns with no header line. All numbers should be between 0 and 1. (Default = will be determined based on presetForFlagger)"
        includeContigListText : "(Optional) Only used for filtering the truth bed file before augmenting the coverage file. The coverage file itself is not subset. (Default: all contigs)"
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
        flaggerDockerImage : "Docker image for HMM-Flagger (Default : mobinasri/flagger:v1.2.1)"
        enableOutputtingBigWig: "If true it will make bigwig files from cov files and output them. bigwig files can be easily imported into IGV sessions (Default: true)"
        enableCreatingConservativeBed: "If true it will map assembly contigs to themselves to create self-homology mappings and it will use them for filtering HMM-Flagger calls. Among outputs there will be a conservative bed file and also its related summary tables. (Default: true)"
        truthBedForMisassemblies : "(Optional) A BED file containing the coordinates and labels of the truth misassemblies. It can be useful when the misassemblies are simulated (e.g. with Falsifier) (Default: None)"
    }
    input{
        String sampleName
        String suffix

        File hap1AssemblyFasta
        File hap2AssemblyFasta
        File coverage
        File? alphaTsv

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
        String flaggerDockerImage="mobinasri/flagger:v1.2.1"

        Boolean enableOutputtingBigWig = true
        Boolean enableCreatingConservativeBed = true
        File? truthBedForMisassemblies
    }

    Int mapqThreshold = 10

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
                refAssemblyFasta = createDipAsm.diploidAssemblyFastaGz,
                memSize = 64
        }
    }

    # Remove truth/prediction columns (if any) from the input
    # coverage file so that it looks like a fresh bam2cov output
    call removeLabelsFromCoverage {
        input:
            coverage = coverage,
            dockerImage = flaggerDockerImage
    }

    call hmm_flagger_t.hmmFlagger {
        input:
            coverage = removeLabelsFromCoverage.coverageGz,
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
            coverage = removeLabelsFromCoverage.coverageGz,
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
                selfAsmMapBamIndex = select_first([selfHomologyMapping.sortedBamIndexFile]),
                flaggerBed = hmmFlagger.predictionBed,
                dockerImage = flaggerDockerImage
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
                coverage = removeLabelsFromCoverage.coverageGz,
                fai = produceFai.fai,
                numberOfLabels = 4,
                truthBed = labelTruth.labeledBed,
                predictionBed = labelPredictionConservative.labeledBed,
                includeContigListText = includeContigListText,
                suffix="augmented.conservative",
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
         call hmm_flagger_end_to_end_t.getFinalBed as getFinalBedConservative {
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
                coverage = removeLabelsFromCoverage.coverageGz,
                windowLen = 1000,
                minHighMapqCov = 4,
                minMapq = mapqThreshold,
                hap1Fai = produceFaiHap1.fai,
                hap2Fai = produceFaiHap2.fai,
                trackName = trackName,
                dockerImage = flaggerDockerImage,
        }
    }

    # add labels for N bases and remove Hap labels from bed file
    call hmm_flagger_end_to_end_t.getFinalBed as getFinalBed {
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

        # bigwig files if user enabled outputting them
        Array[File]? bigwigArray = cov2bigwig.bigwigArray
        File? mappableHap1Bed = cov2bigwig.mappableHap1Bed
        File? mappableHap2Bed = cov2bigwig.mappableHap2Bed
    }
}

task removeLabelsFromCoverage {
    input {
        File coverage
        # runtime configurations
        Int memSize=8
        Int threadCount=4
        Int diskSize=ceil(3 * size(coverage, "GB")) + 32
        String dockerImage="mobinasri/flagger:v1.2.1"
        Int preemptible=2
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        COV_FILE_PATH="~{coverage}"
        # strip ".augmented" (and ".conservative") so output names look like the original bam2cov output
        PREFIX=$(echo $(basename ${COV_FILE_PATH%.gz}) | sed -e 's/\.cov$//' -e 's/\.bed$//' -e 's/\.augmented\.conservative$//' -e 's/\.augmented$//')

        if [[ "${COV_FILE_PATH}" == *.gz ]]; then
            CAT_CMD="zcat"
        else
            CAT_CMD="cat"
        fi

        mkdir output
        # Columns of a cov file (non start-only mode):
        #   start, end, coverage, high_mapq_cov, high_clip_cov, annotations, region, [truth], [prediction]
        # Keep only the first 7 columns and reset truth/prediction tags in the header
        ${CAT_CMD} ${COV_FILE_PATH} | \
            awk 'BEGIN{OFS="\t"}
                 /^#label:len:/ {next}
                 /^#truth:/ {print "#truth:false"; next}
                 /^#prediction:/ {print "#prediction:false"; next}
                 /^#start-only:true/ {print "ERROR: start-only coverage files are not supported" > "/dev/stderr"; exit 1}
                 /^#/ || /^>/ {print; next}
                 {print $1,$2,$3,$4,$5,$6,$7}' | \
            pigz -p ~{threadCount} -c > output/${PREFIX}.cov.gz
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File coverageGz = glob("output/*.cov.gz")[0]
    }
}
