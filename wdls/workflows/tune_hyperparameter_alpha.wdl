version 1.0 

import "../tasks/other/tune_alpha.wdl" as tune_alpha_t
import "../tasks/coverage/make_summary_table.wdl" as make_summary_table_t
import "../tasks/coverage/bam_coverage.wdl" as cov_t
import "../tasks/coverage/cov2bin.wdl" as cov2bin_t
import "../tasks/other/misc.wdl" as misc_t
import "../tasks/hmm_flagger/hmm_flagger.wdl" as hmm_flagger_t
import "../tasks/alignment/produce_fai.wdl" as fai_t
import "../tasks/coverage/cov2wig.wdl" as cov2wig_t
import "../tasks/coverage/augment_coverage_by_labels.wdl" as augment_cov_t


workflow runTuneHyperparameterAlpha{
    meta {
        author: "Mobin Asri"
        email: "masri@ucsc.edu"
        description: "Tuning hyperparameter of alpha for HMM-Flagger. More information at https://github.com/mobinasri/flagger"
    }
    parameter_meta {
        bamArray: "(Required) An array of BAM files, each containing read alignments to falsified assemblies."
        baiArray: "(Required) An array of BAM index files, each corresponding to a BAM file in the bamArray, in the same order."
        faiArray: "(Required) An array of FASTA files, each corresponding to a BAM file in the bamArray, in the same order."
        misassembliesBedArray: "(Required) An array of BED files, each containing the coordinates of induced misassemblies, associated with the BAM files in the bamArray, in the same order."
        validationContigPatternArray: "(Required) An array of strings, where at least one of the strings should be present in the name of a contig for it to be considered in the validation data set. (Default: ['chr20', 'chr21', 'chr22'])"
        annotationBedArray2D: "(Optional) A 2D array of annotation bed files. The length of array should be equal to the number of bam files. Each element is an array of BED files, each containing the coordinates of an annotation e.g. HSat1A or HSat2 (Default: None)"
        biasAnnotationNameArray2D: "(Optional) A 2D array of the names of the annotations potentially having coverage biases. The length of array should be equal to the number of bam files. Each element is an array of names related to the annotations that should be considered in the bias detection step. The given names should match the names of the files in annotationBedArray2D without the suffix '.bed'. The length of each internal array cannot be larger than its equavalent array in annotationBedArray2D. (Default: None)"
        downsampleRateArray: "An array of down-sampling rates (Default: [1.0])"
        flaggerDockerImage: "Docker image for running HMM-Flagger and related programs (Default: 'mobinasri/flagger:v1.0.0_alpha--73711ad369ab956a33be9d14edf285a2191c2a94')"
        threadsPerFlaggerRun: "The number of threads for HMM-Flagger (Default: 8)"
        maxFlaggerRunsPerIteration: "The maximum number of HMM-Flagger runs per EGO iteration (Default: 8)"
        windowLen: "The size of the window for HMM-Flagger (Default: 4000)"
        chunkLen: "The size of the chunk for HMM-Flagger (Default: 20000000)"
        modelType: "Model type for HMM-Flagger (Default: 'gaussian')"
        otherFlaggerParamsForTuning: "The remaining parameter values for HMM-Flagger (Default: '--convergenceTol 1e-2  --iterations 10')"
        binArrayTsv: "(Optional)  A tsv file (tab-delimited) that contains bin arrays for stratifying results by event size. Bin intervals can have overlap. It should contain three columns. 1st column is the closed start of the bin and the 2nd column is the open end. The 3rd column has a name for each bin. (Default: all sizes in a single bin named ALL_SIZES)"
        numberOfEgoStartPoints: "Number of random alpha matrices to be evaluated before performing EGO algorithm (Default: 10)"
        numberOfEgoIterations: "Number of EGO iterations (Default: 50)"
        runFlaggerPostOptForTrain: "If this parameter is true after finding the optimal alpha matrix run HMM-Flagger on all training coverage files with more EM iterations and then make summary tables with the resolution of bases not windows (HMM-Flagger outputs summary tables with the resolution of windows) (Default = false)"
        runFlaggerPostOptForValidation: "If this parameter is true after finding the optimal alpha matrix run HMM-Flagger on all validation coverage files with more EM iterations and then make summary tables with the resolution of bases not windows (HMM-Flagger outputs summary tables with the resolution of windows) (Default = true)"
        flaggerPostOptConvergenceTol: "Convergence tolerance for running HMM-Flagger post optimization (Default = 0.001)"
        flaggerPostOptIterations : "Maxmimum number of EM iterations for running HMM-Flagger post optimization (Default = 50)"
    }
    input{
        Array[File] bamArray
        Array[File] baiArray
        Array[File] fastaArray
        Array[File] misassembliesBedArray
        Array[String] validationContigPatternArray=["chr20","chr21","chr22"]
        Array[Array[File]]? annotationBedArray2D
        Array[Array[String]]? biasAnnotationNameArray2D
        Array[Float] downsampleRateArray=[1.0]
        Boolean runFlaggerPostOptForTrain = false
        Boolean runFlaggerPostOptForValidation = true
        String flaggerDockerImage="mobinasri/flagger:v1.0.0_alpha--73711ad369ab956a33be9d14edf285a2191c2a94"
        String otherFlaggerParamsForTuning="--convergenceTol 1e-2  --iterations 10"
        Float flaggerPostOptConvergenceTol = 0.001
        Int flaggerPostOptIterations = 50
        String modelType="gaussian"
        File? binArrayTsv
        Int numberOfEgoStartPoints = 10
        Int numberOfEgoIterations = 50
        Int chunkLen = 20000000
        Int windowLen = 4000
        String labelNames = "Err,Dup,Hap,Col"
        Int maxFlaggerRunsPerIteration = 8
        Int threadsPerFlaggerRun = 8
    }
 
    # create two arrays of coverage files; one for training and one for validation
    scatter (indexAndRate in cross(range(length(bamArray)), downsampleRateArray)){
        Int fileIndex = indexAndRate.left
        Float rate = indexAndRate.right

        # Index diploid assembly
        call fai_t.produceFai {
            input:
                fasta = fastaArray[fileIndex]
        }

        ## create two lists of contigs one for training and one for validation
        call getTrainAndValidationContigText as getContigList{
            input:
                fai = produceFai.fai,
                validationContigPatternArray = validationContigPatternArray,
                dockerImage = flaggerDockerImage
        }

        if (defined(annotationBedArray2D)){
            Array[Array[File]] annotationBedArray2DValid = select_first([annotationBedArray2D])
            Array[File] annotationBedArray = annotationBedArray2DValid[fileIndex]
        }
        if (defined(biasAnnotationNameArray2D)){
            Array[Array[String]] biasAnnotationNameArray2DValid = select_first([biasAnnotationNameArray2D])
            Array[String] biasAnnotationNameArray = biasAnnotationNameArray2DValid[fileIndex]
        }

        ## TRAINING ##
        # convert bam to cov and add annotation indices to the cov file
        call cov_t.bam2cov as bam2covForTrain{
            input:
                bam = bamArray[fileIndex],
                bai = baiArray[fileIndex],
                fasta = fastaArray[fileIndex],
                suffix = "training",
                mapqThreshold = 20,
                clipRatioThreshold = 0.1,
                downsampleRate = rate,
                annotationBedArray = select_first([annotationBedArray, []]),
                biasAnnotationNameArray = select_first([biasAnnotationNameArray,[]]),
                baselineAnnotationName = "WHOLE_GENOME_DEFAULT",
                includeContigListText = getContigList.trainContigListText,
                runBiasDetection = (length(select_first([biasAnnotationNameArray,[]])) > 0),
                format = "all",
                memSize = 32,
                threadCount = 16,
                dockerImage = flaggerDockerImage
        }

        # Get coordinates of canonical bases only (no "N" which may come from scaffolding)
        call misc_t.getCanonicalBasesBed as dipCanonical{
            input:
                assemblyFasta = fastaArray[fileIndex],
                dockerImage = flaggerDockerImage
        }

        # make a truth bed file that contains
        # state indices instead of the names of the states
        call misc_t.getIndexLabeledBed as labelTruth{
                input:
                    bed = misassembliesBedArray[fileIndex],
                    canonicalBasesDiploidBed = dipCanonical.canonicalBasesBed,
                    addHapCoordinates = true
        }

        # Augment coverage file with prediction labels
        # it will add truth labels if available
        call augment_cov_t.augmentCoverageByLabels as augmentCovByTruthForTrain{
            input:
                coverage = bam2covForTrain.coverageGz,
                fai = produceFai.fai,
                numberOfLabels = 4,
                truthBed = labelTruth.labeledBed,
                includeContigListText = getContigList.trainContigListText,
                suffix="truth_added",
                dockerImage = flaggerDockerImage
        }

        ## VALIDATION ##
        # convert bam to cov and add annotation indices to the cov file
        call cov_t.bam2cov as bam2covForValidation{
            input:
                bam = bamArray[fileIndex],
                bai = baiArray[fileIndex],
                fasta = fastaArray[fileIndex],
                suffix = "validation",
                mapqThreshold = 20,
                clipRatioThreshold = 0.1,
                downsampleRate = rate,
                annotationBedArray = select_first([annotationBedArray, []]),
                biasAnnotationNameArray = select_first([biasAnnotationNameArray,[]]),
                baselineAnnotationName = "WHOLE_GENOME_DEFAULT",
                includeContigListText = getContigList.validationContigListText,
                runBiasDetection = (length(select_first([biasAnnotationNameArray,[]])) > 0),
                format = "all",
                memSize = 32,
                threadCount = 16,
                dockerImage = flaggerDockerImage
        }
        
        # Augment coverage file with truth labels for validation
        call augment_cov_t.augmentCoverageByLabels as augmentCovByTruthForValidation{
            input:
                coverage = bam2covForValidation.coverageGz,
                fai = produceFai.fai,
                numberOfLabels = 4,
                truthBed = labelTruth.labeledBed,
                includeContigListText = getContigList.validationContigListText,
                suffix="truth_added",
                dockerImage = flaggerDockerImage
        }

        call cov2bin_t.cov2bin as cov2binForTrain{
            input:
                coverage = augmentCovByTruthForTrain.augmentedCoverageGz,
                windowLen = windowLen,
                chunkLen = chunkLen,
                fai = produceFai.fai
        }
        
        call cov2bin_t.cov2bin as cov2binForValidation{
            input:
                coverage = augmentCovByTruthForValidation.augmentedCoverageGz,
                windowLen = windowLen,
                chunkLen = chunkLen,
                fai = produceFai.fai
        }
    }

    Int numberOfCoverageFiles = length(augmentCovByTruthForTrain.augmentedCoverageGz)

    Int memSizeForTuneAlpha = if (maxFlaggerRunsPerIteration < numberOfCoverageFiles) then (maxFlaggerRunsPerIteration * 8) else (numberOfCoverageFiles * 8)
    Int threadCountForTuneAlpha = if (maxFlaggerRunsPerIteration < numberOfCoverageFiles) then (maxFlaggerRunsPerIteration * threadsPerFlaggerRun) else (numberOfCoverageFiles * threadsPerFlaggerRun)  

    # find optimal alpha matrix
    call tune_alpha_t.tuneAlpha {
        input:
            trainBinArray = cov2binForTrain.bin, 
            validationBinArray = cov2binForValidation.bin,
            numberOfEgoStartPoints = numberOfEgoStartPoints,
            numberOfEgoIterations = numberOfEgoIterations,
            chunkLen = chunkLen,
            windowLen = windowLen,
            labelNames = labelNames,
            maxFlaggerRunsPerIteration = maxFlaggerRunsPerIteration,
            threadsPerFlaggerRun = threadsPerFlaggerRun,
            otherFlaggerParams = otherFlaggerParamsForTuning,
            modelType = modelType,
            memSize = memSizeForTuneAlpha,
            threadCount = threadCountForTuneAlpha,
            annotationLabel = "WHOLE_GENOME_DEFAULT",
            sizeLabel = "ALL_SIZES",
            dockerImage = flaggerDockerImage
    }

    #########################
    ## Post Optimization   ##
    ## Final Run for Train ##
    #########################
    if (runFlaggerPostOptForTrain){
        scatter(coverageIndex in range(numberOfCoverageFiles)){
            File trainCoverageGz = augmentCovByTruthForTrain.augmentedCoverageGz[coverageIndex]

            # We have optimal alpha matrix now
            # run HMM-Flagger on all coverage files
            ## TRAINING ##
            call hmm_flagger_t.hmmFlagger as hmmFlaggerPostOptTrain {
                input:
                    coverage = trainCoverageGz,
                    binArrayTsv = binArrayTsv,
                    chunkLen = chunkLen,
                    windowLen = windowLen,
                    labelNames = labelNames,
                    trackName = "hmm_flagger",
                    numberOfIterations = flaggerPostOptIterations,
                    convergenceTolerance = flaggerPostOptConvergenceTol,
                    alphaTsv = tuneAlpha.optimalAlphaTsv,
                    modelType = modelType,
                    memSize = 32,
                    threadCount = 16,
                    dockerImage = flaggerDockerImage
            }

            # make a prediction bed file that contains
            # state indices instead of the names of the states
            call misc_t.getIndexLabeledBed as labelPredictionTrain{
                input:
                    bed = hmmFlaggerPostOptTrain.predictionBed
            }

            # Augment coverage file with prediction labels
            # it will add truth labels if available
            call augment_cov_t.augmentCoverageByLabels as augmentCovByPredictionForTrain{
                input:
                    coverage = trainCoverageGz,
                    fai = produceFai.fai[coverageIndex],
                    numberOfLabels = 4,
                    predictionBed = labelPredictionTrain.labeledBed,
                    suffix = "pred_added",
                    dockerImage = flaggerDockerImage
            }

            call make_summary_table_t.makeSummaryTable as summaryForTrain{
                input:
                   coverage = augmentCovByPredictionForTrain.augmentedCoverageGz,
                   binArrayTsv = binArrayTsv,
                   dockerImage = flaggerDockerImage
            }
        }
    }


    ##############################
    ## Post Optimization        ##
    ## Final Run for Validation ##
    ##############################
    if (runFlaggerPostOptForValidation){
        scatter(coverageIndex in range(numberOfCoverageFiles)){
            File validationCoverageGz = augmentCovByTruthForValidation.augmentedCoverageGz[coverageIndex]

            ## VALIDATION ##
            call hmm_flagger_t.hmmFlagger as hmmFlaggerPostOptValidation {
                input:
                    coverage = validationCoverageGz,
                    binArrayTsv = binArrayTsv,
                    chunkLen = chunkLen,
                    windowLen = windowLen,
                    labelNames = labelNames,
                    trackName = "hmm_flagger",
                    numberOfIterations = flaggerPostOptIterations,
                    convergenceTolerance = flaggerPostOptConvergenceTol,
                    alphaTsv = tuneAlpha.optimalAlphaTsv,
                    modelType = modelType,
                    memSize = 32,
                    threadCount = 16,
                    dockerImage = flaggerDockerImage
            }
            # make a prediction bed file that contains
            # state indices instead of the names of the states
            call misc_t.getIndexLabeledBed as labelPredictionValidation{
                input:
                    bed = hmmFlaggerPostOptValidation.predictionBed
            }

            # Augment coverage file with prediction labels
            # it will add truth labels if available
            call augment_cov_t.augmentCoverageByLabels as augmentCovByPredictionForValidation{
                input:
                    coverage = validationCoverageGz,
                    fai = produceFai.fai[coverageIndex],
                    numberOfLabels = 4,
                    predictionBed = labelPredictionValidation.labeledBed,
                    suffix="pred_added",
                    dockerImage = flaggerDockerImage
            }

            call make_summary_table_t.makeSummaryTable as summaryForValidation{
                input:
                   coverage = augmentCovByPredictionForValidation.augmentedCoverageGz,
                   binArrayTsv = binArrayTsv,
                   dockerImage = flaggerDockerImage
            }
        }
    }


    output {

        # Output files realted to optimizing alpha matrix

        File optimalAlphaTsv = tuneAlpha.optimalAlphaTsv
        File optimizationOutputDataTarGz = tuneAlpha.optimizationOutputDataTarGz
        Float trainOptimalScore = tuneAlpha.trainOptimalScore
        Float validationOptimalScore = tuneAlpha.validationOptimalScore
        File trainScoresTsv = tuneAlpha.trainScoresTsv
        File validationScoresTsv = tuneAlpha.validationScoresTsv
        File trainFilesListText = tuneAlpha.trainFilesListText
        File validationFilesListText = tuneAlpha.validationFilesListText
        
        # HMM-Flagger results for training data after optimizing alpha
        # They will be present only if runFlaggerPostOptTrain is true

        Array[File]? benchmarkingSummaryTsvArrayForTrain = summaryForTrain.benchmarkingSummaryTsv
        Array[File]? contiguitySummaryTsvArrayForTrain = summaryForTrain.contiguitySummaryTsv
        Array[File]? fullStatsTsvArrayForTrain = summaryForTrain.fullStatsTsv
        Array[File]? predictionBedArrayForTrain = hmmFlaggerPostOptTrain.predictionBed

        # HMM-Flagger results for validation data after optimizing alpha
        # They will be present only if runFlaggerPostOptValidation is true

        Array[File]? benchmarkingSummaryTsvArrayForValidation = summaryForValidation.benchmarkingSummaryTsv
        Array[File]? contiguitySummaryTsvArrayForValidation = summaryForValidation.contiguitySummaryTsv
        Array[File]? fullStatsTsvArrayForValidation = summaryForValidation.fullStatsTsv
        Array[File]? predictionBedArrayForValidation = hmmFlaggerPostOptValidation.predictionBed

    }
}


task getTrainAndValidationContigText{
    input{
        File fai
        Array[String] validationContigPatternArray
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=8
        String dockerImage="mobinasri/flagger:v1.0.0_alpha"
        Int preemptible=2
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace
        
        for pattern in ~{sep=" " validationContigPatternArray}
        do
            echo ${pattern} >> validation_patterns.txt
        done

        PREFIX=$(echo $(basename ~{fai}) | sed -e 's/\.fa.fai$//' -e 's/\.fasta.fai$//')

        mkdir output

        cat ~{fai} | \
            cut -f1 | \
            grep -v -f validation_patterns.txt > output/${PREFIX}.train_contigs.txt

        cat ~{fai} | \
            cut -f1 | \
            grep -f validation_patterns.txt > output/${PREFIX}.validation_contigs.txt

    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output{
        File trainContigListText = glob("output/*.train_contigs.txt")[0]
        File validationContigListText = glob("output/*.validation_contigs.txt")[0]
    }
}



