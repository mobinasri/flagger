version 1.0

workflow runTuneAlpha{
    call tuneAlpha
    output{
        File optimalAlphaTsv = tuneAlpha.optimalAlphaTsv
        File optimizationOutputDataTarGz = tuneAlpha.optimizationOutputDataTarGz
        Float trainOptimalScore = tuneAlpha.trainOptimalScore
        Float validationOptimalScore = tuneAlpha.validationOptimalScore
        File trainScoresTsv = tuneAlpha.trainScoresTsv
        File validationScoresTsv = tuneAlpha.validationScoresTsv
    }
}


task tuneAlpha{
    input{
        Array[File] trainBinArray
        Array[File] validationBinArray
        File? binArrayTsv
        Int numberOfEgoStartPoints = 10
        Int numberOfEgoIterations = 50
        Int chunkLen = 20000000
        Int windowLen = 4000
        String labelNames = "Err,Dup,Hap,Col"
        # maximum number of hmm-flagger jobs that can be executed per EGO iteration
        Int maxFlaggerRunsPerIteration = 8
        # number of threads for each hmm-flagger job
        Int threadsPerFlaggerRun = 8
        String otherFlaggerParams = "--convergenceTol 1e-2  --iterations 10"
        File? candidateAlphaTsv
        String modelType = "gaussian"
        String? annotationLabel = "WHOLE_GENOME_DEFAULT"
        String? sizeLabel = "ALL_SIZES"
        # runtime configurations
        Int memSize=64
        Int threadCount=64
        Int diskSize=ceil(size(trainBinArray, "GB"))  + ceil(size(validationBinArray, "GB")) + 64
        String dockerImage="mobinasri/flagger:v1.0_alpha"
        Int preemptible=2
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace


        mkdir -p train_bin_files
        mkdir -p validation_bin_files

        for TRAIN_BIN_FILE in ~{sep=" " trainBinArray}
        do
            FILENAME=$(basename ${TRAIN_BIN_FILE})
            ln -s ${TRAIN_BIN_FILE} train_bin_files/${FILENAME}
        done

        for VALIDATION_BIN_FILE in ~{sep=" " validationBinArray}
        do
            FILENAME=$(basename ${VALIDATION_BIN_FILE})
            ln -s ${VALIDATION_BIN_FILE} validation_bin_files/${FILENAME}
        done

        function join_by { local IFS="$1"; shift; echo "$*"; }

        TRAIN_BIN_FILES=$(find train_bin_files/ | grep ".bin$" | sort)
        VALIDATION_BIN_FILES=$(find  validation_bin_files/ | grep ".bin$" | sort)

        TRAIN_BIN_FILES_STR=$(join_by , ${TRAIN_BIN_FILES})
        VALIDATION_BIN_FILES_STR=$(join_by , ${VALIDATION_BIN_FILES})

        touch train_files.txt
        touch validation_files.txt

        for TRAIN_BIN_FILE in ${TRAIN_BIN_FILES[@]}
        do
            echo $(basename ${TRAIN_BIN_FILE%%.bin}) >> train_files.txt
        done

        for VALIDATION_BIN_FILE in ${VALIDATION_BIN_FILES[@]}
        do
            echo $(basename ${VALIDATION_BIN_FILE%%.bin}) >> validation_files.txt
        done

        OUTPUT_DIR="tune_alpha_~{modelType}_w_~{windowLen}_n_~{numberOfEgoIterations}"
        mkdir -p ${OUTPUT_DIR}

        
        echo -n "--chunkLen ~{chunkLen}  --windowLen ~{windowLen}  --labelNames ~{labelNames} " > ${OUTPUT_DIR}/other_params.txt
        echo -n " --threads ~{threadsPerFlaggerRun}  ~{otherFlaggerParams} " >> ${OUTPUT_DIR}/other_params.txt
        if [ -n "~{binArrayTsv}" ]
        then 
            echo -n " --binArrayFile ~{binArrayTsv}" > ${OUTPUT_DIR}/other_params.txt >> ${OUTPUT_DIR}/other_params.txt
        fi
        # just for a new line
        echo "" >> ${OUTPUT_DIR}/other_params.txt

        ADDITIONAL_ARGS=""
        if [ -n "~{candidateAlphaTsv}" ]
        then
            ADDITIONAL_ARGS="--candidateAlphaTsv ~{candidateAlphaTsv}"
        fi

        if [ -n "~{annotationLabel}" ]
        then
            ADDITIONAL_ARGS="${ADDITIONAL_ARGS} --annotationLabel ~{annotationLabel}"
        fi

        if [ -n "~{sizeLabel}" ]
        then
            ADDITIONAL_ARGS="${ADDITIONAL_ARGS} --sizeLabel ~{sizeLabel}"
        fi

        python3 /home/programs/src/tune_alpha_hmm_flagger.py \
            --inputFilesTrain ${TRAIN_BIN_FILES_STR} \
            --inputFilesValidation ${VALIDATION_BIN_FILES_STR} \
            --outputDir ${OUTPUT_DIR} \
            --otherParamsText ${OUTPUT_DIR}/other_params.txt \
            --numberOfStartPoints ~{numberOfEgoStartPoints} \
            --modelType ~{modelType} \
            --iterations ~{numberOfEgoIterations} \
            --maxJobs ~{maxFlaggerRunsPerIteration} ${ADDITIONAL_ARGS}
        
        mkdir -p output
        cp ${OUTPUT_DIR}/alpha_optimum.tsv output/alpha_optimum_~{modelType}_w_~{windowLen}_n_~{numberOfEgoIterations}.tsv
        cp ${OUTPUT_DIR}/scores_validation.tsv output/scores_validation_~{modelType}_w_~{windowLen}_n_~{numberOfEgoIterations}.tsv
        cp ${OUTPUT_DIR}/scores_train.tsv output/scores_train_~{modelType}_w_~{windowLen}_n_~{numberOfEgoIterations}.tsv
       
        cat ${OUTPUT_DIR}/scores_validation.tsv | awk '$4=="YES"{print $3}' > output/validation_optimum_score.txt
        cat ${OUTPUT_DIR}/scores_train.tsv | awk '$4=="YES"{print $3}' > output/train_optimum_score.txt

        tar -cf  ${OUTPUT_DIR}.tar ${OUTPUT_DIR}
        pigz -p8 ${OUTPUT_DIR}.tar 
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output{
        File optimalAlphaTsv = glob("output/alpha_optimum_*.tsv")[0]
        File optimizationOutputDataTarGz = glob("tune_alpha_*.tar.gz")[0]
        Float trainOptimalScore = read_float("output/train_optimum_score.txt")
        Float validationOptimalScore = read_float("output/validation_optimum_score.txt")
        File trainScoresTsv = glob("output/scores_train_*.tsv")[0]
        File validationScoresTsv = glob("output/scores_validation_*.tsv")[0]
        File trainFilesListText = glob("train_files.txt")[0]
        File validationFilesListText = glob("validation_files.txt")[0]
    }
}


