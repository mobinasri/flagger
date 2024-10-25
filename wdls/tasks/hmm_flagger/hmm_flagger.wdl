version 1.0

workflow runHmmFlagger{
    call hmmFlagger
    output{
        File predictionBed = hmmFlagger.predictionBed
        File predictionSummaryTsv = hmmFlagger.predictionSummaryTsv
        File loglikelihoodTsv = hmmFlagger.loglikelihoodTsv
        File outputTarGz = hmmFlagger.outputTarGz
    }
}


task hmmFlagger{
    input{
        File coverage
        File? binArrayTsv
        Int chunkLen = 20000000
        Int windowLen = 4000
        String labelNames = "Err,Dup,Hap,Col"
        String trackName = "hmm_flagger_v1.0"
        Int numberOfIterations = 100
        Float convergenceTolerance = 0.001
        Float maxHighMapqRatio=0.25
        Float minHighMapqRatio=0.5
        String? moreOptions
        File? alphaTsv
        String modelType = "gaussian"
        Array[Int] minimumBlockLenArray = []
        # runtime configurations
        Int memSize=32
        Int threadCount=8
        Int diskSize=ceil(size(coverage, "GB")) + 64
        String dockerImage="mobinasri/flagger:v1.1.0"
        Int preemptible=2
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace


        COV_FILE_PATH="~{coverage}"
        PREFIX=$(echo $(basename ${COV_FILE_PATH%.gz}) | sed -e 's/\.cov$//' -e 's/\.bed$//')

        ADDITIONAL_ARGS=""
        if [ -n "~{alphaTsv}" ]
        then
            ADDITIONAL_ARGS="--alpha ~{alphaTsv}"
        fi

        if [ -n "~{binArrayTsv}" ]
        then
            ADDITIONAL_ARGS="${ADDITIONAL_ARGS} --binArrayFile ~{binArrayTsv}"
        fi 
        
        if [ -n "~{moreOptions}" ]
        then
            ADDITIONAL_ARGS="${ADDITIONAL_ARGS} ~{moreOptions}"
        fi
        
        if [ ~{length(minimumBlockLenArray)} -gt "0" ]
        then
            ADDITIONAL_ARGS="${ADDITIONAL_ARGS} --minimumLengths ~{sep=',' minimumBlockLenArray}"
        fi
 
        OUTPUT_DIR="${PREFIX}"
        mkdir -p ${OUTPUT_DIR}

        hmm_flagger \
            --input ~{coverage} \
            --outputDir ${OUTPUT_DIR}  \
            --chunkLen ~{chunkLen} \
            --windowLen ~{windowLen} \
            --maxHighMapqRatio ~{maxHighMapqRatio} \
            --minHighMapqRatio ~{minHighMapqRatio} \
            --modelType ~{modelType} \
            --iterations ~{numberOfIterations} \
            --trackName ~{trackName} \
            --convergenceTol ~{convergenceTolerance} \
            --labelNames ~{labelNames} \
            --threads ~{threadCount} ${ADDITIONAL_ARGS}
        
        mkdir -p output
        cp ${OUTPUT_DIR}/*.bed output/${PREFIX}.hmm_flagger_prediction.bed
        cp ${OUTPUT_DIR}/prediction_summary_final.tsv output/${PREFIX}.prediction_summary_final.tsv
        cp ${OUTPUT_DIR}/loglikelihood.tsv output/

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
        File predictionBed = glob("output/*.bed")[0]
        File predictionSummaryTsv = glob("output/*.prediction_summary_final.tsv")[0]
        File loglikelihoodTsv = glob("output/loglikelihood.tsv")[0]
        File outputTarGz = glob("*.tar.gz")[0]
    }
}


