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
        String preset
        String suffix
        File? binArrayTsv
        String labelNames = "Err,Dup,Hap,Col"
        String trackName = "hmm_flagger"
        Int numberOfIterations = 100
        Float convergenceTolerance = 0.001
        Float maxHighMapqRatio=0.25
        Float minHighMapqRatio=0.75
        String? moreOptions
        File? alphaTsv
        String? modelType
        Int? chunkLen
        Int? windowLen
        Array[Int] minimumBlockLenArray = []
        # runtime configurations
        Int memSize=32
        Int threadCount=8
        Int diskSize=ceil(size(coverage, "GB")) + 64
        String dockerImage="mobinasri/flagger:v1.2.0"
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

        if [ -n "~{chunkLen}" ]
        then
            ADDITIONAL_ARGS="${ADDITIONAL_ARGS} --chunkLen ~{chunkLen}"
        fi

        if [ -n "~{windowLen}" ]
        then
            ADDITIONAL_ARGS="${ADDITIONAL_ARGS} --windowLen ~{windowLen}"
        fi

        if [ -n "~{modelType}" ]
        then
            ADDITIONAL_ARGS="${ADDITIONAL_ARGS} --modelType ~{modelType}"
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
            --preset ~{preset} \
            --outputDir ${OUTPUT_DIR}  \
            --maxHighMapqRatio ~{maxHighMapqRatio} \
            --minHighMapqRatio ~{minHighMapqRatio} \
            --iterations ~{numberOfIterations} \
            --trackName ~{trackName} \
            --convergenceTol ~{convergenceTolerance} \
            --labelNames ~{labelNames} \
            --threads ~{threadCount} ${ADDITIONAL_ARGS}
        
        mkdir -p output
        cp ${OUTPUT_DIR}/*.bed output/${PREFIX}.~{suffix}_prediction.bed
        cp ${OUTPUT_DIR}/prediction_summary_final.tsv output/${PREFIX}.~{suffix}.prediction_summary_final.tsv
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



task filterHmmFlaggerCalls{
    input{
        File selfAsmMapBam
        File flaggerBed
        Int minAlignmentLen=10000
        Float maxDivergence=0.005
        String? moreOptions
        # runtime configurations
        Int memSize=32
        Int threadCount=8
        Int diskSize=32
        String dockerImage="mobinasri/flagger:v1.2.0"
        Int preemptible=2
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        IN_BED_PATH="~{flaggerBed}"
        PREFIX=$(basename ${IN_BED_PATH%.bed})

        mkdir -p output

        python3 /home/programs/src/filter_hmm_flagger_calls.py \
            --inputBam ~{selfAsmMapBam} \
            --inputBed ~{flaggerBed} \
            --outputBed output/${PREFIX}.conservative.bed \
            --minAlignmentLen ~{minAlignmentLen} \
            --maxDivergence ~{maxDivergence} \
            --threads ~{threadCount} ~{moreOptions}
        
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output{
        File conservativeBed = glob("output/*.conservative.bed")[0]
    }
}


