version 1.0

workflow runAugmentCoverageByLabels{
    call augmentCoverageByLabels
    output{
        File augmentedCoverageGz = augmentCoverageByLabels.augmentedCoverageGz
    }
}


task augmentCoverageByLabels{
    input{
        File coverage
        File fai
        Int numberOfLabels
        File? truthBed
        File? predictionBed
        String suffix="augmented"
        File? includeContigListText
        # runtime configurations
        Int memSize=32
        Int threadCount=8
        Int diskSize=ceil(size(coverage, "GB"))  + 64
        String dockerImage="mobinasri/flagger:v1.0.0"
        Int preemptible=2
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        COV_FILE_PATH="~{coverage}"
        PREFIX=$(echo $(basename ${COV_FILE_PATH%.gz}) | sed -e 's/\.cov$//' -e 's/\.bed$//')
        EXTENSION=${COV_FILE_PATH##*.}

        INPUT_BED_ARGS=""
        if [ -n "~{truthBed}" ]; then
            # filter contigs
            if [ -n "~{includeContigListText}" ]; then
                cat ~{truthBed} | grep -F -f ~{includeContigListText} > truth.bed
            else
                cp ~{truthBed} truth.bed 
            fi
            INPUT_BED_ARGS="${INPUT_BED_ARGS} --truthBed truth.bed"     
        fi
        if [ -n "~{predictionBed}" ]; then
            # filter contigs
            if [ -n "~{includeContigListText}" ]; then
                cat ~{predictionBed} | grep -F -f ~{includeContigListText} > prediction.bed
            else
                cp ~{predictionBed} prediction.bed
            fi
            INPUT_BED_ARGS="${INPUT_BED_ARGS} --predictionBed prediction.bed"
        fi

        mkdir output
        OUTPUT="output/${PREFIX}.~{suffix}.cov.gz"

        if [ -n "${INPUT_BED_ARGS}" ]
        then
            augment_coverage_by_labels \
                -i ${COV_FILE_PATH} \
                -o ${OUTPUT} \
                ${INPUT_BED_ARGS} \
                --numberOfLabels ~{numberOfLabels} \
                --fai ~{fai} \
                --threads ~{threadCount}
        else
           if [ "${EXTENSION}" == "gz" ]
           then
               cp ~{coverage} output/${PREFIX}.cov.gz
           else
               pigz -c -p ~{threadCount} ~{coverage} > output/${PREFIX}.cov.gz
           fi
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
        File augmentedCoverageGz = glob("output/*.cov.gz")[0]
    }
}


