version 1.0

workflow runMakeSummaryTable{
    call makeSummaryTable
    output{
        File benchmarkingSummaryTsv = makeSummaryTable.benchmarkingSummaryTsv
        File contiguitySummaryTsv = makeSummaryTable.contiguitySummaryTsv
        File fullStatsTsv = makeSummaryTable.fullStatsTsv
    }
}


task makeSummaryTable{
    input{
        File coverage
        String labelNames="Err,Dup,Hap,Col"
        File? binArrayTsv
        # runtime configurations
        Int memSize=32
        Int threadCount=8
        Int diskSize=ceil(size(coverage, "GB"))  + 64
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
        if [ -n "~{binArrayTsv}" ]; then
            ADDITIONAL_ARGS="--binArrayFile ~{binArrayTsv}"
        fi


        mkdir output
        OUTPUT="output/${PREFIX}.stats.tsv"

        # make empty tsv files just to avoid not-existing-file error
        # stats.benchmarking.tsv and stats.benchmarking.auN_ratio.tsv are
        # created only when truth and prediction labels are both available
        # in the given coverage file 
        touch output/${PREFIX}.stats.tsv
        touch output/${PREFIX}.stats.benchmarking.tsv
        touch output/${PREFIX}.stats.benchmarking.auN_ratio.tsv

        make_summary_table \
            --input ${COV_FILE_PATH} \
            --output ${OUTPUT} \
            --labelNames ~{labelNames} \
            -@ ~{threadCount} ${ADDITIONAL_ARGS}

    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output{
        File benchmarkingSummaryTsv = glob("output/*.benchmarking.tsv")[0]
        File contiguitySummaryTsv = glob("output/*.benchmarking.auN_ratio.tsv")[0]
        File fullStatsTsv = glob("output/*.stats.tsv")[0]
    }
}


