version 1.0

workflow runShardReads {
    call shardReads
}

task shardReads {
    input {
        File readFile
        Int linesPerFile = 256000000
        Int memSizeGB = 4
        Int threadCount = 1
        Int diskSizeGB = 128
        String dockerImage = "tpesout/hpp_base:latest"
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

        ORIG_FILENAME=$(basename -- "~{readFile}")
        PREFIX="${ORIG_FILENAME%.*}"
        SUFFIX="${ORIG_FILENAME##*.}"

        # remove gz
        if [[ $SUFFIX == "gz" ]] ; then
            SUFFIX="${PREFIX##*.}"
            PREFIX="${PREFIX%.*}"
        fi

        # prep output
        mkdir -p split
        OUTPUT_PREFIX="split/$PREFIX."

        # extract with split
        if [[ ${ORIG_FILENAME##*.} == "gz" ]]; then
            zcat ~{readFile} | split -a 4 -d -l ~{linesPerFile} --additional-suffix=".$SUFFIX" - $OUTPUT_PREFIX
        else
            split -a 4 -d -l ~{linesPerFile} --additional-suffix=".$SUFFIX" ~{readFile} $OUTPUT_PREFIX
        fi

        # get size of first file
        du -s -BG split/$(ls split/ | head -n1) | sed 's/G.*//' >outputsize
    >>>

    output {
        Array[File] shardedReads = glob("split/*")
        Int fileSizeGB = read_int("outputsize")
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
