version 1.0

workflow runBusco {
	call busco
}

task busco {
    input {
        File assemblyFasta
        String? customDatasetURL
        String extraArguments="-m geno -sp human"
        Int memSizeGB = 32
        Int threadCount = 16
        Int diskSizeGB = 64
        String dockerImage = "tpesout/hpp_busco:latest"
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

        # file initialization
        FILENAME=$(basename -- "~{assemblyFasta}")
        if [[ $FILENAME =~ \.gz$ ]]; then
            cp ~{assemblyFasta} .
            gunzip $FILENAME
            FILENAME="${FILENAME%\.gz}"
        fi
        PREFIX="${FILENAME%.*}"
        SUFFIX="${FILENAME##*.}"

        # environment variables for config files needed by BUSCO
        export BUSCO_CONFIG_FILE="/root/tools/BUSCO/busco/config/config.ini"
        export AUGUSTUS_CONFIG_PATH="/root/tools/Augustus/Augustus-3.3.1/config/"

        # initialize script
        cmd=(  python3 /root/tools/BUSCO/busco/scripts/run_BUSCO.py )
        cmd+=( -i $FILENAME )
        cmd+=( -o $PREFIX.busco )
        cmd+=( -c ~{threadCount} )

        # include reference fasta if supplied
        if [[ ! -z "~{customDatasetURL}" ]]; then
            DATASET_FILE=$(basename ~{customDatasetURL})
            DATASET="${DATASET_FILE%.tar.gz}"
            mkdir dataset
            cd dataset
            wget ~{customDatasetURL}
            tar xvf $DATASET_FILE
            rm $DATASET_FILE
            cd ..
            cmd+=( -l dataset/$DATASET )
        else
            cmd+=( -l /root/tools/BUSCO/dataset/vertebrata_odb9 )
        fi

        # include extra arguments
        if [[ ! -z "~{extraArguments}" ]] ; then
            cmd+=( ~{extraArguments} )
        fi

        # run command
        "${cmd[@]}"

        tar czvf $PREFIX.busco.tar.gz run_$PREFIX.busco

	>>>
	output {
		File outputTarball = glob("*.busco.tar.gz")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
