version 1.0

workflow runMerqury {
    call merqury
}

task merqury {
    input {
        File assemblyFasta
        File? altHapFasta
        File kmerTarball
        File? matKmerTarball
        File? patKmerTarball
        Int memSize = 12
        Int threadCount = 16
        Int diskSizeGB = 256
        String dockerImage = "juklucas/hpp_merqury:latest"
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
        OMP_NUM_THREADS=~{threadCount}

        # get filename
        ASM_ID=$(basename ~{assemblyFasta} | sed 's/.gz$//' | sed 's/.fa\(sta\)*$//' | sed 's/[._][pm]at\(ernal\)*//')

        # extract kmers
        tar xvf ~{kmerTarball} &
        if [[ -f "~{matKmerTarball}" && -f "~{patKmerTarball}" ]]; then
            tar xvf ~{matKmerTarball} &
            tar xvf ~{patKmerTarball} &
        fi
        wait

        # initilize command
        cmd=( merqury.sh )
        cmd+=( $(basename ~{kmerTarball} | sed 's/.gz$//' | sed 's/.tar$//') )
        if [[ -f "~{matKmerTarball}" && -f "~{patKmerTarball}" ]]; then
            cmd+=( $(basename ~{matKmerTarball} | sed 's/.gz$//' | sed 's/.tar$//') )
            cmd+=( $(basename ~{patKmerTarball} | sed 's/.gz$//' | sed 's/.tar$//') )
        fi

        # link primary asm file
        FILENAME=$(basename -- "~{assemblyFasta}")
        if [[ $FILENAME =~ \.gz$ ]]; then
            cp ~{assemblyFasta} .
            gunzip $FILENAME
            mv ${FILENAME%\.gz} asm.fasta
        else
            ln -s ~{assemblyFasta} asm.fasta
        fi
        cmd+=( asm.fasta )

        # link secondary asm file
        if [[ -f "~{altHapFasta}" ]]; then
            FILENAME=$(basename -- "~{altHapFasta}")
            if [[ $FILENAME =~ \.gz$ ]]; then
                cp ~{altHapFasta} .
                gunzip $FILENAME
                mv ${FILENAME%\.gz} altHap.fasta
            else
                ln -s ~{altHapFasta} altHap.fasta
            fi
            cmd+=( altHap.fasta )
        fi

        # prep output
        cmd+=( $ASM_ID.merqury )

        # run command
        ${cmd[@]}

        # get output
        tar czvf $ASM_ID.merqury.tar.gz $ASM_ID.merqury*

	>>>
	output {
		File QV = glob("*.merqury.qv")[0]
		File outputTarball = glob("*.merqury.tar.gz")[0]
	}
    runtime {
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}

