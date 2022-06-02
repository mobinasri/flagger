version 1.0

workflow runDipcall {
	call dipcall
}

task dipcall {
    input {
        File assemblyFastaPat
        File assemblyFastaMat
        File referenceFasta
        Boolean isMaleSample
        Boolean referenceIsHS38 = true
        Int memSizeGB = 64
        Int threadCount = 16
        Int diskSizeGB = 64
        String dockerImage = "humanpangenomics/hpp_dipcall_v0.3:latest"
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
        PATH="/root/bin/samtools_1.9:$PATH"

        # get output base
        PREFIX=$(basename ~{assemblyFastaPat} | sed 's/.gz$//' | sed 's/.fa\(sta\)*$//' | sed 's/[._][pm]at\(ernal\)*//')
        mkdir $PREFIX.dipcall

        # prep paternal
        PAT_FILENAME=$(basename -- "~{assemblyFastaPat}")
        if [[ $PAT_FILENAME =~ \.gz$ ]]; then
            cp ~{assemblyFastaPat} .
            gunzip $PAT_FILENAME
            PAT_FILENAME="${PAT_FILENAME%.gz}"
        else
            ln -s ~{assemblyFastaPat}
        fi

        # prep maternal
        MAT_FILENAME=$(basename -- "~{assemblyFastaMat}")
        if [[ $MAT_FILENAME =~ \.gz$ ]]; then
            cp ~{assemblyFastaMat} .
            gunzip $MAT_FILENAME
            MAT_FILENAME="${MAT_FILENAME%.gz}"
        else
            ln -s ~{assemblyFastaMat}
        fi

        # prep reference
        REF_FILENAME=$(basename -- "~{referenceFasta}")
        if [[ $REF_FILENAME =~ \.gz$ ]]; then
            cp ~{referenceFasta} .
            gunzip $REF_FILENAME
            REF_FILENAME="${REF_FILENAME%.gz}"
        else
            ln -s ~{referenceFasta}
        fi
        samtools faidx $REF_FILENAME

        # initialize script
        cmd=( /opt/dipcall/dipcall.kit/run-dipcall )

        # male samples need PAR region excluded
        if [[ ~{isMaleSample} == true ]]; then
            if [[ ~{referenceIsHS38} ]]; then
                cmd+=( -x /opt/dipcall/dipcall.kit/hs38.PAR.bed )
            else
                cmd+=( -x /opt/dipcall/dipcall.kit/hs37d5.PAR.bed )
            fi

        fi

        # finalize script
        cmd+=( $PREFIX.dipcall/$PREFIX )
        cmd+=( $REF_FILENAME )
        cmd+=( $PAT_FILENAME )
        cmd+=( $MAT_FILENAME )

        # generate makefile
        "${cmd[@]}" >$PREFIX.mak

        # run dipcall
        make -j 2 -f $PREFIX.mak

        # finalize
        rm $PREFIX.dipcall/*sam.gz
        tar czvf $PREFIX.dipcall.tar.gz $PREFIX.dipcall/
        cp $PREFIX.dipcall/$PREFIX.dip.bed $PREFIX.dipcall.bed
        cp $PREFIX.dipcall/$PREFIX.dip.vcf.gz $PREFIX.dipcall.vcf.gz

        # cleanup
        rm $REF_FILENAME
        rm $MAT_FILENAME
        rm $PAT_FILENAME

	>>>
	output {
		File outputTarball = glob("*.dipcall.tar.gz")[0]
		File outputVCF = glob("*.dipcall.vcf.gz")[0]
		File outputBED = glob("*.dipcall.bed")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
