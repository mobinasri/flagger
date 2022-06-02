version 1.0

workflow runAsmgene {
	call asmgene
}

task asmgene {
    input{
        File assemblyFasta
        File genesFasta
        File? referenceFasta
        File? genesToReferencePaf
        # runtime configurations
        Int threadCount = 32
        Int memSizeGB = 64
        Int diskSizeGB = 32
        String dockerImage = "tpesout/hpp_minimap2:latest"
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

        # name prep
        FILENAME=$(basename -- "~{assemblyFasta}" | sed 's/.gz$//' )
        PREFIX="${FILENAME%.*}"

        # determine if we need to align genes to reference
        if [[ -f "~{genesToReferencePaf}" ]] ; then
            ln -s ~{genesToReferencePaf} genesToRef.paf
        elif [[ -f "~{referenceFasta}" ]] ; then
            minimap2 -cx splice:hq -t ~{threadCount} ~{referenceFasta} ~{genesFasta} > genesToRef.paf
        else
            echo "Either referenceFasta or genesToReferencePaf needs to be supplied"
            exit 1
        fi

        # Aligning genes to assembly
        minimap2 -cx splice:hq -t ~{threadCount} ~{assemblyFasta} ~{genesFasta} > $PREFIX.paf

        # Computing statistics for gene completeness
        paftools.js asmgene -a genesToRef.paf $PREFIX.paf > $PREFIX.gene_stats.txt
    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible: 1
    }

    output {
        File geneStats = glob("*.gene_stats.txt")[0]
    }
}

