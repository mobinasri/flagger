version 1.0 

workflow runCov2CountsByWindow{
    call cov2countsByWindow
    output {
        File windowCountsTarGz = cov2countsByWindow.windowCountsTarGz
        File windowCovsTarGz = cov2countsByWindow.windowCovsTarGz
        File windowsText = cov2countsByWindow.windowsText
    }
}

task cov2countsByWindow {
    input {
        File coverageGz
        Int windowSize
        File fai
        Array[File] excludeBedArray
        # runtime configurations
        Int memSize=8
        Int threadCount=8
        Int diskSize=64
        String dockerImage="mobinasri/flagger:dev-v0.1"
        Int preemptible=2
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

        FILENAME=$(basename ~{fai})
        PREFIX_FAI=${FILENAME%.fa.fai}

        FILENAME=$(basename ~{coverageGz})
        PREFIX_COV=${FILENAME%.cov.gz}
        
        gunzip -c ~{coverageGz} > ${PREFIX_COV}.cov

        # Make a bed file of the included regions
        cat ~{fai} | awk '{print $1"\t0\t"$2}' | sort -k1,1V -k2,2n > asm.bed
        cat ~{sep=" " excludeBedArray} | sort -k1,1V -k2,2n | bedtools merge -i - > exclude.bed || true
        bedtools subtract -a asm.bed -b exclude.bed > asm.excluded.bed

        # Remove excluded regions from cov file
        cat ${PREFIX_COV}.cov | \
            awk '{if(substr($1,1,1) == ">") {contig=substr($1,2); len_contig=$2} else {print contig"\t"$1-1"\t"$2"\t"$3"\t"len_contig}}' | \
            bedtools intersect -a - -b asm.excluded.bed | \
            awk '{if(contig != $1){contig=$1; print ">"contig" "$5}; print $2+1"\t"$3"\t"$4}' > ${PREFIX_COV}.excluded.cov
        
        # make a tab-delimited file including the contig names and their effective length (not excluded)
        cat asm.excluded.bed | awk '{ctg_len[$1] += $3-$2}END{for (c in ctg_len){print c"\t"ctg_len[c]}}' > ctg_lens.txt
        mkdir covs counts
        # Make a separate cov file for each window
        split_cov_by_window_test -c ${PREFIX_COV}.excluded.cov -f ctg_lens.txt -p covs/${PREFIX_COV} -s ~{windowSize} > ${PREFIX_FAI}.windows.txt
        # Count each window-specific cov file
        for c in $(ls covs);do cov2counts -i covs/$c -o counts/${c/.cov/.counts}; echo $c" finished";done

        # Compress Counts files
        tar -cf ${PREFIX_COV}.counts.tar counts
        gzip ${PREFIX_COV}.counts.tar

        # Compress Cov files
        tar -cf ${PREFIX_COV}.covs.tar covs
        gzip ${PREFIX_COV}.covs.tar
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File windowCountsTarGz = glob("*.counts.tar.gz")[0]
        File windowCovsTarGz = glob("*.covs.tar.gz")[0]
        File windowsText = glob("*.windows.txt")[0]
    }
}

