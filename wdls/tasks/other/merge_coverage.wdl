version 1.0 

workflow runMergeCoverage{
    call mergeCoverage
    output{
        File mergedCovGz = mergeCoverage.mergedCovGz
    }
}
task mergeCoverage {
    input {
        Array[File] beds
        File fai
        File covGz
        File highMapqCovGz
        # runtime configurations
        Int memSize=8
        Int threadCount=4
        Int diskSize=512
        String dockerImage="quay.io/masri2019/hpp_coverage:latest"
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
        

        touch union.bed
        counter=1
        for bed in ~{sep=" " beds}
        do
            bedtools subtract -a ${bed} -b union.bed | \
                bedtools sort -i - | \
                bedtools merge -i - | \
                awk -vc=${counter} '{print $0"\t"c}' > r_${counter}.bed
            cat ${bed} union.bed | \
                bedtools sort -i - | \
                bedtools merge -i - > tmp.bed && mv tmp.bed union.bed
            ((++counter))
        done
        cat ~{fai} | awk '{print $1"\t0\t"$2}' > asm.bed
        bedtools subtract -a asm.bed -b union.bed | \
            bedtools sort -i - | \
            bedtools merge -i - | \
            awk '{print $0"\t"0}'> r_0.bed
        cat r_*.bed | sort -k1,1V -k2,2n > regions.bed

        # put chrM last
        cat regions.bed | grep -v "chrM" > regions_no_M.bed
        cat regions.bed | grep "chrM" > regions_only_M.bed
        cat regions_no_M.bed regions_only_M.bed > regions.bed

        FILENAME=$(basename ~{covGz})
        PREFIX=${FILENAME%%.cov.gz}

        gunzip -c ~{covGz} > all.cov
        gunzip -c ~{highMapqCovGz} > high_mapq.cov
        merge_blocks -a all.cov -b high_mapq.cov -o merged.cov
        merge_blocks -a merged.cov -b regions.bed -o ${PREFIX}.merged.cov
        pigz -p~{threadCount} ${PREFIX}.merged.cov
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File mergedCovGz = glob("*.merged.cov.gz")[0]
    }
}

