version 1.0 

workflow runFlaggerStats{
    call flaggerStats
    output{
        File flaggerStatsTsv = flaggerStats.flaggerStatsTsv
    }
}
task flaggerStats {
    input {
        File fai
        File flaggerBed
        File difficultBed
        String difficultString
        File sexBed
        String sample
        String prefix
        Int minContigSize = 1000000
        # runtime configurations
        Int memSize=4
        Int threadCount=4
        Int diskSize=128
        String dockerImage="mobinasri/bio_base:v0.1"
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

        cat ~{fai} | awk '{print $1"\t0\t"$2}' | bedtools sort -i - > asm.bed
        cat ~{fai} | awk '$2>~{minContigSize}{print $1"\t0\t"$2}' > asm_long.bed
        bedtools subtract -a asm.bed -b ~{sexBed} > autosome.bed
        bedtools subtract -a asm.bed -b ~{difficultBed} > easy.bed
        bedtools intersect -a autosome.bed -b easy.bed > autosome_easy.bed
        bedtools intersect -a autosome_easy.bed -b asm_long.bed > autosome_easy_long.bed
        

        values="sample\tinfo"
        columns="~{sample}\t~{prefix}"
        for x in asm.bed,All ~{sexBed},sex autosome.bed,Autosome ~{difficultBed},~{difficultString} easy.bed,Non_~{difficultString} autosome_easy.bed,Autosome_Non_~{difficultString} autosome_easy_long,Autosome_Non_~{difficultString}_Long
        do
            IFS=, read bed name <<< "$x"
            err=$(bedtools intersect -a ~{flaggerBed} -b ${bed} | grep "Err" | awk '{s+=$3-$2}END{printf("%.2f", s/1e6)}')
            dup=$(bedtools intersect -a ~{flaggerBed} -b ${bed} | grep "Dup" | awk '{s+=$3-$2}END{printf("%.2f", s/1e6)}')
            hap=$(bedtools intersect -a ~{flaggerBed} -b ${bed} | grep "Hap" | awk '{s+=$3-$2}END{printf("%.2f", s/1e6)}')
            col=$(bedtools intersect -a ~{flaggerBed} -b ${bed} | grep "Col" | awk '{s+=$3-$2}END{printf("%.2f", s/1e6)}')
            tot=$(bedtools intersect -a ~{flaggerBed} -b ${bed} | awk '{s+=$3-$2}END{printf("%.2f", s/1e6)}')
            values_curr=$(echo ${err} ${dup} ${hap} ${col} ${tot} | awk {printf $1"\t"$2"\t"$3"\t"$4"\t"$1+$2+$4"\t"$1+$2+$4/$5 * 100"\t"$5})
            columns_curr="Err_${name}\tDup_${name}\tHap_${name}\tCol_${name}\tUnreliable_${name}\tUnreliable_${name}_Percent\tTotal_${name}"
            values="${values}\t${values_curr}"
            columns="${columns}\t${columns_curr}"
        done

        echo ${columns} > ~{sample}.~{prefix}.flagger_stats.tsv
        echo ${values} >> ~{sample}.~{prefix}.flagger_stats.tsv
 
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File flaggerStatsTsv = glob("*tsv")[0]
    }
}

