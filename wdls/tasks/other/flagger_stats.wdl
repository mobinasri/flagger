version 1.0 

workflow runFlaggerStats{
    call flaggerStats
    output{
        File flaggerStatsTsv = flaggerStats.flaggerStatsTsv
        File flaggerStatsPercOnlyTsv = flaggerStats.flaggerStatsPercOnlyTsv
        File flaggerHapNgxText = flaggerStats.flaggerHapNgxText
    }
}
task flaggerStats {
    input {
        File fastaGz
        File flaggerBed
        File difficultBed_1
        String difficultString_1
        File difficultBed_2
        String difficultString_2
        Array[File]? additionalBeds
        Array[String]? additionalStrings
        File sexBed
        String sample
        String prefix
        Int minContigSize = 1000000
        # runtime configurations
        Int memSize=4
        Int threadCount=4
        Int diskSize=128
        String dockerImage="mobinasri/flagger:v0.4.0--98d66028a969211773077f005511f3d78afdc21c"
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

        ln -s ~{fastaGz} asm.fa.gz
        gunzip -c asm.fa.gz > asm.fa
        # ignore Ns
        python3 /home/scripts/get_contig_coords.py --inputFasta asm.fa | bedtools sort -i - > asm.bed
        cat asm.bed | awk '($3-$2)>~{minContigSize}{print $0}' > asm_long.bed

        bedtools intersect -a asm_long.bed -b ~{difficultBed_1} > diff_long_1.bed
        bedtools intersect -a asm_long.bed -b ~{difficultBed_2} > diff_long_2.bed

        bedtools subtract -a asm.bed -b ~{sexBed} > autosome.bed
        bedtools subtract -a asm.bed -b ~{difficultBed_1} > easy_1.bed
        bedtools subtract -a asm.bed -b ~{difficultBed_2} > easy_2.bed
        bedtools intersect -a easy_1.bed -b easy_2.bed > easy_all.bed

        bedtools intersect -a autosome.bed -b easy_1.bed > autosome_easy_1.bed
        bedtools intersect -a autosome.bed -b easy_2.bed > autosome_easy_2.bed
        bedtools intersect -a autosome.bed -b easy_all.bed > autosome_easy_all.bed

        bedtools intersect -a autosome_easy_all.bed -b asm_long.bed > autosome_easy_all_long.bed

        columns="sample\tinfo"
        values="~{sample}\t~{prefix}"

        columns_2="sample\tinfo"
        values_2="~{sample}\t~{prefix}"

        BED_AND_NAME_TUPLES=()
        BED_AND_NAME_TUPLES[0]="asm.bed,All"
        BED_AND_NAME_TUPLES[1]="asm_long.bed,Long"
        BED_AND_NAME_TUPLES[2]="~{sexBed},sex"
        BED_AND_NAME_TUPLES[3]="autosome.bed,Autosome"
        BED_AND_NAME_TUPLES[4]="~{difficultBed_1},~{difficultString_1}"
        BED_AND_NAME_TUPLES[5]="diff_long_1.bed,~{difficultString_1}_Long"
        BED_AND_NAME_TUPLES[6]="~{difficultBed_2},~{difficultString_2}"
        BED_AND_NAME_TUPLES[7]="diff_long_2.bed,~{difficultString_2}_Long"
        BED_AND_NAME_TUPLES[8]="easy_1.bed,Non_~{difficultString_1}"
        BED_AND_NAME_TUPLES[9]="easy_2.bed,Non_~{difficultString_2}"
        BED_AND_NAME_TUPLES[10]="autosome_easy_1.bed,Autosome_Non_~{difficultString_1}"
        BED_AND_NAME_TUPLES[11]="autosome_easy_2.bed,Autosome_Non_~{difficultString_2}"
        BED_AND_NAME_TUPLES[12]="autosome_easy_all.bed,Autosome_Non_~{difficultString_1}_Non_~{difficultString_2}"
        BED_AND_NAME_TUPLES[13]="autosome_easy_all_long.bed,Autosome_Non_~{difficultString_1}_Non_~{difficultString_2}_Long"

        ADDITIONAL_BED_ARRAY=(~{sep=" " additionalBeds})
        ADDITIONAL_NAME_ARRAY=(~{sep=" " additionalStrings})
        for i in $(seq 0 $(( ${#ADDITIONAL_BED_ARRAY[@]} - 1 )))
        do
            ADDITIONAL_NAME=${ADDITIONAL_NAME_ARRAY[$i]}
            BED_AND_NAME_TUPLES[$i + 14]="${ADDITIONAL_BED_ARRAY[$i]},${ADDITIONAL_NAME}"
        done

        for x in ${BED_AND_NAME_TUPLES[@]}
        do
            IFS=, read bed name <<< "$x"
            err=$(bedtools intersect -a ~{flaggerBed} -b ${bed} | grep "Err" | awk '{s+=$3-$2}END{printf("%.2f", s/1e6)}') || true
            dup=$(bedtools intersect -a ~{flaggerBed} -b ${bed} | grep "Dup" | awk '{s+=$3-$2}END{printf("%.2f", s/1e6)}') || true
            hap=$(bedtools intersect -a ~{flaggerBed} -b ${bed} | grep "Hap" | awk '{s+=$3-$2}END{printf("%.2f", s/1e6)}') || true
            col=$(bedtools intersect -a ~{flaggerBed} -b ${bed} | grep "Col" | awk '{s+=$3-$2}END{printf("%.2f", s/1e6)}') || true
            unk=$(bedtools intersect -a ~{flaggerBed} -b ${bed} | grep "Unk" | awk '{s+=$3-$2}END{printf("%.2f", s/1e6)}') || true
            tot=$(bedtools intersect -a ~{flaggerBed} -b ${bed} | awk '{s+=$3-$2}END{printf("%.2f", s/1e6)}')
            unreliable_perc=$(echo ${err} ${dup} ${hap} ${col} ${unk} ${tot} | awk '{printf "%0.2f", ($1+$2+$4+$5)/($6+1e-10) * 100}')
            values_curr=$(echo ${err} ${dup} ${hap} ${col} ${unk} ${tot} ${unreliable_perc} | awk '{printf $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"($1+$2+$4+$5)"\\t"$7"\\t"$6}')
            values_curr_2="${unreliable_perc}"
            columns_curr="Err_${name}\tDup_${name}\tHap_${name}\tCol_${name}\tUnk_${name}\tUnreliable_${name}\tUnreliable_${name}_Percent\tTotal_${name}"
            columns_curr_2="${name}"
            values="${values}\t${values_curr}"
            columns="${columns}\t${columns_curr}"
            values_2="${values_2}\t${values_curr_2}"
            columns_2="${columns_2}\t${columns_curr_2}"
            # store sorted sizes of the Hap blocks
            printf ${name}"\t" >> ~{sample}.~{prefix}.flagger.hap_ngx.txt
            bedtools intersect -a ~{flaggerBed} -b ${bed} | grep "Hap" | awk '{print $3-$2}' | sort -nk1,1 -r | tr '\n' ',' >> ~{sample}.~{prefix}.flagger.hap_ngx.txt || true
            printf "\n" >> ~{sample}.~{prefix}.flagger.hap_ngx.txt
        done

        printf ${columns}"\n" > ~{sample}.~{prefix}.flagger_stats.tsv
        printf ${values}"\n" >> ~{sample}.~{prefix}.flagger_stats.tsv

        printf ${columns_2}"\n" > ~{sample}.~{prefix}.flagger_stats.perc_only.tsv
        printf ${values_2}"\n" >> ~{sample}.~{prefix}.flagger_stats.perc_only.tsv
 
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File flaggerStatsTsv = glob("*flagger_stats.tsv")[0]
        File flaggerStatsPercOnlyTsv = glob("*perc_only.tsv")[0]
        File flaggerHapNgxText = glob("*flagger.hap_ngx.txt")[0]
    }
}
