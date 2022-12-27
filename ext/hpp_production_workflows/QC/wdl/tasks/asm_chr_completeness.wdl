version 1.0

workflow run_calc_completeness {
    input {
        File ref_fasta
        File asm_fasta
        String name
    }

    call calc_completeness {
        input:
            ref_fasta  = ref_fasta,
            asm_fasta  = asm_fasta,
            name       = name
    }

   output {
        File chr_translation      = calc_completeness.chr_translation
        File chr_completeness     = calc_completeness.chr_completeness
        File chr_completeness_max = calc_completeness.chr_completeness_max
    }  
}


task calc_completeness {
    input {
        File ref_fasta
        File asm_fasta
        String name

        Int memSizeGB   = 24
        Int threadCount = 16
        Int diskSize    = 100
        Int preempts    = 2
    }


    command <<<
        set -eux -o pipefail

        mashmap \
            -r ~{ref_fasta} \
            -q ~{asm_fasta} \
            -f one-to-one \
            --pi 95 \
            -s 10000


        cat mashmap.out \
            | awk '{if ($NF > 99 && $4-$3 > 500000) print $1"\t"$6"\t"$2"\t"$7}' \
            | sort | uniq \
            > ~{name}.translation.txt

        cat ~{name}.translation.txt \
            | sort -k2,2 \
            | awk '{if ($3 > 15000000) print $0}' \
            | awk -v LAST="" -v S="" '{if (LAST != $2) { if (S > 0) print LAST"\t"C"\t"SUM/S*100"\t"TIG; SUM=0; C=0; } LAST=$2; S=$NF; SUM+=$3; C+=1; TIG=$1} END \
            {print LAST"\t"C"\t"SUM/S*100"\t"TIG;}' \
            | awk '{if ($2 == 1 && $(NF-1) > 95) print $0}' \
            | sort -k4 \
            > ~{name}.chr_completeness.txt

        cat ~{name}.translation.txt \
            | sort -k2,2 \
            | awk '{if ($3 > 15000000) print $0}' \
            | awk -v LAST="" -v S="" '{if (LAST != $2) { if (S > 0) print LAST"\t"C"\t"SUM/S*100"\t"MAX/S*100"\t"TIG; SUM=0; MAX=0; C=0; } LAST=$2; S=$NF; SUM+=$3; if (MAX < $3) MAX=$3; C+=1; TIG=$1} END {print LAST"\t"C"\t"SUM/S*100"\t"MAX/S*100"\t"TIG;}' | awk '{print $1"\t"$4}' \
            | sort -nk1,1 -s > ~{name}.chr_completeness_max.txt
   
    >>>

    output {
        File chr_translation = "~{name}.translation.txt"
        File chr_completeness = "~{name}.chr_completeness.txt"
        File chr_completeness_max = "~{name}.chr_completeness_max.txt"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        docker: "quay.io/biocontainers/mashmap:2.0--h543ed2d_4"
        preemptible: preempts
    }
}
