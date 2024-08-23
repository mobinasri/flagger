version 1.0

workflow run_calc_ont_stats {
    meta {
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
        description: "Calculates summary statistics for ONT runs from sequencing_summary files."
    }

    parameter_meta {
        sequencing_summary: "Sequencing summary file output by Dorado/Guppy."
    }


    input {
        File sequencing_summary
    }


    call calc_ont_summary_stats {
        input:
            sequencing_summary  = sequencing_summary
    }


   output {
        File pass_summary_stats = calc_ont_summary_stats.pass_output_stats
        File fail_summary_stats = calc_ont_summary_stats.fail_output_stats
    }  
}


task calc_ont_summary_stats {
    input {

        File sequencing_summary

        Int memSizeGB   = 8
        Int threadCount = 2
        Int diskSize    = 50       
        Int preempts    = 3
    }

    command <<<
        set -eux -o pipefail

        BASENAME=$(basename -- "~{sequencing_summary}" _sequencing_summary.txt.gz)
        
        # Create summary for just pass reads
        # grab header, then reads over Q10
        #gzip -cd ~{sequencing_summary} | sed -n 1p > ${BASENAME}_pass.txt
        zcat ~{sequencing_summary} | awk '{ if ($15 >= 10) { print } }' >> ${BASENAME}_pass.txt

        # Create summary for just fail reads
        gzip -cd ~{sequencing_summary} | sed -n 1p > ${BASENAME}_fail.txt
        zcat ~{sequencing_summary} | awk '{ if ($15 < 10) { print } }' >> ${BASENAME}_fail.txt

        python3 /opt/calculate_summary_stats.py ${BASENAME}_pass.txt > ${BASENAME}.pass_summary_stats.txt
        python3 /opt/calculate_summary_stats.py ${BASENAME}_fail.txt > ${BASENAME}.fail_summary_stats.txt

    >>>

    output {
        File pass_output_stats = glob("*.pass_summary_stats.txt")[0]
        File fail_output_stats = glob("*.fail_summary_stats.txt")[0]
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        docker: "humanpangenomics/ont_summary_stats@sha256:ec6eba759f8f91cf461262ba5155aa951b77acf241f764dda4b58f7d5247d96b"
        preemptible: preempts
    }
}
