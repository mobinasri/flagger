version 1.0

workflow len_filter_fasta {

    call filter_fasta

    output {
        File output_fa_gz = filter_fasta.filteredFasta
    }
}


task filter_fasta {

    input {
        File inputFasta
        String sampleName
        String outputFileTag

        Int min_size   = 100000

        Int memSizeGB  = 4
        Int diskSizeGB = 64
        String dockerImage = "quay.io/biocontainers/seqkit:0.15.0--0"
    }

    String outputFastaFN = "${sampleName}.${outputFileTag}.fa.gz"

    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace


        seqkit seq ~{inputFasta} \
            --min-len ~{min_size} \
            | gzip > ~{outputFastaFN}

    >>>

    output {

        File filteredFasta = outputFastaFN
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}