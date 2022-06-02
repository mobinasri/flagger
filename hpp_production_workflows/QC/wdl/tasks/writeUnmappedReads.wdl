version 1.0

workflow writeUnmappedReads {

    call writeUnmapped 

    output {
        File unmappedReadsFasta  = writeUnmapped.unmappedReadsFasta
        File unmappedReadsList   = writeUnmapped.unmappedReadsList
    }
}



task writeUnmapped {

    input {
        File inputBam

        Int memSizeGB = 4
        Int diskSizeGB = 32
        String dockerImage = "biocontainers/samtools:v1.9-4-deb_cv1"
    }

    String suffix_stripped_basename  = basename(inputBam, ".bam")
    String outputFasta    = "${suffix_stripped_basename}_unmapped_contigs.fa"
    String outputReadList = "${suffix_stripped_basename}_unmapped_contig_list.txt"

    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace

        ## Pull all unmapped reads
        samtools view \
            -f 4 \
            ~{inputBam} \
            > unmapped_reads.sam

        ## Write read names (contig names) to file
        cat \
            unmapped_reads.sam \
            | cut -f1 \
            > ~{outputReadList}

        ## Write unmapped reads (contigs) as fasta
        samtools fasta unmapped_reads.sam > ~{outputFasta}

    >>>

    output {

        File unmappedReadsFasta  = outputFasta
        File unmappedReadsList   = outputReadList
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}