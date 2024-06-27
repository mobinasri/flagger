version 1.0


workflow verkko_extract_run_folder {

    call extract_run_folder
    
    output{
        File haplotype1_fa  = extract_run_folder.haplotype1_fa
        File haplotype2_fa  = extract_run_folder.haplotype2_fa
        File unassigned_fa  = extract_run_folder.unassigned_fa
        File diploid_asm_fa = extract_run_folder.diploid_asm_fa
        File mito_fa        = extract_run_folder.mito_fa
        File mito_all_fa    = extract_run_folder.mito_all_fa
    }
}

task extract_run_folder {
    
    input {
        File run_folder
        String sample_name = "assembly"
        String tag         = "verkko_asm"

        Int threadCount = 6
        Int memSizeGB   = 24
        Int diskSizeGB  = 1000
        Int preemptible = 1
    }

    command <<<

        set -eux -o pipefail

        ## extract tar w/ snakemake run to cwd
        tar xvf ~{run_folder} --directory ./


        cp */final_*/assembly.fasta assembly.fasta

        grep "haplotype1.*" assembly.fasta | sed "s/.*>//" > haplotype1.names
        grep "haplotype2.*" assembly.fasta | sed "s/.*>//" > haplotype2.names
        grep "unassigned.*" assembly.fasta | sed "s/.*>//" > unassigned.names

        seqtk subseq assembly.fasta haplotype1.names | gzip > ~{sample_name}_~{tag}_haplotype1.fasta.gz &
        seqtk subseq assembly.fasta haplotype2.names | gzip > ~{sample_name}_~{tag}_haplotype2.fasta.gz &
        seqtk subseq assembly.fasta unassigned.names | gzip > ~{sample_name}_~{tag}_unassigned.fasta.gz &

        gzip -cvf assembly.fasta > ~{sample_name}_~{tag}_diploid.fasta.gz &
        gzip -cvf */final_*/7-consensus/assembly.mito.fasta > ~{sample_name}_~{tag}_mito.fasta.gz &
        gzip -cvf */final_*/7-consensus/assembly.mito.exemplar.fasta > ~{sample_name}_~{tag}_mito_exemplar.fasta.gz &

        wait

    >>>

    output {

        File haplotype1_fa  = "~{sample_name}_~{tag}_haplotype1.fasta.gz"
        File haplotype2_fa  = "~{sample_name}_~{tag}_haplotype2.fasta.gz"
        File unassigned_fa  = "~{sample_name}_~{tag}_unassigned.fasta.gz"
        File diploid_asm_fa = "~{sample_name}_~{tag}_diploid.fasta.gz"
        File mito_fa        = "~{sample_name}_~{tag}_mito_exemplar.fasta.gz"
        File mito_all_fa    = "~{sample_name}_~{tag}_mito.fasta.gz"

    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/biocontainers/seqtk:1.3--hed695b0_2"
        preemptible: preemptible
    }
}
