version 1.0


workflow verkko_rerun_rukki {

    call rerun_rukki
    
    output{
        File hap1_fa = rerun_rukki.hap1_fa
        File hap2_fa = rerun_rukki.hap2_fa
    }
}

task rerun_rukki {
    
    input {
        File run_folder

        String name = "assembly"
        String tag  = "verkko_gfase"

        Int threadCount = 4
        Int memSizeGB   = 16
        Int diskSizeGB  = 1000
        Int preemptible = 0
    }

    command <<<

        set -eux -o pipefail
        
        ## Neccesary so conda environment will activate...
        source ~/.bashrc


        ## extract tar w/ snakemake run to cwd
        tar xvf ~{run_folder} --directory ./

        cd assembly/6-rukki


        params=""
        params="$params --init-assign new_out_init_ann.csv"
        params="$params --refined-assign new_out_refined_ann.csv"
        params="$params --final-assign new_out_final_ann.csv"
        params="$params --marker-sparsity 5000"
        params="$params --issue-sparsity 1000"
        params="$params --try-fill-bubbles"

        params="$params  --solid-ratio 1.5 --issue-ratio 1. --fillable-bubble-len 500000 --fillable-bubble-diff 1000 --solid-homozygous-cov-coeff 1.1"

        /opt/rukki/target/release/rukki trio -g unitig-popped-unitig-normal-connected-tip.noseq.gfa -m unitig-popped-unitig-normal-connected-tip.colors.csv              -p new-unitig-popped-unitig-normal-connected-tip.paths.tsv $params
        /opt/rukki/target/release/rukki trio -g unitig-popped-unitig-normal-connected-tip.noseq.gfa -m unitig-popped-unitig-normal-connected-tip.colors.csv --gaf-format -p new-unitig-popped-unitig-normal-connected-tip.paths.gaf $params


        grep ISSUE out_final_ann.csv | awk '{print $1}' > issue.ids
        grep -w -f issue.ids new_out_final_ann.csv > to_add.txt
        # utig4-1205      PAT     15258547        m62674:p139679  #8888FF
        # utig4-1206      MAT     15198107        m143101:p65692  #FF8888

        cat to_add.txt | grep "PAT" | cut -f 1 > pat_unitigs.txt
        cat to_add.txt | grep "MAT" | cut -f 1 > mat_unitigs.txt

        ## Now look up contig IDs to pull
        grep -f pat_unitigs.txt ../6-layoutContigs/unitig-popped.layout.scfmap | awk '{print $2}' > pat_contigs_to_add.txt
        grep -f mat_unitigs.txt ../6-layoutContigs/unitig-popped.layout.scfmap | awk '{print $2}' > mat_contigs_to_add.txt

        grep -f pat_contigs_to_add.txt -A1 ../assembly.unassigned.fasta > pat_unassigned.fa
        grep -f mat_contigs_to_add.txt -A1 ../assembly.unassigned.fasta > mat_unassigned.fa


        ## Create a fasta with the original assembly + the unnasigned sequence
        cat ../assembly.haplotype1.fasta mat_unassigned.fa | gzip > ../../assembly.haplotype1.fa.gz
        cat ../assembly.haplotype2.fasta pat_unassigned.fa | gzip > ../../assembly.haplotype2.fa.gz

    >>>

    output {

        File hap1_fa = "assembly.haplotype1.fa.gz"
        File hap2_fa = "assembly.haplotype2.fa.gz"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        cpuPlatform: "Intel Cascade Lake"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "humanpangenomics/verkko_for_gfase@sha256:9beccde8dc9416c34e69f1e9989b6eceddc6f983bbc0459d0fb20f0d95868924"
        preemptible: preemptible
    }
}