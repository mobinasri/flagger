version 1.0


workflow verkko_create_unphased_gfa {

    call create_unphased_gfa
    
    output{
        File unphased_gfa = create_unphased_gfa.unphased_gfa
    }
}

task create_unphased_gfa {
    
    input {
        File run_folder
        String name = "assembly"
        String tag  = "verkko_homopolymer_uncompr"

        Int threadCount = 12
        Int memSizeGB   = 64
        Int diskSizeGB  = 1000
        Int preemptible = 1
    }

    command <<<

        set -eux -o pipefail

        ## Neccesary so conda environment will activate...
        source ~/.bashrc

        ## extract tar w/ snakemake run to cwd
        tar xvf ~{run_folder} --directory ./

        ## Creates rows like:
        ## unassigned-0000001      utig4-0
        cat assembly/6-*/*scfmap \
            | grep utig4 \
            | awk '{print $2"\t"$NF}' \
            > rename.map


        /root/miniconda3/envs/hprc_verkko/lib/verkko/scripts/fasta_combine.py rename \
            assembly_renamed.fasta \
            rename.map \
            assembly/assembly.fasta    


        ## Create homopolymer uncompressed version of gfa
        ## Need to use the alignGFA version that is built (the conda version is out of date)
        /root/miniconda3/envs/hprc_verkko/lib/verkko/bin/alignGFA \
            -V -e 0.30 \
            -gfa \
            -i assembly/assembly.homopolymer-compressed.gfa \
            -T assembly_renamed.fasta 0 \
            -t 8 \
            -o ~{name}_~{tag}.gfa
    >>>

    output {

        File unphased_gfa = "~{name}_~{tag}.gfa"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        cpuPlatform: "Intel Cascade Lake"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "humanpangenomics/hprc_verkko@sha256:7cad26822fd3c1382982cd00fd8f5bb7082909245fd16ee16f222bb4739b54ba"
        preemptible: preemptible
    }
}