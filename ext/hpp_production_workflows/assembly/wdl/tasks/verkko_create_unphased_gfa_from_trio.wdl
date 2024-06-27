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
        Array[File] input_hifi
        Array[File] input_nanopore
        String name = "assembly"
        String tag  = "verkko_homopolymer_uncompr"

        Int threadCount = 24
        Int memSizeGB   = 180
        Int diskSizeGB  = 3500
        Int preemptible = 0
    }

    command <<<

        set -eux -o pipefail
        ## Neccesary so conda environment will activate...
        source ~/.bashrc

        ## localize nanopore reads to one directory
        ont_files=(~{sep=" " input_nanopore})

        mkdir ont 
        
        for ont_file in ${ont_files[@]};
        do
          cp $ont_file ont/
        done

        ## localize HiFi reads to one directory
        hifi_files=(~{sep=" " input_hifi})
        
        mkdir hifi 

        for hifi_file in ${hifi_files[@]};
        do
            ## Have to copy so I can set timestamp
            cp $hifi_file hifi/
        done


        ## Set modtime to a few years ago so it doesn't trigger snakemake
        touch -a -m -t 202001011205.02 hifi/*
        touch -a -m -t 202001011205.02 ont/*


        ## extract tar w/ snakemake run to cwd
        tar xvf ~{run_folder} --directory ./

        ## Just run consensus w/ paths in --paths using info from prior run in --assembly
        ## This creates an unphased assembly run
        sh $VERKKO/bin/verkko \
            --paths assembly/6-layoutContigs/consensus_paths.txt \
            --assembly assembly \
            -d consensus_unitigs/ \
            --local-memory ~{memSizeGB} \
            --local-cpus ~{threadCount} \
            --hifi hifi/*fast*.gz \
            --nano ont/*fast*.gz


        ## Creates rows like:
        ## unassigned-0000001      utig4-0
        cat consensus_unitigs/6-*/*scfmap \
            | grep utig4 \
            | awk '{print $2"\t"$NF}' \
            > rename.map

        ## Rename contigs in assembly
        python3 $VERKKO/lib/verkko/scripts/process_reads.py rename \
            assembly_renamed.fasta \
            rename.map \
            consensus_unitigs/assembly.fasta
            

        ## Create homopolymer uncompressed version of gfa
        ## Need to use the alignGFA version that is built (the conda version is out of date)
        /opt/verkko/src/build/bin/alignGFA \
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
        docker: "humanpangenomics/verkko_for_gfase@sha256:9c16cd3a2ff3f32b1c87ffe2787701941df434e63271269a10b7a3aded556346"
        preemptible: preemptible
    }
}