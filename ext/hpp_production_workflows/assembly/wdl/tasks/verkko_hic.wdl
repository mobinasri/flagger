version 1.0


workflow run_verkko_hic {

    call verkko_hic
    
    output{
        File hic_phased_run_folder = verkko_hic.hic_phased_run_folder
    }
}

task verkko_hic {
    
    input {
        File run_folder
        Array[File] input_hifi
        Array[File] input_nanopore
        Array[File] aligned_hic      
        String name = "assembly"
        String tag  = "verkko_hic"

        Int threadCount = 72
        Int memSizeGB   = 400
        Int diskSizeGB  = 2500
        Int preemptible = 1
    }

    command <<<

        set -eux -o pipefail
        
        ## Neccesary so conda environment will activate...
        source ~/.bashrc

        ## localize nanopore reads to one directory
        ont_files=(~{sep=" " input_nanopore})

        mkdir -p inputs/ont 
        
        for ont_file in ${ont_files[@]};
        do
          cp $ont_file inputs/ont/
        done

        ## localize HiFi reads to one directory
        hifi_files=(~{sep=" " input_hifi})
        
        mkdir  -p inputs/hifi 

        for hifi_file in ${hifi_files[@]};
        do
            ## Have to copy so I can set timestamp
            cp $hifi_file inputs/hifi/
        done


        ## Set modtime to a few years ago so it doesn't trigger snakemake
        touch -a -m -t 202001011205.02 inputs/hifi/*
        touch -a -m -t 202001011205.02 inputs/ont/*


        ## merge then sort hic bams
        samtools merge hic_merged.bam ~{sep=" " aligned_hic}
        samtools sort -n hic_merged.bam -o hic_merged_sorted.bam
        rm hic_merged.bam
    

        ## extract tar w/ snakemake run to cwd
        tar xvf ~{run_folder} --directory ./

        mkdir verkko_hic

        norealign_wrapper.sh \
            assembly \
            verkko_hic \
            inputs \
            hic_merged_sorted.bam

        tar -cvf ~{name}_~{tag}.tar verkko_hic
    >>>

    output {

        File hic_phased_run_folder = "~{name}_~{tag}.tar"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "humanpangenomics/verkko_hic@sha256:a0842a0cfcac363dd80333dd60c4794914f2bc4ff0dd8bc6ce84a1086dcbbc25"
        preemptible: preemptible
    }
}