version 1.0

import "https://raw.githubusercontent.com/human-pangenomics/hpp_production_workflows/master/QC/wdl/tasks/extract_reads_toGZ.wdl" as extractReadstoGZ

workflow verkko_wf {
    
    meta {
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
        description: "Calls Verkko for trio assemblies"
    }
    
    input {
        Array[File] input_hifi
        Array[File] input_nanopore
        File unphased_asm_tar
        File mat_hapmer_tar
        File pat_hapmer_tar
        String name = "assembly"

        File? referenceFasta
        Int fileExtractionDiskSizeGB = 512
    }


    # extract HiFi reads
    scatter (readFile in input_hifi) {
        call extractReadstoGZ.extractReadstoGZ as hifi_reads_extracted {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=fileExtractionDiskSizeGB
        }
    }

    ## extract nanopore reads
    scatter (readFile in input_nanopore) {
        call extractReadstoGZ.extractReadstoGZ as ont_reads_extracted {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=fileExtractionDiskSizeGB
        }
    }


    ## Trio phasing
    call trio_phase {
        input:
            run_folder     = unphased_asm_tar,
            input_hifi     = hifi_reads_extracted.extractedRead,            
            input_nanopore = ont_reads_extracted.extractedRead,
            mat_hapmer_tar = mat_hapmer_tar,
            pat_hapmer_tar = pat_hapmer_tar,
            name           = name
    }

    output {
        File mat_fasta  = trio_phase.mat_fasta
        File pat_fasta  = trio_phase.pat_fasta
        File asm_files  = trio_phase.asm_files
        File noseq_gfa  = trio_phase.noseq_gfa
        File colors_csv = trio_phase.colors_csv
    }
}



task trio_phase {
    input {
        File run_folder
        Array[File] input_hifi
        Array[File] input_nanopore
        File mat_hapmer_tar
        File pat_hapmer_tar
        String name = "assembly"

        Int threadCount = 64
        Int memSizeGB   = 400
        Int diskSizeGB  = 2000
        Int preemptible = 0
    }

    command <<<

        set -eux -o pipefail

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


        ## extract tar w/ unphased verkko run
        tar xvf ~{run_folder} --directory ./

        ## clean up to trigger full rerun
        rm assembly/assembly.*
        rm -rf assembly/6-*
        mkdir tmp 
        cp -p assembly/7-consensus/ont* tmp/
        rm -rf assembly/7-consensus/*
        cp -p tmp/* assembly/7-consensus/

        ## Extract meryl hapmer DBs
        mkdir -p /meryl/mat/
        mkdir -p /meryl/pat/

        tar xvf ~{mat_hapmer_tar} --directory /meryl/mat/
        tar xvf ~{pat_hapmer_tar} --directory /meryl/pat/

        mat_hapmer_db="$(basename ~{mat_hapmer_tar} .tar)"
        pat_hapmer_db="$(basename ~{pat_hapmer_tar} .tar)"


        ## Call Verkko 
        verkko \
            -d assembly \
            --screen human \
            --local-memory ~{memSizeGB} \
            --local-cpus ~{threadCount} \
            --hifi hifi/* \
            --nano ont/* \
            --hap-kmers  /meryl/mat/$mat_hapmer_db /meryl/pat/$pat_hapmer_db trio \
            --snakeopts "-R rukki"


        ## compress assemblies for export
        gzip -cvf assembly/assembly.haplotype1.fasta > ~{name}_verkko_v1.3_trio.maternal.fa.gz
        gzip -cvf assembly/assembly.haplotype2.fasta > ~{name}_verkko_v1.3_trio.paternal.fa.gz

        ## compress all assembly files into an archive (useful for reruns locally)
        tar -cvf ~{name}_verkko_v1.3_trio.tar assembly

    >>>

    output {
        File mat_fasta  = glob("*.maternal.fa.gz")[0]
        File pat_fasta  = glob("*.paternal.fa.gz")[0]
        File asm_files  = glob("*_verkko_v1.3_trio.tar")[0]
        File noseq_gfa  = "assembly/assembly.homopolymer-compressed.noseq.gfa"
        File colors_csv = "assembly/6-rukki/unitig-popped-unitig-normal-connected-tip.colors.csv"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        cpuPlatform: "Intel Cascade Lake"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/biocontainers/verkko:1.3.1--h64afbab_0"
        preemptible: preemptible
    }
}