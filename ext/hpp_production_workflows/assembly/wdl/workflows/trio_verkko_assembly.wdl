version 1.0

import "https://raw.githubusercontent.com/human-pangenomics/hpp_production_workflows/master/QC/wdl/tasks/extract_reads.wdl" as extract_reads

workflow verkko_wf {
    
    meta {
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
        description: "Calls Verkko for trio assemblies"
    }
    
    input {
        Array[File] input_hifi
        Array[File] input_nanopore
        File mat_hapmer_tar
        File pat_hapmer_tar

        File? referenceFasta
        Int fileExtractionDiskSizeGB = 256
    }


    # extract HiFi reads
    scatter (readFile in input_hifi) {
        call extract_reads.extractReads as hifi_reads_extracted {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSize=4,
                threadCount=4,
                diskSizeGB=fileExtractionDiskSizeGB
        }
    }

    ## extract nanopore reads
    scatter (readFile in input_nanopore) {
        call extract_reads.extractReads as ont_reads_extracted {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSize=4,
                threadCount=4,
                diskSizeGB=fileExtractionDiskSizeGB
        }
    }

    call verkko {
        input:
            input_hifi     = hifi_reads_extracted.extractedRead,
            input_nanopore = ont_reads_extracted.extractedRead,
            mat_hapmer_tar = mat_hapmer_tar,
            pat_hapmer_tar = pat_hapmer_tar
    }

    output {
        File mat_fasta  = verkko.mat_fasta
        File pat_fasta  = verkko.pat_fasta
        File asm_files  = verkko.asm_files
        File noseq_gfa  = verkko.noseq_gfa
        File colors_csv = verkko.colors_csv
    }
}


task verkko {
    input {
        Array[File] input_hifi
        Array[File] input_nanopore
        File mat_hapmer_tar
        File pat_hapmer_tar
        String name = "assembly"

        String? extra_args

        Int threadCount = 80
        Int memSize   = 240
        Int diskSizeGB  = 2500
    }
    
    parameter_meta {
        input_hifi: "Pacbio HiFi reads for sample. Files must be in fasta or fasta.gz format."
        input_nanopore: "Oxford nanopore reads for sample. Files must be in fasta or fasta.gz format."
        mat_hapmer_tar: "Meryl Hapmer DB as tar file."
        pat_hapmer_tar: "Meryl Hapmer DB as tar file."
        extra_args: "For human assemblies, recommend using: --cns-run 32 0 48 --ali-run 12 32 48 --ovb-run 8 20 48 "
        threadCount: "This WDL requests N2 custom instance types, do not request more than 80 threads."
        memSize: "recommend a ration of 2x-4x of CPU value"
    }


    command <<<

        set -eux -o pipefail

        ## localize HiFi reads to one directory
        hifi_files=(~{sep=" " input_hifi})
        
        mkdir hifi 

        for hifi_file in ${hifi_files[@]};
        do
          ln -s $hifi_file hifi/
        done


        ## localize nanopore reads to one directory
        ont_files=(~{sep=" " input_nanopore})

        mkdir ont 
        
        for ont_file in ${ont_files[@]};
        do
          ln -s $ont_file ont/
        done


        ## Extract meryl hapmer DBs
        mkdir -p /meryl/mat/
        mkdir -p /meryl/pat/

        tar xvf ~{mat_hapmer_tar} --directory /meryl/mat/
        tar xvf ~{pat_hapmer_tar} --directory /meryl/pat/

        mat_hapmer_db="$(basename ~{mat_hapmer_tar} .tar)"
        pat_hapmer_db="$(basename ~{pat_hapmer_tar} .tar)"


        mkdir assembly

        ## Pass extraArgs if set
        if [ -z "~{extra_args}" ]
        then
            EXTRA_ARGS=""
        else
            EXTRA_ARGS="~{extra_args}"
        fi


        ## Run Verkko
        verkko \
            -d assembly \
            --local-memory ~{memSize} \
            --local-cpus ~{threadCount} \
            --hifi hifi/* \
            --nano ont/* \
            --hap-kmers  /meryl/mat/$mat_hapmer_db /meryl/pat/$pat_hapmer_db trio \
            ${EXTRA_ARGS}


        ## compress assemblies for export
        gzip -cvf assembly/assembly.haplotype1.fasta > ~{name}_verkko_v1.2_trio.maternal.fa.gz
        gzip -cvf assembly/assembly.haplotype2.fasta > ~{name}_verkko_v1.2_trio.paternal.fa.gz

        ## compress all assembly files into an archive (useful for reruns locally)
        tar -cvf ~{name}_verkko_v1.2_trio.tar assembly
    >>>

    output {
        File mat_fasta  = glob("*.maternal.fa.gz")[0]
        File pat_fasta  = glob("*.paternal.fa.gz")[0]
        File asm_files  = glob("*_verkko_v1.2_trio.tar")[0]
        File noseq_gfa  = "assembly/assembly.homopolymer-compressed.noseq.gfa"
        File colors_csv = "assembly/6-rukki/unitig-popped-unitig-normal-connected-tip.colors.csv"
    }

    runtime {
        memory: memSize
        cpu: threadCount
        cpuPlatform: "Intel Cascade Lake"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/biocontainers/verkko:1.2--h64afbab_0"
        preemptible: 0
    }
}