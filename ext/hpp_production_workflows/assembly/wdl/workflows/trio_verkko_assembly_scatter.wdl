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
        String name = "assembly"
        Int preemptible=1

        File? referenceFasta
        Int fileExtractionDiskSizeGB = 512
    }


    # extract HiFi reads
    scatter (readFile in input_hifi) {
        call extract_reads.extractReads as hifi_reads_extracted {
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
        call extract_reads.extractReads as ont_reads_extracted {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=fileExtractionDiskSizeGB
        }
    }

    call configure_overlap {
        input:
            input_hifi     = hifi_reads_extracted.extractedRead,
            preemptible    = preemptible
    }

    ## compute HiFi overlaps
    scatter(ovlp_job_id in configure_overlap.job_ids) {
        call overlap {
            input:
                configure_ovl_tar = configure_overlap.run_folder_tar,
                job_id            = ovlp_job_id,
                preemptible       = preemptible
        }
    }

    ## Combine overlaps --> graphindexing (plus splitONT)
    call create_graph {
        input:
            run_folder     = configure_overlap.run_folder_tar,
            overlap_tars   = overlap.ovlp_tar,
            input_hifi     = hifi_reads_extracted.extractedRead,            
            input_nanopore = ont_reads_extracted.extractedRead,
            preemptible    = preemptible
    }
    
    ## Align ONT reads
    scatter(aln_job_id in create_graph.job_ids) {
        call graph_aligner {
            input:
                create_graph_tar  = create_graph.run_folder_tar,
                aln_job_id        = aln_job_id,
                preemptible       = preemptible
        }
    }

    ## Combine ONT --> end
    call complete_asm {
        input:
            run_folder     = create_graph.run_folder_tar,
            aln_tars       = graph_aligner.aln_tar,
            input_hifi     = hifi_reads_extracted.extractedRead,            
            input_nanopore = ont_reads_extracted.extractedRead,
            mat_hapmer_tar = mat_hapmer_tar,
            pat_hapmer_tar = pat_hapmer_tar,
            name           = name
    }

    output {
        File mat_fasta  = complete_asm.mat_fasta
        File pat_fasta  = complete_asm.pat_fasta
        File asm_files  = complete_asm.asm_files
        File noseq_gfa  = complete_asm.noseq_gfa
        File colors_csv = complete_asm.colors_csv
    }
}


task configure_overlap {
    input {
        Array[File] input_hifi
        String name     = "assembly"

        Int threadCount = 24
        Int memSizeGB   = 64
        Int diskSizeGB  = 2500
        Int preemptible = 1
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


        mkdir assembly

        ## Run up until, but including, configure overlaps
        verkko \
            -d assembly \
            --local-memory ~{memSizeGB} \
            --local-cpus ~{threadCount} \
            --snakeopts "-U configureOverlaps" \
            --hifi hifi/* 

        ## Run until, but including, countKmers (it's not triggered by previous job)
        verkko \
            -d assembly \
            --local-memory ~{memSizeGB} \
            --local-cpus ~{threadCount} \
            --snakeopts "-U countKmers" \
            --hifi hifi/* 


        ## Get a list of the file names...
        find assembly/0-correction/overlap-jobs/*.params \
            -type f \
            | sed 's/assembly\/0-correction\/overlap-jobs\///' \
            | sed 's/\.params//' \
            > overlap_job_ids.txt


        ## compress all assembly files into an archive
        tar -cvf ~{name}_conf_overlap.tar assembly
    >>>

    output {
        File run_folder_tar   = glob("*_conf_overlap.tar")[0]
        Array[String] job_ids = read_lines("overlap_job_ids.txt")
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        cpuPlatform: "Intel Cascade Lake"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/biocontainers/verkko:1.2--h64afbab_0"
        preemptible: preemptible
    }
}


task overlap {
    input {
        File configure_ovl_tar
        String job_id
        String name = "assembly"

        # String? extra_args

        Int threadCount = 8
        Int memSizeGB   = 24
        Int diskSizeGB  = 250
        Int preemptible = 2
    }

    command <<<

        set -eux -o pipefail

        ## extract tar w/ snakemake run up to now to cwd
        tar xvf ~{configure_ovl_tar} --directory ./

        ## Read the job's parameters to a variable
        job_params=`cat assembly/0-correction/overlap-jobs/~{job_id}.params`
        
        cd assembly/0-correction/overlap-jobs

        ## Compute the overlap...
        /usr/local/lib/verkko/bin/overlapInCore \
            -t 8 \
            -k 28 \
            -k ../../0-correction/hifi.ignoremers \
            --hashbits 25 \
            --hashload 0.8 \
            --maxerate  0.01 \
            --minlength 2000 \
            $job_params \
            -o ../../0-correction/overlap-jobs/~{job_id}.ovb.WORKING \
            -s ../../0-correction/overlap-jobs/~{job_id}.stats \
            ../../0-correction/hifi.seqStore

        ## Mark as done
        mv -f \
            ../../0-correction/overlap-jobs/~{job_id}.ovb.WORKING \
            ../../0-correction/overlap-jobs/~{job_id}.ovb

        cd ../../../

        ## Copy all files from overlap job to a directory and archive
        mkdir ~{job_id}_overlaps

        cp -p \
            assembly/0-correction/overlap-jobs/~{job_id}* \
            ~{job_id}_overlaps/

        tar -cvf ~{job_id}_overlaps.tar ~{job_id}_overlaps

    >>>

    output {
        File ovlp_tar = glob("*_overlaps.tar")[0]
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        cpuPlatform: "Intel Cascade Lake"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/biocontainers/verkko:1.2--h64afbab_0"
        preemptible: preemptible
    }
}

task create_graph {
    input {
        File run_folder
        Array[File] overlap_tars
        Array[File] input_hifi
        Array[File] input_nanopore
        String name = "assembly"

        # String? extra_args

        Int threadCount = 80
        Int memSizeGB   = 400
        Int diskSizeGB  = 2500
        Int preemptible = 1
    }

    command <<<

        set -eux -o pipefail

        ## localize nanopore reads to one directory
        ont_files=(~{sep=" " input_nanopore})

        mkdir ont 
        
        for ont_file in ${ont_files[@]};
        do
          ln -s $ont_file ont/
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

        ## extract tar w/ snakemake run up to configureOverlaps to cwd
        tar xvf ~{run_folder} --directory ./

        ## localize extract overlap files and then copy to snakemake directory
        ovl_tars=(~{sep=" " overlap_tars})

        mkdir tmp

        for ovl_tar in "${ovl_tars[@]}"
        do
            tar xvf $ovl_tar --directory tmp/
        done

        cp -p tmp/*/* assembly/0-correction/overlap-jobs/

        verkko \
            -d assembly \
            --local-memory ~{memSizeGB} \
            --local-cpus ~{threadCount} \
            --red-run 8 32 4 \
            --snakeopts "-U indexGraph" \
            --hifi hifi/* \
            --nano ont/*

        verkko \
            -d assembly \
            --local-memory ~{memSizeGB} \
            --local-cpus ~{threadCount} \
            --red-run 8 32 4 \
            --snakeopts "-U splitONT" \
            --hifi hifi/* \
            --nano ont/*

        ## compress all assembly files into an archive
        tar -cvf ~{name}_indexGraph.tar assembly
    
        
        ## Get a list of the ONT file names to include in each shard...
        find assembly/3-align/split/ont*.fasta.gz \
            -type f \
            | sed 's/assembly\/3-align\/split\/ont//' \
            | sed 's/\.fasta\.gz//' \
            > align_job_ids.txt

    >>>

    output {
        File run_folder_tar    = glob("*_indexGraph.tar")[0]
        Array[String] job_ids = read_lines("align_job_ids.txt")
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        cpuPlatform: "Intel Cascade Lake"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/biocontainers/verkko:1.2--h64afbab_0"
        preemptible: preemptible
    }
}

task graph_aligner {
    input {
        File create_graph_tar
        String aln_job_id

        Int threadCount = 12
        Int memSizeGB   = 32
        Int diskSizeGB  = 1000
        Int preemptible = 1
    }

    command <<<

        set -eux -o pipefail

        ## extract tar w/ snakemake run up to now to cwd
        tar xvf ~{create_graph_tar} --directory ./

        cd assembly/3-align

        ## Align the ONT data...
        /usr/local/bin/GraphAligner \
            -t ~{threadCount} \
            -g ../2-processGraph/unitig-unrolled-hifi-resolved.gfa \
            -f ../3-align/split/ont~{aln_job_id}.fasta.gz \
            -a ../3-align/aligned~{aln_job_id}.WORKING.gaf \
            --seeds-mxm-cache-prefix ../3-align/graph \
            --seeds-mxm-length 30 \
            --seeds-mem-count 10000 \
            --bandwidth 15 \
            --multimap-score-fraction 0.99 \
            --precise-clipping 0.85 \
            --min-alignment-score 5000 \
            --hpc-collapse-reads \
            --discard-cigar \
            --clip-ambiguous-ends 100 \
            --overlap-incompatible-cutoff 0.15 \
            --max-trace-count 5

        mv -f ../3-align/aligned~{aln_job_id}.WORKING.gaf ../3-align/aligned~{aln_job_id}.gaf

        cd ../../


        mkdir ont_~{aln_job_id}_shard
        cp -p assembly/3-align/aligned~{aln_job_id}.gaf ont_~{aln_job_id}_shard/aligned~{aln_job_id}.gaf
        tar -cvf ont_~{aln_job_id}_shard.tar ont_~{aln_job_id}_shard

    >>>

    output {
        File aln_tar = glob("ont_*_shard.tar")[0]
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        cpuPlatform: "Intel Cascade Lake"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/biocontainers/verkko:1.2--h64afbab_0"
        preemptible: preemptible
    }
}

task complete_asm {
    input {
        File run_folder
        Array[File] aln_tars
        Array[File] input_hifi
        Array[File] input_nanopore
        File mat_hapmer_tar
        File pat_hapmer_tar
        String name = "assembly"

        # String? extra_args

        Int threadCount = 64
        Int memSizeGB   = 400
        Int diskSizeGB  = 2500
        Int preemptible = 1
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



        ## extract tar w/ snakemake run up to configureOverlaps to cwd
        tar xvf ~{run_folder} --directory ./

        ## localize extract overlap files and then copy to snakemake directory
        align_tars=(~{sep=" " aln_tars})

        mkdir tmp

        for align_tar in "${align_tars[@]}"
        do
            tar xvf $align_tar --directory tmp/
        done

        cp -p tmp/*/* assembly/3-align/


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
            --local-memory ~{memSizeGB} \
            --local-cpus ~{threadCount} \
            --cns-run 32 0 48 \
            --hifi hifi/* \
            --nano ont/* \
            --hap-kmers  /meryl/mat/$mat_hapmer_db /meryl/pat/$pat_hapmer_db trio


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
        memory: memSizeGB + " GB"
        cpu: threadCount
        cpuPlatform: "Intel Cascade Lake"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/biocontainers/verkko:1.2--h64afbab_0"
        preemptible: preemptible
    }
}