version 1.0

workflow ntsm_workflow {
    input {
        Array[File] input_reads_1
        Array[File] input_reads_2
        String sample_id   = "sample"
        String read_1_type = "type1"
        String read_2_type = "type2"

        File? cram_reference
    }

    # Type 1: NTSM count
    call ntsm_count as ntsm_count_1 {
        input:
            input_reads    = input_reads_1,
            cram_reference = cram_reference,
            sample_id      = sample_id,
            read_type      = read_1_type
    }

    # Type 2: NTSM count
    call ntsm_count as ntsm_count_2 {
        input:
            input_reads    = input_reads_2,
            cram_reference = cram_reference,
            sample_id      = sample_id,
            read_type      = read_2_type
    }

    call ntsm_eval {
        input:
            count_1_file = ntsm_count_1.ntsm_counts,
            count_2_file = ntsm_count_2.ntsm_counts,
            sample_id    = sample_id,
            read_1_type  = read_1_type,
            read_2_type  = read_2_type
    }

   output {
        File ntsv_count_1  = ntsm_count_1.ntsm_counts
        File ntsv_count_2  = ntsm_count_2.ntsm_counts
        File ntsm_eval_out = ntsm_eval.ntsm_eval_out
    }

    meta {
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
        description: "Calls [ntsm](https://github.com/JustinChu/ntsm) to identify mismatched read files through k-mer analysis"
    }
}

task ntsm_count {
    input {
        Array[File] input_reads
        String sample_id
        String read_type

        File? cram_reference
        Int count_kmer_size = 19

        Int memSizeGB   = 4
        Int threadCount = 4
        Int addldisk    = 10
        Int preempts    = 2
    }
    
    parameter_meta {
        input_reads: "Files must be in fastq, fastq.gz (preferred), BAM, or CRAM format. CRAM files must be passed with cram_reference"
        count_kmer_size: "k-mer size to use (sliding window is applied: highest count is stored)"
    }

    # Estimate disk size required
    Int input_read_size  = ceil(size(input_reads, "GB"))       
    Int final_disk_dize   = input_read_size * 5 + addldisk

    # Create output file name
    String output_counts_fn = "${sample_id}_${read_type}_counts.txt"

    command <<<
        set -eux -o pipefail

        READS=(~{sep=" " input_reads})

        ## dump converted reads into output folder
        mkdir output

        ## convert all reads to fastq.gz format
        for READFILE in "${READS[@]}"
            do
                echo "$READFILE"

                FILENAME="$(basename -- $READFILE)"
                PREFIX="${FILENAME%.*}"
                SUFFIX="${FILENAME##*.}"

                if [[ "$SUFFIX" == "bam" ]] ; then
                    samtools fastq -@~{threadCount} $READFILE | gzip --fast --stdout > output/${PREFIX}.fastq.gz

                elif [[ "$SUFFIX" == "cram" ]] ; then
                    if [[ ! -f "~{cram_reference}" ]] ; then
                        echo "Could not extract $FILENAME, reference file not supplied"
                        exit 1
                    fi
                    ln -s ~{cram_reference}
                    # samtools fastq -@~{threadCount} -c 1 -1 output/${PREFIX}_1.fastq.gz -2 output/${PREFIX}_2.fastq.gz --reference `basename ~{cram_reference}` $READFILE
                    samtools fastq -@~{threadCount} --reference `basename ~{cram_reference}` $READFILE | gzip --fast --stdout >  output/${PREFIX}.fastq.gz              
                elif [[ "$SUFFIX" == "fastq" ]] ; then
                    gzip $READFILE > output/${PREFIX}.fastq.gz
                
                elif [[ "$SUFFIX" != "gz" ]] ; then
                    echo "Unsupported file type: ${PREFIX}.${SUFFIX}"
                    exit 1
                else
                    ln $READFILE output/$FILENAME
                fi

            done


        ## localize kmers
        ln -s /opt/ntsm/data/31_AT.fa .
        ln -s /opt/ntsm/data/31_CG.fa .

        ## Count number of k-mers found 
        /opt/ntsm/src/ntsmCount \
            -r 31_AT.fa \
            -a 31_CG.fa \
            -k ~{count_kmer_size} \
            -t ~{threadCount} \
            <(pigz -cd output/*) \
            > ~{output_counts_fn}
    >>>

    output {
        File ntsm_counts = output_counts_fn
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + final_disk_dize + " SSD"
        docker: "humanpangenomics/ntsm:latest"
        preemptible: preempts
    }
}

task ntsm_eval {
    input {
        File count_1_file
        File count_2_file       
        String sample_id
        String read_1_type
        String read_2_type

        Int memSizeGB   = 4
        Int threadCount = 4
        Int diskSize    = 32
        Int preempts    = 2
    }

    # Create output file name
    String output_eval_fn = "${sample_id}_${read_1_type}_vs_${read_2_type}.txt"

    command <<<
        set -eux -o pipefail

        /opt/ntsm/src/ntsmEval \
            ~{count_1_file} \
            ~{count_2_file} \
            > ~{output_eval_fn}
    >>>

    output {
        File ntsm_eval_out = output_eval_fn
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        docker: "humanpangenomics/ntsm:latest"
        preemptible: preempts
    }
}