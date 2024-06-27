version 1.0

workflow runExtractReadstoGZ {
    input {
        Array[File] inputFiles
        File? referenceFasta
        Int threadCount=1
        String dockerImage
    }

    scatter (file in inputFiles) {
        call extractReadstoGZ {
            input:
                readFile=file,
                referenceFasta=referenceFasta,
                threadCount=threadCount,
                dockerImage=dockerImage
        }
    }

    output {
        Array[File] reads = extractReadstoGZ.extractedRead
    }
}

task extractReadstoGZ {
    input {
        File readFile
        File? referenceFasta
        Boolean sortByName = false
        String excludeString="" # exclude lines with this string from fastq
        String fastqOptions = ""
        Int memSizeGB = 4
        Int threadCount = 8
        Int diskSizeGB = 128
        String dockerImage = "mobinasri/bio_base:v0.3.0"
    }


    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        FILENAME=$(basename -- "~{readFile}")
        PREFIX="${FILENAME%.*}"
        SUFFIX="${FILENAME##*.}"

        
        mkdir output
        SORT_BY_NAME=~{true="yes" false="no" sortByName}

        if [[ "$SUFFIX" == "bam" ]] ; then
            if [[ ${SORT_BY_NAME} == "yes" ]]; then
                samtools sort -@~{threadCount} -n ~{readFile} | \
                         samtools fastq ~{fastqOptions} -@~{threadCount} - | \
                         pigz -p~{threadCount} > output/${PREFIX}.fq.gz
            else
                samtools fastq ~{fastqOptions} -@~{threadCount} ~{readFile} | \
                         pigz -p~{threadCount} > output/${PREFIX}.fq.gz
            fi
        elif [[ "$SUFFIX" == "cram" ]] ; then
            if [[ ! -f "~{referenceFasta}" ]] ; then
                echo "Could not extract $FILENAME, reference file not supplied"
                exit 1
            fi
            ln -s ~{referenceFasta}
            if [[ ${SORT_BY_NAME} == "yes" ]]; then
                samtools sort -@~{threadCount} -n --reference `basename ~{referenceFasta}` -O bam ~{readFile} | \
                         samtools fastq ~{fastqOptions} -@~{threadCount} - | \
                         pigz -p~{threadCount} > output/${PREFIX}.fq.gz 
            else
                samtools fastq ~{fastqOptions} -@~{threadCount} --reference `basename ~{referenceFasta}` ~{readFile} | \
                         pigz -p~{threadCount} > output/${PREFIX}.fq.gz
            fi
        elif [[ "$SUFFIX" == "gz" ]] ; then
            ln -s ~{readFile} output/${PREFIX}.gz
        elif [[ "$SUFFIX" == "fastq" ]] || [[ "$SUFFIX" == "fq" ]] ; then
            cat ~{readFile} | pigz -p~{threadCount} > output/${PREFIX}.fq.gz
        elif [[ "$SUFFIX" != "fastq" ]] && [[ "$SUFFIX" != "fq" ]] && [[ "$SUFFIX" != "fasta" ]] && [[ "$SUFFIX" != "fa" ]] ; then
            echo "Unsupported file type: ${SUFFIX}"
            exit 1
        fi


        mkdir output_final
        OUTPUT_NAME=$(ls output)

        if [ "~{excludeString}" != "" ]; then
            zcat output/${OUTPUT_NAME} | grep -v "~{excludeString}" | pigz -p~{threadCount} > output_final/${OUTPUT_NAME}
        else
            ln output/${OUTPUT_NAME} output_final/${OUTPUT_NAME}
        fi
        

        OUTPUTSIZE=`du -s -BG output_final/ | sed 's/G.*//'`
        if [[ "0" == $OUTPUTSIZE ]] ; then
            OUTPUTSIZE=`du -s -BG ~{readFile} | sed 's/G.*//'`
        fi
        echo $OUTPUTSIZE >outputsize
    >>>

    output {
        File extractedRead = flatten([glob("output_final/*"), [readFile]])[0]
        Int fileSizeGB = read_int("outputsize")
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }

    parameter_meta {
        readFile: {description: "Reads file in BAM, FASTQ, or FASTA format (optionally gzipped)"}
    }
}
