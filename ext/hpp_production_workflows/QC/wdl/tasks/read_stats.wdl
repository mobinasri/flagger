version 1.0

import "extract_reads.wdl" as extractReads_t
import "shard_reads.wdl" as shardReads_t
import "arithmetic.wdl" as arithmetic_t

workflow runReadStats {

    input {
        Array[File] reads
        File? referenceFasta
        String identifier="sample"
        Int histogramMinLength=0
        Int histogramMaxLength=0
        Int fileExtractionDiskSizeGB = 256
        String dockerImage = "humanpangenomics/fai_read_stats:latest"

        String identifierAll = "${identifier}_all"
    }

    # extract reads
    scatter (readFile in reads) {
        call extractReads_t.extractReads as readsExtracted {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=fileExtractionDiskSizeGB,
                dockerImage=dockerImage
        }
    }

    # get file size of results
    call arithmetic_t.max as maxReadSize {
        input:
            integers=readsExtracted.fileSizeGB
    }

    # index reads
    scatter (readFile in readsExtracted.extractedRead) {
        call indexReads {
            input:
                readFile=readFile,
                diskSizeGB=maxReadSize.value + 50,
                dockerImage=dockerImage
        }
    }

    # run ReadStats on individual files
    scatter (indexFile in indexReads.indexFile) {
        call readStats {
            input:
                indexFile=indexFile,
                histogramMinLength=histogramMinLength,
                histogramMaxLength=histogramMaxLength,
                dockerImage=dockerImage
        }
    }

    # concatenate fai files (to give per-sample distribution)
    call concatFais {
        input:
            indexFiles=indexReads.indexFile,
            identifier=identifier,

    }

    # run readStats on concatenated fais
    call readStats as concatReadStats {
        input:
            indexFile=concatFais.indexFile,
            histogramMinLength=histogramMinLength,
            histogramMaxLength=histogramMaxLength,
            dockerImage=dockerImage,
        }

    call consolidateReadStats {
        input:
            readStatsTarballs=flatten([readStats.outputTarball, 
                                      [concatReadStats.outputTarball]]),
            readStatsReports=flatten([readStats.reportForConsolidation, 
                                     [concatReadStats.reportForConsolidation]]),
            identifier=identifierAll,

    }

	output {
		File ReadStatsTarball = consolidateReadStats.outputTarball
		File ReadStatsReport = consolidateReadStats.outputReport
        

        ## File ReadStatsTarballAll = consolidateReadStatsAll.outputTarball
        ## File ReadStatsReportAll = consolidateReadStatsAll.outputReport
	}
}


task indexReads {
    input {
        File readFile
        Int memSizeGB = 4
        Int threadCount = 4
        Int diskSizeGB = 64
        String dockerImage = "humanpangenomics/fai_read_stats:latest"
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

        ln -s ~{readFile}
                
        FILE=$(basename ~{readFile})
        OUTPUT="$FILE.fai"

        cat $FILE | awk '{if(NR%4==2) print length($1)}' > $OUTPUT

    >>>

    output {
        File indexFile = glob("*.fai")[0]

    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}


task readStats {
    input {
        File indexFile
        Int histogramMinLength = 0
        Int histogramMaxLength = 0
        Int memSizeGB = 2
        Int threadCount = 2
        Int diskSizeGB = 64
        String dockerImage = "humanpangenomics/fai_read_stats:latest"
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

        # localize file
        ln -s ~{indexFile}
        FILE=$(basename ~{indexFile})

        # get output name
        OUTPUT=$(basename ~{indexFile} | sed -E 's/.(fastq|fq|fasta|fa).fai*$//')

        # hist parameters
        if [[ ~{histogramMinLength} -eq 0 && ~{histogramMaxLength} -eq 0 ]] ; then
            HIST_PARAM="--hist_auto_bounds"
        else
            HIST_PARAM="--hist_min ~{histogramMinLength} --hist_max ~{histogramMaxLength}"
        fi

        # sketch
        fai_read_stats.py -i $FILE -o $OUTPUT $HIST_PARAM

        # for consolidation
        cat $OUTPUT/report.tsv | sed "s/^/$OUTPUT\t/" >$OUTPUT.report.tsv

        # rename and tar
        ls $OUTPUT/ | xargs -n 1 -I{} mv $OUTPUT/{} $OUTPUT/$OUTPUT.{}
        tar czvf $OUTPUT.tar.gz $OUTPUT/

	>>>

	output {
		File outputTarball = glob("*.tar.gz")[0]
		File reportForConsolidation = glob("*.report.tsv")[0]
	}

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}

task consolidateReadStats {
    input {
        Array[File] readStatsTarballs
        Array[File] readStatsReports
        String identifier="sample"
        Int memSizeGB = 4
        Int threadCount = 4
        Int diskSizeGB = 64
        String dockerImage = "humanpangenomics/fai_read_stats:latest"
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
        OMP_NUM_THREADS=~{threadCount}

        # extract all tarballs
        mkdir ~{identifier}_fai_lengths/
        cd ~{identifier}_fai_lengths/
        for tb in ~{sep=" " readStatsTarballs} ; do
            tar xvf $tb &
        done
        wait
        cd ..
        tar czvf ~{identifier}_fai_lengths.tar.gz ~{identifier}_fai_lengths/

        # concat all reports
        for rp in ~{sep=" " readStatsReports} ; do
            cat $rp >>~{identifier}.report.tsv
        done
	>>>

	output {
		File outputTarball = glob("*.tar.gz")[0]
		File outputReport = glob("*.report.tsv")[0]
	}

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}


task concatFais {
    input {
        Array[File] indexFiles
        String identifier="sample"
        Int memSizeGB = 4
        Int threadCount = 1
        Int diskSizeGB = 64
        String dockerImage = "humanpangenomics/fai_read_stats:latest"
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

        # concat all reports
        for fai_file in ~{sep=" " indexFiles} ; do
            cat $fai_file >>~{identifier}.fastq.fai
        done
    >>>

    output {
        File indexFile = glob("~{identifier}.fastq.fai")[0]
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}