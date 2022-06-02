version 1.0

import "../../../QC/wdl/tasks/extract_reads.wdl" as extractReads_t
import "tar.wdl" as tar_t
import "gzip.wdl" as gzip_t

workflow filterReads{
    input {
        String sampleName
        Array[File] readsHiFi
        File? referenceFasta
        Int fileExtractionDiskSizeGB = 256
        String zones = "us-west2-a"
    }

    scatter (readFile in readsHiFi) {
        call extractReads_t.extractReads as readsExtracted {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=fileExtractionDiskSizeGB,
                dockerImage="tpesout/hpp_base:latest"
        }
        call filterHiFiAdapter as filterAdapter {
            input:
                readFastq = readsExtracted.extractedRead,
                diskSizeGB = fileExtractionDiskSizeGB
        }
        call gzip_t.gzip as gzip {
            input:
                fileInput = filterAdapter.filteredReadFastq,
                diskSize = fileExtractionDiskSizeGB
        }
    }
    call tar_t.tarGz as blastoutTar{
        input:
            tarGzName = "${sampleName}.HiFiAapterFilt.blastout",
            files = filterAdapter.blastout
    }
    call tar_t.tarGz as blocklistTar{
        input:
            tarGzName = "${sampleName}.HiFiAdapterFilt.blocklist",
            files = filterAdapter.blocklist
    }
    call tar_t.tarGz as countReadsTar{
        input:
            tarGzName = "${sampleName}.countReads",
            files = filterAdapter.countReads
    }
    output {
        Array[File] filteredReadFastqGz = gzip.fileGz
        File adapterBlastOutTarGz = blastoutTar.fileTarGz
        File adapterBlockListTarGz = blocklistTar.fileTarGz
        File countReadsTarGz = countReadsTar.fileTarGz
    }
}


task filterHiFiAdapter {
    input{
        File readFastq
        # runtime configurations
        Int memSizeGB=32
        Int threadCount=16
        Int diskSizeGB=128
        Int preemptible=1
        String dockerImage="quay.io/masri2019/hpp_hifi_adapter_filt:latest"
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

        mkdir data
        cd data
        FILENAME=$(basename -- "~{readFastq}")
        PREFIX="${FILENAME%.*}"
        ln ~{readFastq} ${PREFIX}.fastq
        wc -l ${PREFIX}.fastq | awk '{print $1/4}' > ${PREFIX}.countReads
        bash ${HIFI_ADAPTER_FILTER_BASH} -t ~{threadCount}
        OUTPUTSIZE=`du -s -BG *.filt.fastq | sed 's/G.*//'`
        echo $OUTPUTSIZE > outputsize
    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible: preemptible
    }

    output {
        File filteredReadFastq = glob("data/*.filt.fastq")[0]
        File blastout = glob("data/*.contaminant.blastout")[0]
        File blocklist = glob("data/*.blocklist")[0]
        File countReads = glob("data/*.countReads")[0] 
        Int fileSizeGB = read_int("data/outputsize")
    }
}

task cutadapt {
    input{
        File readFastq
        # runtime configurations
        Int memSizeGB=32
        Int threadCount=16
        Int diskSizeGB=128
        Int preemptible=1
        String dockerImage="quay.io/masri2019/hpp_hifi_adapter_filt:latest"
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

        mkdir data
        cd data
        FILENAME=$(basename -- "~{readFastq}")
        PREFIX="${FILENAME%.*}"
        ln ~{readFastq} ${PREFIX}.fastq
        cutadapt -b "AAAAAAAAAAAAAAAAAATTAACGGAGGAGGAGGA;min_overlap=35" \
         -b "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT;min_overlap=45" \
         --discard-trimmed -o ${PREFIX}.filt.fastq ${PREFIX}.fastq -j ~{threadCount} --revcomp -e 0.05
        OUTPUTSIZE=`du -s -BG *.filt.fastq | sed 's/G.*//'`
        echo $OUTPUTSIZE > outputsize
    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible: preemptible
    }

    output {
        File filteredReadFastq = glob("data/*.filt.fastq")[0]
        Int fileSizeGB = read_int("data/outputsize")
    }
}

