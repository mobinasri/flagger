version 1.0

import "../../ext/hpp_production_workflows/QC/wdl/tasks/extract_reads.wdl" as extractReads_t
import "../../ext/hpp_production_workflows/QC/wdl/tasks/arithmetic.wdl" as arithmetic_t
import "../tasks/alignment/long_read_aligner.wdl" as longReadAligner_t
import "../tasks/alignment/merge_bams.wdl" as mergeBams_t

workflow AlignSubsetReads {
    input {
        String aligner="winnowmap"
        String preset
        String sampleName
        String sampleSuffix
        Array[File] readFiles
        File assembly
        File readListTarGz
        File? referenceFasta
        Int preemptible=2
        Int extractReadsDiskSize=256
        String zones
    }

    scatter (readFile in readFiles) {
        call extractReads_t.extractReads as extractReads {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=extractReadsDiskSize,
                dockerImage="tpesout/hpp_base:latest"
        }
        
        call subset {
           input:
                fastq = extractReads.extractedRead,
                readListTarGz = readListTarGz,
                diskSize = 32 + extractReads.fileSizeGB,
         } 
         ## align reads to the assembly
         call longReadAligner_t.alignmentBam as alignment{
             input:
                 aligner =  aligner,
                 preset = preset,
                 refAssembly=assembly,
                 readFastq_or_queryAssembly = subset.subsetReads,
                 diskSize = extractReads.fileSizeGB * 3,
                 threadCount = 8,
                 memSize = 16,
                 preemptible = preemptible,
                 zones = zones
        }
    }

    call arithmetic_t.sum as bamSize {
        input:
            integers=alignment.fileSizeGB
    }

    ## merge the bam files
    call mergeBams_t.merge as mergeBams{
        input:
            sampleName = sampleName,
            sampleSuffix = sampleSuffix,
            sortedBamFiles = alignment.sortedBamFile,
            # runtime configurations
            diskSize = floor(bamSize.value * 2.5),
            preemptible = preemptible,
            zones = zones
    }
    output {
        File sortedBamFile = mergeBams.mergedBam
        File baiFile = mergeBams.mergedBai
    }

}

task subset {
    input {
        File fastq
        File readListTarGz
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=32
        String dockerImage="mobinasri/flagger:v0.3.0"
        Int preemptible=2
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

        mkdir read_lists
        tar -xzvf ~{readListTarGz} --strip-components 1 --directory read_lists
        cat read_lists/* > whole_list.txt

        FILENAME=$(basename ~{fastq})
        mkdir output
        # -A3 is for grepping sequences and base qualities
        grep -F -f whole_list.txt ~{fastq} -A3 | grep -v "^\-\-$" > output/${FILENAME}

    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File subsetReads = glob("output/*")[0]
    }
}

