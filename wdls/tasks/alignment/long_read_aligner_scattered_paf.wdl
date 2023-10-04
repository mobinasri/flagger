version 1.0

import "../../../ext/hpp_production_workflows/QC/wdl/tasks/extract_reads.wdl" as extractReads_t
import "../../../ext/hpp_production_workflows/QC/wdl/tasks/arithmetic.wdl" as arithmetic_t
import "merge_pafs.wdl" as mergePafs_t
import "read_set_splitter.wdl" as readSetSplitter_t
import "long_read_aligner.wdl" as longReadAligner_t


workflow longReadAlignmentSplit {
    input {
        String aligner="winnowmap"
        String preset
        String sampleName
        String sampleSuffix
        Array[File] readFiles
        Int splitNumber = 16
        File assembly
        File? referenceFasta
        Int preemptible=2
        Int extractReadsDiskSize=512
        String zones="us-west2-a"
    }

    scatter (readFile in readFiles) {
        call extractReads_t.extractReads as extractReads {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSize=4,
                threadCount=4,
                diskSizeGB=extractReadsDiskSize,
                dockerImage="tpesout/hpp_base:latest"
        }
    }
    call arithmetic_t.sum as readSize {
        input:
            integers=extractReads.fileSizeGB
    }

    call readSetSplitter_t.readSetSplitter {
        input:
            readFastqs = extractReads.extractedRead,
            splitNumber = splitNumber,
            diskSize = floor(readSize.value * 2.5)
    }
   
    scatter (readFastq in readSetSplitter.splitReadFastqs) {
         ## align reads to the assembly
         call longReadAligner_t.alignmentPaf as alignment{
             input:
                 aligner =  aligner,
                 preset = preset,
                 refAssembly=assembly,
                 readFastq_or_queryAssembly = readFastq,
                 diskSize = 8 + floor(readSize.value / splitNumber * 2.5),
                 preemptible = preemptible,
                 zones = zones
        }
    }
    call arithmetic_t.sum as pafSize {
        input:
            integers=alignment.fileSizeGB
    }

    ## merge the bam files
    call mergePafs_t.merge {
        input:
            sampleName = sampleName,
            sampleSuffix = sampleSuffix,
            pafFiles = alignment.pafFile,
            # runtime configurations
            diskSize = floor(pafSize.value * 2.5),
            preemptible = preemptible,
            zones = zones
    }
    output {
        File pafFile = merge.mergedPaf
    }
}

