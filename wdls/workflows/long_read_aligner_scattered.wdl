version 1.0

import "../../ext/hpp_production_workflows/QC/wdl/tasks/extract_reads.wdl" as extractReads_t
import "../../ext/hpp_production_workflows/QC/wdl/tasks/arithmetic.wdl" as arithmetic_t
import "../tasks/alignment/merge_bams.wdl" as mergeBams_t
import "../tasks/alignment/read_set_splitter.wdl" as readSetSplitter_t
import "../tasks/alignment/long_read_aligner.wdl" as longReadAligner_t
import "../tasks/alignment/calmd.wdl" as calmd_t
import "../tasks/alignment/veritymap.wdl" as verityMap_t

workflow longReadAlignmentScattered {
    input {
        String aligner="winnowmap" # can be either minimap2, winnowmap and veritymap
        String winnowmapOrMinimap2Preset
        String veritymapMode
        String sampleName
        String sampleSuffix
        Array[File] readFiles
        Int splitNumber = 16
        File assembly
        File? referenceFasta
        String winnowmapOrMinimap2Options=""
        String fastqOptions=""
        Int kmerSize = 15
        Int preemptible=2
        Int extractReadsDiskSize=512
        String zones="us-west2-a"
    }

    scatter (readFile in readFiles) {
        call extractReads_t.extractReads as extractReads {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                fastqOptions = fastqOptions,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=extractReadsDiskSize,
                dockerImage="mobinasri/bio_base:dev-v0.1"
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
   
    
   scatter (readFastqAndSize in zip(readSetSplitter.splitReadFastqs, readSetSplitter.readSizes)) {
         ## align reads to the assembly
         if ("${aligner}" == "minimap2" || "${aligner}" == "winnowmap"){
             call longReadAligner_t.alignmentBam as winnowmapOrMinimap2{
                 input:
                     aligner =  aligner,
                     preset = winnowmapOrMinimap2Preset,
                     refAssembly=assembly,
                     readFastq_or_queryAssembly = readFastqAndSize.left,
                     diskSize = 8 + floor(readFastqAndSize.right) * 6,
                     preemptible = preemptible,
                     zones = zones,
                     options = winnowmapOrMinimap2Options,
                     kmerSize = kmerSize
            }
        }
        if ("${aligner}" == "veritymap"){
            call verityMap_t.verityMap as verityMap{
                input:
                    assemblyFastaGz = assembly, 
                    readsFastq = readFastqAndSize.left,
                    suffix = sampleSuffix, 
                    mode = veritymapMode,
                    preemptible = preemptible,
                    zones = zones,
                    diskSize= 8 + floor(readFastqAndSize.right) * 6
            }
        }
    }

    Array[File?] alignmentBams = if ("${aligner}" == "veritymap") then verityMap.sortedBamFile else winnowmapOrMinimap2.sortedBamFile
    Array[Int?] alignmentSizesGB = if ("${aligner}" == "veritymap") then verityMap.fileSizeGB else winnowmapOrMinimap2.fileSizeGB

    call arithmetic_t.sum as bamSize {
        input:
            integers = alignmentSizesGB
    }
    

    ## merge the bam files
    call mergeBams_t.merge as mergeBams{
        input:
            sampleName = sampleName,
            sampleSuffix = sampleSuffix,
            sortedBamFiles = alignmentBams,
            # runtime configurations
            diskSize = floor(bamSize.value * 2.5) + 32,
            preemptible = preemptible,
            zones = zones
    }
    
    ## add Md tag
    call calmd_t.calmd {
        input:
            bamFile = mergeBams.mergedBam,
            assemblyFastaGz = assembly,
            diskSize = floor(bamSize.value * 2.5) + 32,
            preemptible = preemptible,
            zones = zones
    }
    output {
        File bamFile = calmd.outputBamFile
        File baiFile = calmd.outputBaiFile
    }
}

