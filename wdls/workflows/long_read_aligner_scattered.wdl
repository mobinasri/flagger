version 1.0

import "../../ext/hpp_production_workflows/QC/wdl/tasks/extract_reads_toGZ.wdl" as extractReadsToGZ_t
import "../../ext/hpp_production_workflows/QC/wdl/tasks/arithmetic.wdl" as arithmetic_t
import "../tasks/alignment/merge_bams.wdl" as mergeBams_t
import "../tasks/alignment/read_set_splitter.wdl" as readSetSplitter_t
import "../tasks/alignment/long_read_aligner.wdl" as longReadAligner_t
import "../tasks/alignment/calmd.wdl" as calmd_t
import "../tasks/other/misc.wdl" as misc_t
import "../../ext/secphase/wdls/workflows/correct_bam.wdl" as correct_bam_t
import "../../ext/secphase/wdls/workflows/secphase.wdl" as secphase_t
import "../../ext/hpp_production_workflows/assembly/wdl/tasks/filter_short_reads.wdl" as filter_short_reads_t

workflow longReadAlignmentScattered {
    meta {
        author: "Mobin Asri"
        email: "masri@ucsc.edu"
        description: "Map reads to a diploid or haploid assembly using different mappers like minimap2, winnowmap and veritymap."
    }
    parameter_meta {
        readfiles: "An array of read files. Their format can be either fastq, fq, fastq.gz, fq.gz, bam or cram. For cram format referenceFastaForReadExtraction should also be passed."
        assemblyFasta: "Path to uncompressed or gzip-compressed fasta file of the diploid or haploid assembly" 
        aligner: "Name of the aligner. It can be either minimap2, winnowmap or veritymap. (Default = winnowmap)"
        preset: "Paremeter preset should be selected based on aligner and sequencing platform. Common presets are map-pb/map-hifi/map-ont for minimap2, map-pb/map-ont for winnowmap and hifi-haploid/hifi-haploid-complete/hifi-diploid/ont-haploid-complete for veritymap"
        kmerSize: "The kmer size for using minimap2 or winnowmap. With winnowmap kmer size should be 15 and with minimap2 kmer size should be 17 and 19 for using the presets map-ont and map-hifi/map-pb respectively."
	alignerOptions: "Aligner options. It can be something like '--eqx --cs -Y -L' for minimap2/winnowmap. Note that if assembly is diploid and aligner is either minimap2 or winnowmap '-I8g' is necessary. If the reads contain modification tags and these tags are supposed to be present in the final alignment file, alignerOptions should contain '-y' and the aligner should be either minimap2 or winnowmap. If running secphase is enabled it is recommended to add '-p0.5' to alignerOptions; it will keep more secondary alignments so secphase will have more candidates per read. For veritymap '--careful' can be used but not recommended for whole-genome assembly since it increases the runtime dramatically."
        readExtractionOptions: "The options to be used while converting bam to fastq with samtools fastq. If the reads contain epigentic modification tags it is necessary to use '-TMm,Ml'"
        sampleName: "Name of the sample or assembly. For example 'HG002', 'GRCh38' or 'HG002_hifiasm_v0.19.5'"
        suffix: "Suffix string that contains information about this alignment. It will be appended to the name of the final alignment. For example 'hifi_winnowmap_v2.03_hprc_y2'"
        referenceFastaForReadExtraction: "If reads are in CRAM format then the related reference should be passed to this parameter for read extraction."
        enableAddingMDTag: "If true it will call samtools calmd to add MD tags to the final bam file. [Default = true]"
        splitNumber: "The number of chunks which the input reads should be equally split into. Note that enableSplittingReadsEqually should be set to true if user wants to split reads into equally sized chunks. [Default = 16]"
        enableSplittingReadsEqually: "If true it will merge all reads together and then split them into multiple chunks of roughly equal size. Each chunk will then be aligned via a separate task. This feature is useful for running alignment on cloud/slurm systems where there are  multiple nodes available with enough computing power and having alignments tasks distributed among small nodes is more efficient or cheaper than running a single alignment task in a large node. If the  whole workflow is being on a single node it is not recommened to use this feature since mergin and splitting reads takes its own time. [Default = false]"
        minReadLength: "If it is greater than zero, a task will be executed for filtering reads shorter than this value before alignment. [Default = 0]"
        enableRunningSecphase: "If true it will run Secphase and apply the corrections reported by Secphase to the final output. [Default = false]"
        secphaseOptions: "--hifi for hifi reads and --ont for ont reads [Default = '--hifi'] "
        secphaseDockerImage: "The secphase docker image. [Default = 'mobinasri/secphase:v0.4.3']"
        secphaseVersion: "Secphase version [Default = 'v0.4.3']"
        correctBamOptions: "Options for the correct_bam program that applies the secphase output. [Default = '--primaryOnly --minReadLen 5000 --minAlignment 5000 --maxDiv 0.1']"
        preemptible: "Number of retries to use preemptible nodes on Terra/GCP. [Default = 2]"
        zones: "Name of the zone for taking nodes on Terra/GCP. (Default = us-west2-a)"
    }
    input {
        Array[File] readFiles
        File assemblyFasta
        String aligner="winnowmap"
        String preset
        Int kmerSize = 15
        String alignerOptions="--eqx -Y -L"
        String readExtractionOptions=""
        String sampleName
        String suffix
        File? referenceFastaForReadExtraction
        Boolean enableAddingMDTag=true
        Int splitNumber = 16
        Boolean enableSplittingReadsEqually=false
        Int minReadLength = 0
        Boolean enableRunningSecphase=false
        String secphaseOptions="--hifi"
        String secphaseDockerImage="mobinasri/secphase:v0.4.3"
        String secphaseVersion="v0.4.3"
        String correctBamOptions="--primaryOnly --minReadLen 5000 --minAlignment 5000 --maxDiv 0.1"
        Int preemptible=2
        String zones="us-west2-a"
    }

    # convert reads file into gzip-compressed fastq
    scatter (readFile in readFiles) {
        call extractReadsToGZ_t.extractReadstoGZ {
            input:
                readFile=readFile,
                referenceFasta=referenceFastaForReadExtraction,
                fastqOptions=readExtractionOptions,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB= floor(size(readFile, 'GB')) * 4 + 32,
                dockerImage="mobinasri/bio_base:v0.4.0"
        }
    }
    
    # Distribute reads into multiple files with roughly equal sizes.
    # The number of files is as many as requested by splitNumber parameter
    if (enableSplittingReadsEqually){
        call readSetSplitter_t.readSetSplitter {
            input:
                readFastqs = extractReadstoGZ.extractedRead,
                splitNumber = splitNumber,
                diskSize = floor(size(extractReadstoGZ.extractedRead, 'GB') * 2.5)
        }
    }
    Array[File] readsArrayToMap = select_first([readSetSplitter.splitReadFastqs, extractReadstoGZ.extractedRead])

    # Map reads
    scatter (readFile in readsArrayToMap) {
        # filter short reads if minReadLength is greater than 0
        if (0 < minReadLength){
            call filter_short_reads_t.filterShortReads {
                input:
                    readFastq = readFile,
                    minReadLength = minReadLength,
                    diskSizeGB= floor(size(readFile, 'GB')) * 3 + 32, 
            }
        }
        call longReadAligner_t.alignmentBam as alignment{
            input:
                aligner =  aligner,
                preset = preset,
                kmerSize = kmerSize,
                refAssembly = assemblyFasta,
                readFastq_or_queryAssembly = select_first([filterShortReads.longReadFastqGz, readFile]),
                options = alignerOptions,
                diskSize = 64 + floor(size(readFile, 'GB')) * 6,
                preemptible = preemptible,
                zones = zones,
         }
    }

    ## Merge bam files
    call mergeBams_t.merge as mergeBams{
        input:
            sampleName = sampleName,
            sampleSuffix = suffix,
            sortedBamFiles = alignment.sortedBamFile,
            diskSize = floor(size(alignment.sortedBamFile, 'GB') * 2.5) + 32,
            preemptible = preemptible,
            zones = zones
    }
    
    ## If it is requested add MD tag with samtools calmd
    if (enableAddingMDTag) {
        call calmd_t.calmd {
            input:
                bamFile = mergeBams.mergedBam,
                assemblyFasta = assemblyFasta,
                diskSize = floor(size(mergeBams.mergedBam, 'GB') * 2.5) + 32,
                preemptible = preemptible,
                zones = zones
        }
    }

    # If it is requested run secphase and make a new bam file with corrected alignments
    if (enableRunningSecphase) {
        # gz-compress assembly if it's not since secphase wdl works only with gzipped assembly
        call misc_t.gzipCompress as compressAssembly{
            input:
                inputFile = assemblyFasta
        } 

        # Running Secphase to find the alignments that need correction
        call secphase_t.runSecPhase as secphase {
            input:
                inputBam = select_first([calmd.outputBamFile, mergeBams.mergedBam]),
                diploidAssemblyFastaGz = compressAssembly.gzCompressedFile,
                secphaseOptions = secphaseOptions,
                secphaseDockerImage = secphaseDockerImage,
                version = secphaseVersion
        }

        # Running correctBam to swap primary/secondary tag for detected alignments
        call correct_bam_t.correctBam {
            input:
                phasingLogText = secphase.outLog,
                bam = select_first([calmd.outputBamFile, mergeBams.mergedBam]),
                options = correctBamOptions,
                suffix = "secphase_${secphaseVersion}",
                flagRemoveMultiplePrimary = false,
                flagRemoveSupplementary = false,
                memSize=8,
                threadCount=8,
                dockerImage=secphaseDockerImage,
                preemptible=preemptible
        }
    }


    # Tasks for running secphase and adding MD tags make index files by default
    # if none of them is enabled then make an index using this task
    if ((enableAddingMDTag == false) && (enableRunningSecphase == false)) {
        call longReadAligner_t.indexBam{
            input:
                bam = mergeBams.mergedBam
        }
    }
    File finalBam = select_first([correctBam.correctedBam, calmd.outputBamFile, mergeBams.mergedBam])
    File finalBamIndex = select_first([correctBam.correctedBamIndex, calmd.outputBaiFile, indexBam.bamIndex])

    output {
        File bamFile = finalBam
        File baiFile = finalBamIndex

        # secphase output
        File? secphaseOutputLog = secphase.outLog
        File? secphaseModifiedReadBlocksMarkersBed = secphase.modifiedReadBlocksMarkersBed
        File? secphaseMarkerBlocksBed = secphase.markerBlocksBed
    }
}

