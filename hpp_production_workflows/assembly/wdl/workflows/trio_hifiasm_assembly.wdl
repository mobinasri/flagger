version 1.0 

import "../tasks/yak_no_stats.wdl" as yak_t
import "../tasks/hifiasm.wdl" as hifiasm_t
import "../tasks/gfatools.wdl" as gfatools_t

workflow trioHifiasmAssembly {
    input {
        String childID
        String paternalID
        String maternalID
        Array[File] childReadsHiFi
        Array[File] paternalReadsILM
        Array[File] maternalReadsILM
        String? hifiasmExtraOptions
        File? inputBinFilesTarGz
        File? referenceFasta
        # runtime configurations for running hifiasm
        Int threadCount=48
        Int memSizeGB=256
        Int preemptible=1
    }


    ### Yak ###
    call yak_t.runYak as paternalYakCount {
        input:
            sampleReadsILM=paternalReadsILM,
            sampleName="pat.${paternalID}",
            referenceFasta=referenceFasta,
    }
    call yak_t.runYak as maternalYakCount {
        input:
            sampleReadsILM=maternalReadsILM,
            sampleName="mat.${maternalID}",
            referenceFasta=referenceFasta,
    }

    ### Trio Hifiasm ###
    call hifiasm_t.runTrioHifiasm as trioHifiasm{
        input:
            paternalYak = paternalYakCount.outputYak,
            maternalYak = maternalYakCount.outputYak,
            childReadsHiFi = childReadsHiFi,
            childID = childID,
            hifiasmExtraOptions = hifiasmExtraOptions,
            inputBinFilesTarGz = inputBinFilesTarGz,
            memSizeGB = memSizeGB,
            threadCount = threadCount,
            preemptible = preemptible
    }

    ### Convert GFA to FASTA ###
    call gfatools_t.phasedGFAs2Fasta as gfa2fasta{
        input:
            paternalGfa = trioHifiasm.outputPaternalGfa,
            maternalGfa = trioHifiasm.outputMaternalGfa,
            childID = childID
    }

    output {
        File paternalFastaGz = gfa2fasta.outputPaternalFastaGz
        File maternalFastaGz = gfa2fasta.outputMaternalFastaGz
        File paternalYak = paternalYakCount.outputYak
        File maternalYak = maternalYakCount.outputYak
        File paternalContigGfaTarGz = trioHifiasm.outputPaternalContigGfa 
        File maternalContigGfaTarGz = trioHifiasm.outputMaternalContigGfa 
        File rawUnitigGfaTarGz = trioHifiasm.outputRawUnitigGfa
        File binFilesTarGz = trioHifiasm.outputBinFiles
    }
    parameter_meta {
        childID: " NIST ID (or Coriell ID) of the child sample whose reads are going to be assembled"
        paternalID: " NIST ID (or Coriell ID) of the paternal sample"
        maternalID: "NIST ID (or Coriell ID) of the maternal sample"
        childReadsHiFi: "An array of files (or a single file) that contain the HiFi reads of the child sample ( Acceptable formats are fastq (or fq), fastq.gz (or fq.gz), bam and cram)"
        paternalReadsILM: "An array of files (or a single file) that contain the Illumina short reads of the child sample ( Acceptable formats are fastq (or fq), fastq.gz (or fq.gz), bam and cram)"
        maternalReadsILM: "An array of files (or a single file) that contain the Illumina short reads of the child sample ( Acceptable formats are fastq (or fq), fastq.gz (or fq.gz), bam and cram)"
        inputBinFilesTarGz: "(optional) The hifiasm produces some bin files which can be saved and used for re-running the assembly process. By having these bin files hifiasm can skip the time-consuming process of finding overlaps (Acceptable format is .tar.gz)"
        referenceFasta: "(optional) If any of the read files (can be either for child, father or mother) are having .cram format, the reference genome should be provided in .fasta format"
        threadCount: "(default=48) The number of cores for running hifiasm"
        memSize: "(default=256) The memory size (GB) for running hifiasm"
        preemptible: "(default=1) The number of tries for using a preemptible node for running hifiasm. Note that if your child data has a coverage of more than 40X, hifiasm (without any given bin files) may take longer than 24 hours. So using a preemptible node is useless beacuse it gets interrupted after 24 hours"
    }
    meta {
        author: "Mobin Asri"
        email: "masri@ucsc.edu"
    }
}

