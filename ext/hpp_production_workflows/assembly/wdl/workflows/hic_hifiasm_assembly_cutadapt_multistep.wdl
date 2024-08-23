version 1.0

import "../tasks/hifiasm_cutadapt_hic_multistep.wdl" as hifiasm_multistep_wf
import "../tasks/gfatools.wdl" as gfatools_wf
import "../tasks/hic_group_xy.wdl" as hic_group_xy_wf

workflow hicHifiasmAssembly {
    input {
        String childID
        Boolean isMaleSample
        Array[File] childReadsHiFi
        Array[File] childReadsHiC1
        Array[File] childReadsHiC2
        Array[File] childReadsONT=[]
        Int? minOntReadLength
        Int? homCov
        Boolean filterAdapters=true
        String? hifiasmExtraOptions
        File? inputBinFilesTarGz
        File? referenceFasta
        File chrY_no_par_yak
        File chrX_no_par_yak
        File par_yak
        # runtime configurations for running hifiasm
        Int threadCount=64
        Int preemptible=2
    }

    ### HiC Hifiasm ###
    call hifiasm_multistep_wf.runHiCHifiasm as hicHifiasm{
        input:
            childReadsHiC1 = childReadsHiC1,
            childReadsHiC2 = childReadsHiC2,
            childReadsHiFi = childReadsHiFi,
            childReadsONT = childReadsONT,
            minOntReadLength = minOntReadLength,
            homCov = homCov,
            childID = childID,
            filterAdapters = filterAdapters,
            hifiasmExtraOptions = hifiasmExtraOptions,
            inputBinFilesTarGz = inputBinFilesTarGz,
            referenceFasta = referenceFasta,
            threadCount = threadCount,
            preemptible = preemptible
    }

    ### Convert GFA to FASTA ###
    call gfatools_wf.phasedGFAs2Fasta as gfa2fasta{
        input:
            paternalGfa = hicHifiasm.outputHaplotype1Gfa,
            maternalGfa = hicHifiasm.outputHaplotype2Gfa,
            patSuffix   = "hap1",
            matSuffix   = "hap2",
            childID     = childID
    }

    ### Ensure that X/Y chromosomes are properly partitioned
    call hic_group_xy_wf.group_xy as asm_group_xy{
        input:
            hap1_gz         = gfa2fasta.outputPaternalFastaGz,
            hap2_gz         = gfa2fasta.outputMaternalFastaGz,
            childID         = childID,
            isMaleSample    = isMaleSample,
            chrY_no_par_yak = chrY_no_par_yak,
            chrX_no_par_yak = chrX_no_par_yak,
            par_yak         = par_yak
    }

    output {
        File rawhap1FastaGz         = gfa2fasta.outputPaternalFastaGz
        File rawhap2FastaGz         = gfa2fasta.outputMaternalFastaGz
        File hap1FastaGz            = asm_group_xy.outputHap1FastaGz
        File hap2FastaGz            = asm_group_xy.outputHap2FastaGz
        File hap1ContigGfaTarGz     = hicHifiasm.outputHaplotype1ContigGfa 
        File hap2ContigGfaTarGz     = hicHifiasm.outputHaplotype2ContigGfa 
        File rawUnitigGfaTarGz      = hicHifiasm.outputRawUnitigGfa
        File binFilesTarGz          = hicHifiasm.outputBinFiles
    }

    parameter_meta {
        childID: "Sample ID of the child sample whose reads are going to be assembled"
        isMaleSample: "(Boolean) used in last step to partition sex chromosomes in male samples"
        childReadsHiFi: "An array of files (or a single file) that contain the HiFi reads of the child sample ( Acceptable formats are fastq (or fq), fastq.gz (or fq.gz), bam and cram)"
        referenceFasta: "(optional) If any of the read files are in cram format, the reference genome should be provided in .fasta format"
        threadCount: "The number of cores for running hifiasm"
        memSize: "The memory size (GB) for running hifiasm"
        preemptible: "The number of tries for using a preemptible node for running hifiasm. Note that if your child data has a coverage of more than 40X, hifiasm (without any given bin files) may take longer than 24 hours. So using a preemptible node is useless beacuse it gets interrupted after 24 hours"
    }
    meta {
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
    }
}

