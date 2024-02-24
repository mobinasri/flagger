version 1.0 

import "../alignment/long_read_aligner.wdl" as aligner_t
import "../alignment/bam2paf.wdl" as bam2paf_t
import "project_blocks.wdl" as project_blocks_t
import "bedtools.wdl" as bedtools_t

workflow runProjectBlocksForFlagger{
    input{
        File hap1AssemblyBam
        File hap2AssemblyBam
        Array[File] refBiasedBlocksBedArray
        File refSexBed
        File refSDBed
        File refCntrBed
        File refCntrCtBed
        # isAssemblySplit should be true if assembly is split before alignment to reference
        Boolean isAssemblySplit = false
        String sampleName
        Int mergingMargin = 100000 # merge projected blocks closer than 100kb
        String zones = "use-west2-a"
    }


    # Convert hap1 bam file to paf file
    call bam2paf_t.bam2paf as bam2pafHap1{
       input:
           bamFile = hap1AssemblyBam,
           minMAPQ = 0,
           primaryOnly = "yes"
    }

    # Convert hap2 bam file to paf file
    call bam2paf_t.bam2paf as bam2pafHap2{
       input:
           bamFile = hap2AssemblyBam,
           minMAPQ = 0,
           primaryOnly = "yes"
    }

    # Project ref biased blocks to hap1 and hap2 assemblies separately
    scatter (bed in refBiasedBlocksBedArray) {
        String bed_suffix = basename(bed, ".txt")
        call project_blocks_t.project as projectHap1{
            input:
                blocksBed = bed,
                asm2refPaf = bam2pafHap1.pafFile,
                sampleName = sampleName,
                suffix = "hap1.${bed_suffix}",
                mode = "ref2asm",
                mergingMargin = mergingMargin,
                isAssemblySplit = isAssemblySplit
        }
        call project_blocks_t.project as projectHap2{
            input:
                blocksBed = bed,
                asm2refPaf = bam2pafHap2.pafFile,
                sampleName = sampleName,
                suffix = "hap2.${bed_suffix}",
                mode = "ref2asm",
                mergingMargin = mergingMargin,
                isAssemblySplit = isAssemblySplit
        }
    }

    # concatenate paf files
    call concatFiles as concatPaf{
        input:
            files = [bam2pafHap1.pafFile, bam2pafHap2.pafFile],
            outputName = "${sampleName}.concatenated",
            extension = "paf"
    }

    # Project SD blocks
    call project_blocks_t.project as projectSD{
        input:
            blocksBed = refSDBed,
            asm2refPaf = concatPaf.outputFile,
            sampleName = sampleName,
            suffix = "SD_Projected",
            mode = "ref2asm",
            mergingMargin = 1,
            isAssemblySplit = isAssemblySplit
    }
    
    # Project Sex blocks
    call project_blocks_t.project as projectSex{
        input:
            blocksBed = refSexBed,
            asm2refPaf = concatPaf.outputFile,
            sampleName = sampleName,
            suffix = "Sex_Projected",
            mode = "ref2asm",
            mergingMargin = 1,
            isAssemblySplit = isAssemblySplit
    }
    
    # Project Cntr blocks
    call project_blocks_t.project as projectCntr{
        input:
            blocksBed = refCntrBed,
            asm2refPaf = concatPaf.outputFile,
            sampleName = sampleName,
            suffix = "Cntr_Projected",
            mode = "ref2asm",
            mergingMargin = 100000,
            isAssemblySplit = isAssemblySplit
    }

    # Project Cntr transition (Ct) blocks
    call project_blocks_t.project as projectCntrCt{
        input:
            blocksBed = refCntrCtBed,
            asm2refPaf = concatPaf.outputFile,
            sampleName = sampleName,
            suffix = "Cntr_Trans_Projected",
            mode = "ref2asm",
            mergingMargin = 10000,
            isAssemblySplit = isAssemblySplit
    }
    
    # Subtract centric transition regions from centromeres
    call bedtools_t.subtract as subtractCntr{
        input:
            firstBed = projectCntr.projectionBed,
            secondBed = projectCntrCt.projectionBed,
            outputPrefix = "${sampleName}" + "censat_no_ct"
    }

    output {
        Array[File] projectionBiasedBedArray = flatten([projectHap1.projectionBed, projectHap2.projectionBed])
        File projectionSDBed = projectSD.projectionBed
        File projectionSexBed = projectSex.projectionBed
        File projectionCntrBed = subtractCntr.subtractBed
    }    
}

task concatFiles {
    input {
        Array[File] files
        String outputName
        String extension
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=128
        String dockerImage="mobinasri/bio_base:v0.1"
        Int preemptible=2
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace
        
        cat ~{sep=" " files} > ~{outputName}.~{extension}
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File outputFile = "~{outputName}.~{extension}"
    }
}


