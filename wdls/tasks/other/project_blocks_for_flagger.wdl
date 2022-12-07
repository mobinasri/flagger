version 1.0 

import "../alignment/long_read_aligner.wdl" as aligner_t
import "../alignment/bam2paf.wdl" as bam2paf_t
import "project_blocks.wdl" as project_blocks_t

workflow runProjectBlocksForFlagger{
    input{
        File hap1AssemblyFastaGz
        File hap2AssemblyFastaGz
        File refAssemblyFastaGz
        Array[File] refBiasedBlocksBedArray
        Array[File] biasedBlocksNameStringArray
        File refSexBed
        File refSDBed
        File refCntrBed
        # splitAssembly is recommended to be true if the assembly is having almost 
        # the same quality as reference (including complete centromeres) since minimap2/winnowmap tests showed
        # that it takes forever to do the alignment of such assemblies to a reference like chm13v2.0
        Boolean splitAssembly = false # split contigs before alignment
        Int splitSize = 20000000 # maximum size (in bases) of split contigs
        String sampleName
        Int mergingMargin = 100000 # merge projected blocks closer than 100kb
        String aligner = "winnowmap"
        String preset = "asm5"
        Int kmerSize = 19
        String zones = "use-west2-a"
    }

    if (splitAssembly){
        call runSplitAssembly as splitHap1{
            input:
                assemblyFastaGz = hap1AssemblyFastaGz,
                splitSize = splitSize
        }
        call runSplitAssembly as splitHap2{
            input:
                assemblyFastaGz = hap2AssemblyFastaGz,
                splitSize = splitSize
        } 
    }

    File hap1AssemblyFastaGzProcessed = select_first([splitHap1.splitAssemblyFastaGz, hap1AssemblyFastaGz])
    File hap2AssemblyFastaGzProcessed = select_first([splitHap2.splitAssemblyFastaGz, hap2AssemblyFastaGz])

    String refSuffix = sub(basename("${refAssemblyFastaGz}"), ".f(ast)?a.gz", "")
    # Align hap1 assembly to reference
    call aligner_t.alignmentBam as alignmentHap1{
        input:
            aligner =  aligner,
            preset = preset,
            suffix = refSuffix,
            refAssembly = refAssemblyFastaGz,
            readFastq_or_queryAssembly = hap1AssemblyFastaGzProcessed,
            kmerSize = kmerSize,
            zones = zones
    }
    # Align hap2 assembly to reference
    call aligner_t.alignmentBam as alignmentHap2{
        input:
            aligner =  aligner,
            preset = preset,
            suffix = refSuffix,
            refAssembly = refAssemblyFastaGz,
            readFastq_or_queryAssembly = hap2AssemblyFastaGzProcessed,
            kmerSize = kmerSize,
            zones = zones
    }

    # Convert hap1 bam file to paf file
    call bam2paf_t.bam2paf as bam2pafHap1{
       input:
           bamFile = alignmentHap1.sortedBamFile,
           minMAPQ = 0,
           primaryOnly = "yes"
    }

    # Convert hap2 bam file to paf file
    call bam2paf_t.bam2paf as bam2pafHap2{
       input:
           bamFile = alignmentHap2.sortedBamFile,
           minMAPQ = 0,
           primaryOnly = "yes"
    }

    # Project ref biased blocks to hap1 and hap2 assemblies separately
    scatter (blocksBed_suffix in zip(refBiasedBlocksBedArray, biasedBlocksNameStringArray)){
        call project_blocks_t.project as projectHap1{
            input:
                blocksBed = blocksBed_suffix.left,
                asm2refPaf = bam2pafHap1.pafFile,
                sampleName = sampleName,
                suffix = "${blocksBed_suffix.right}.hap1",
                mode = "ref2asm",
                mergingMargin = mergingMargin,
                isAssemblySplit = splitAssembly
        }
        call project_blocks_t.project as projectHap2{
            input:
                blocksBed = blocksBed_suffix.left,
                asm2refPaf = bam2pafHap2.pafFile,
                sampleName = sampleName,
                suffix = "${blocksBed_suffix.right}.hap2",
                mode = "ref2asm",
                mergingMargin = mergingMargin,
                isAssemblySplit = splitAssembly
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
            isAssemblySplit = splitAssembly
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
            isAssemblySplit = splitAssembly
    }
    
    # Project Cntr blocks
    call project_blocks_t.project as projectCntr{
        input:
            blocksBed = refCntrBed,
            asm2refPaf = concatPaf.outputFile,
            sampleName = sampleName,
            suffix = "Cntr_Projected",
            mode = "ref2asm",
            mergingMargin = 1,
            isAssemblySplit = splitAssembly
    }
    

    output {
        Array[File] projectionBiasedBedArray = flatten([projectHap1.projectionBed, projectHap2.projectionBed])
        File projectionSDBed = projectSD.projectionBed
        File projectionSexBed = projectSex.projectionBed
        File projectionCntrBed = projectCntr.projectionBed
        File hap1ToRefAlignmentBam = alignmentHap1.sortedBamFile
        File hap2ToRefAlignmentBam = alignmentHap2.sortedBamFile
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


task runSplitAssembly {
    input {
        File assemblyFastaGz
        Int splitSize
        # runtime configurations
        Int memSize=8
        Int threadCount=8
        Int diskSize=128
        String dockerImage="mobinasri/flagger:v0.2"
        Int preemptible=2
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        FILENAME=$(basename ~{assemblyFastaGz})
        PREFIX=${FILENAME%%.f(ast)?a.gz} 

        gunzip -c ~{assemblyFastaGz} > ${PREFIX}.fa
        samtools faidx ${PREFIX}.fa
        python3 /home/programs/src/split_fai_by_length.py --fai asm.fa.fai --splitSize ~{splitSize} > ${PREFIX}.bed
        bedtools getfasta -fi ${PREFIX}.fa -bed ${PREFIX}.bed -fo ${PREFIX}.split.fa
        pigz -p8 ${PREFIX}.split.fa
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File splitAssemblyFastaGz = glob("*.split.fa.gz")
    }
}

