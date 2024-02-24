version 1.0 

import "../alignment/bam2paf.wdl" as bam2paf_t

workflow runProjectBlocks {
    input {
        String assemblyFastaGz
        File asm2refBam
        File blocksBed
        String mode='ref2asm'
        String suffix
    }
    call bam2paf_t.bam2paf {
       input:
           bamFile = asm2refBam,
           minMAPQ = 0,
           primaryOnly = "yes"
    }
    call project {
        input:
            blocksBed = blocksBed,
            asm2refPaf = bam2paf.pafFile,
            sampleName = basename("${assemblyFastaGz}", ".fa.gz"),
            suffix = suffix,
            mode = mode
    }
    output {
       File projectionBed = project.projectionBed
    }
}

task project {
    input {
        File blocksBed
        File asm2refPaf
        String sampleName
        String suffix
        String mode
        Int mergingMargin = 1
        String projectOptions = "--divergence"
        Boolean mergeOutput=true
        # The parameter, isAssemblySplit, is added to handle rare cases where 
        # the assembly contigs are split. For example for HG002_T2T_v0.6 assembly
        # the minimap2/winnowmap took forever for aligning to chm13 so splitting 
        # contigs shorten the runtime dramatically. Each split contig name
        # should be like "${original_contig_name}:${start}-${end}"
        Boolean isAssemblySplit=false
        # runtime configurations
        Int memSize=4
        Int threadCount=8
        Int diskSize=32
        String dockerImage="mobinasri/flagger:dev-v0.1"
        Int preemptible=2
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace
        
        if [ -z ~{suffix} ]; then
            OUTPUT_FILENAME=~{sampleName}.bed
        else
            OUTPUT_FILENAME=~{sampleName}.~{suffix}.bed
        fi

        python3 ${PROJECT_BLOCKS_MULTI_THREADED_PY} ~{projectOptions} --mode ~{mode} --paf ~{asm2refPaf} --blocks ~{blocksBed} --outputProjectable projectable.bed --outputProjection projection.bed --threads ~{threadCount}
        mkdir output
        if [[ ~{isAssemblySplit} == "false" ]]
        then
            bedtools sort -i projection.bed > output/output.tmp.bed
            if [[ ~{mergeOutput} == "false" ]]
            then
                mv output/output.tmp.bed output/${OUTPUT_FILENAME}
            else
                cat output/output.tmp.bed | cut -f1-3 | bedtools merge -d ~{mergingMargin} -i - > output/${OUTPUT_FILENAME}
            fi
        else
            # Convert bed coordinates to the originial one if assembly was split before alignment
            python3 /home/programs/src/convert_bed_coors.py projection.bed > projection_orig_coors.bed
            bedtools sort -i projection_orig_coors.bed > output/output.tmp.bed
            if [[ ~{mergeOutput} == "false" ]]
            then
                mv output/output.tmp.bed output/${OUTPUT_FILENAME}
            else
                cat output/output.tmp.bed | cut -f1-3 | bedtools merge -d ~{mergingMargin} -i - > output/${OUTPUT_FILENAME}
            fi
        fi

    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File projectionBed = glob("output/*.bed")[0]
    }
}
