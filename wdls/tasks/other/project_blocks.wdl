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
        # runtime configurations
        Int memSize=4
        Int threadCount=8
        Int diskSize=32
        String dockerImage="mobinasri/flagger:dev-v0.1"
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
        
        if [ -z ~{suffix} ]; then
            OUTPUT_FILENAME=~{sampleName}.bed
        else
            OUTPUT_FILENAME=~{sampleName}.~{suffix}.bed
        fi

        python3 ${PROJECT_BLOCKS_MULTI_THREADED_PY} --mode ~{mode} --paf ~{asm2refPaf} --blocks ~{blocksBed} --outputProjectable projectable.bed --outputProjection ${OUTPUT_FILENAME} --threads ~{threadCount}
        mkdir output
        bedtools sort -i ${OUTPUT_FILENAME} | bedtools merge -i - > output/${OUTPUT_FILENAME}

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
