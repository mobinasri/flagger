version 1.0

import "../tasks/bam2paf.wdl" as bam2paf_t
import "../tasks/project_blocks.wdl" as projectBlocks_t

workflow runProjectBlocksComps {
    input {
        String assemblyFastaGz
        File asm2refBamMat
        File asm2refBamPat
        File bedsTarGz
    }
    call bam2paf_t.bam2paf as matPaf {
       input:
           bamFile = asm2refBamMat,
           minMAPQ = 0,
           primaryOnly = "yes"
    }
    call bam2paf_t.bam2paf as patPaf {
       input:
           bamFile = asm2refBamPat,
           minMAPQ = 0,
           primaryOnly = "yes"
    }
    call extractBeds {
       input:
           bedsTarGz=bedsTarGz
    }
    scatter (bed in extractBeds.beds) {
        call getHapBed {
            input:
                bed=bed
        }
        call projectBlocks_t.project as projectPat {
            input:
                blocksBed = getHapBed.patBed,
                asm2refPaf = patPaf.pafFile,
                sampleName = basename("${getHapBed.patBed}", ".bed"),
                suffix = "chm13_projected",
                mode="asm2ref"
    	}
        call projectBlocks_t.project as projectMat {
            input:
                blocksBed = getHapBed.matBed,
                asm2refPaf = matPaf.pafFile,
                sampleName = basename("${getHapBed.matBed}", ".bed"),
                suffix = "chm13_projected",
                mode="asm2ref"
        }
    }
    call tarBeds as tarBedsPat{
        input :
            internalName="projected",
            tarGzName=basename("${bedsTarGz}", "beds.filtered.tar.gz") + ".pat.projected",
            beds=projectPat.projectionBed
    }
    call tarBeds as tarBedsMat{
        input :
            internalName="projected",
            tarGzName=basename("${bedsTarGz}", "beds.filtered.tar.gz") + ".mat.projected",
            beds=projectMat.projectionBed
    }

    output {
       File patProjectedBedsTarGz = tarBedsPat.bedsTarGz
       File matProjectedBedsTarGz = tarBedsMat.bedsTarGz
    }
}


task getHapBed {
    input {
        File bed
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=32
        String dockerImage="quay.io/masri2019/hpp_base:latest"
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

        FILENAME=~{bed}
        PREFIX=$(basename ${FILENAME%.bed})

        cat ~{bed} | grep "#1" > ${PREFIX}.pat.bed
        cat ~{bed} | grep "#2" > ${PREFIX}.mat.bed
    >>>

    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File patBed = glob("*.pat.bed")[0]
        File matBed = glob("*.mat.bed")[0]
    }
}


task extractBeds {
    input {
        File bedsTarGz
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=32
        String dockerImage="quay.io/masri2019/hpp_base:latest"
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


        mkdir output
        tar --strip-components 1 -xvzf ~{bedsTarGz} --directory output

    >>>

    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        Array[File] beds = glob("output/*")
    }
}


task tarBeds {
    input {
        Array[File] beds
        String internalName
        String tarGzName
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=32
        String dockerImage="quay.io/masri2019/hpp_base:latest"
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

        
        mkdir ~{internalName}
        cp ~{sep=" " beds} ~{internalName}/


        tar -cf ~{tarGzName}.tar ~{internalName}
        gzip ~{tarGzName}.tar

    >>>

    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File bedsTarGz = glob("*.tar.gz")[0]
    }
}

