version 1.0

import "../tasks/alignment/bam2paf.wdl" as bam2paf_t
import "../tasks/other/project_blocks.wdl" as projectBlocks_t

workflow runProjectFlaggerFinalBed {
    input {
        String sampleName
        File asm2refBamMat
        File asm2refBamPat
        File finalBed
        Array[String] comps = ["Cc", "Hc", "Dd", "Ee", "Dh", "Eh", "Hh", "Ec"]
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
    call extractComps {
        input:
            bed = finalBed,
            comps = comps
    }
    scatter (bed in extractComps.beds) {
        call getHapBed {
            input:
                bed=bed,
                patKeyword = "#1",
                matKeyword = "#2"
        }
        call projectBlocks_t.project as projectPat {
            input:
                blocksBed = getHapBed.patBed,
                asm2refPaf = patPaf.pafFile,
                sampleName = basename("${getHapBed.patBed}", ".bed"),
                suffix = "projected",
                mode="asm2ref"
    	}
        call projectBlocks_t.project as projectMat {
            input:
                blocksBed = getHapBed.matBed,
                asm2refPaf = matPaf.pafFile,
                sampleName = basename("${getHapBed.matBed}", ".bed"),
                suffix = "projected",
                mode="asm2ref"
        }
    }
    # Merge the paternal BED files containing the projected blocks for 
    # different components into a single colored BED file
    call mergeBeds as mergeBedsPat{
        input :
            trackName = "${sampleName}.pat",
            outputName = "${sampleName}.pat.projected",
            beds = projectPat.projectionBed,
            comps = comps
    }
    # Merge the maternal BED files containing the projected blocks for 
    # different components into a single colored BED file
    call mergeBeds as mergeBedsMat{
        input :
            trackName = "${sampleName}.mat",
            outputName = "${sampleName}.mat.projected",
            beds = projectMat.projectionBed,
            comps = comps
    }

    output {
       File patProjectedBed = mergeBedsPat.bed
       File matProjectedBed = mergeBedsMat.bed
    }
}


task extractComps{
    input {
        File bed
        Array[String] comps
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=32
        String dockerImage="mobinasri/bio_base:v0.1"
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
       
        mkdir output 
        for comp in ~{sep=" " comps}
        do 
            cat ~{bed} | awk -v c=${comp} '$4 == c' > output/${comp}.bed
        done
    >>>

    runtime {
        docker: dockerImage
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        Array[File] beds = glob("output/*.bed")
    }
}

task getHapBed {
    input {
        File bed
        String matKeyword
        String patKeyword
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=32
        String dockerImage="mobinasri/bio_base:v0.1"
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

        
        cat ~{bed} | grep "~{patKeyword}" > ${PREFIX}.pat.bed || true
        cat ~{bed} | grep "~{matKeyword}" > ${PREFIX}.mat.bed || true
    >>>

    runtime {
        docker: dockerImage
        memory: memSize
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
        String dockerImage="mobinasri/bio_base:v0.1"
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
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        Array[File] beds = glob("output/*")
    }
}


task mergeBeds {
    input {
        Array[File] beds
        Array[String] comps
        String outputName
        String trackName
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=32
        String dockerImage="mobinasri/flagger:v0.3.1"
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

        mkdir input output
        ln ~{sep=" " beds} input
        for comp in  ~{sep=" " comps}
        do
            BED_COMP=input/$(ls input | grep "${comp}.")
            cat ${BED_COMP} | awk -v comp=${comp} '{print $0"\t"comp}' >> output/all.bed
        done

        echo "track name=\"~{trackName}\" visibility=1 itemRgb="On"" > output/~{outputName}.bed
        bedtools sort -i output/all.bed > output/all.sorted.bed
        awk 'FNR==NR{c[$1]=$2;next}{print $0"\t0\t+\t"$2"\t"$3"\t"c[$4]}' \
            /home/scripts/colors.txt \
            output/all.sorted.bed \
            >> output/~{outputName}.bed
    >>>

    runtime {
        docker: dockerImage
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File bed = "output/${outputName}.bed"
    }
}

