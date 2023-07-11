version 1.0

workflow runSecPhase {
    input {
        File inputBam
        File diploidAssemblyFastaGz
        File? phasedVcf
        File? variantBed
        String secphaseOptions = "--hifi"
        String secphaseDockerImage = "mobinasri/secphase:v0.4.2"
        String version = "v0.4.2"
    }
    call sortByNameAndIndex {
        input:
            bamFile = inputBam,
            dockerImage = secphaseDockerImage,
            diskSize = 7 * ceil(size(inputBam, "GB")) + 64
    }

    call secphase {
        input:
            bam = sortByNameAndIndex.outputBam,
            bamSecphaseIndex = sortByNameAndIndex.outputBamSecphaseIndex,
            diploidAssemblyFastaGz = diploidAssemblyFastaGz,
            phasedVcf = phasedVcf,
            variantBed = variantBed,
            options = secphaseOptions,
            dockerImage = secphaseDockerImage,
            prefix = "secphase_${version}",
            diskSize = ceil(size(sortByNameAndIndex.outputBam, "GB")) + 64
    }

    output {
        File outLog = secphase.outLog
        File variantBlocksBed = secphase.variantBlocksBed
        File markerBlocksBed = secphase.markerBlocksBed
        File modifiedReadBlocksVariantsBed = secphase.modifiedReadBlocksVariantsBed
        File modifiedReadBlocksMarkersBed = secphase.modifiedReadBlocksMarkersBed
        File initalVariantBlocksBed = secphase.initalVariantBlocksBed

    }
}

task sortByNameAndIndex {
    input {
        File bamFile
        String excludeSingleAlignment="no"
        # runtime configurations
        Int memSize=16
        Int threadCount=8
        Int diskSize=1024
        String dockerImage="mobinasri/secphase:v0.4.2"
        Int preemptible=2
        String zones="us-west2-a"
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

        BAM_FILENAME=$(basename ~{bamFile})
        BAM_PREFIX=${BAM_FILENAME%.bam}

        mkdir output
        if [ ~{excludeSingleAlignment} == "yes" ]; then
            samtools view ~{bamFile} | cut -f1 | sort | uniq -c > readnames.txt
            cat readnames.txt | awk '$1 > 1' | tr -s ' ' | cut -d ' ' -f3 > selected_readnames.txt
            extract_reads -i ~{bamFile} -o output/${BAM_PREFIX}.bam -r selected_readnames.txt
        else
            ln ~{bamFile} output/${BAM_PREFIX}.bam
        fi
        samtools sort -n -@8 output/${BAM_PREFIX}.bam > output/${BAM_PREFIX}.sorted_by_qname.bam
        secphase_index -i output/${BAM_PREFIX}.sorted_by_qname.bam
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones : zones
    }
    output {
        File outputBam = glob("output/*.sorted_by_qname.bam")[0]
        File outputBamSecphaseIndex = glob("output/*.bam.secphase.index")[0]
    }
}

task secphaseIndex {
    input {
        File bam
        # runtime configurations
        Int memSize=8
        Int threadCount=4
        Int diskSize=128
        String dockerImage="mobinasri/secphase:v0.4.2"
        Int preemptible=2
        String zones="us-west2-a"
    }
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        ##set -o pipefail
        # cause a bash script to exit immediately when a command fails
        ##set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        ##set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        ##set -o xtrace

        BAM_FILENAME=$(basename ~{bam})
        BAM_PREFIX=${BAM_FILENAME%.bam}

        ln ~{bam} ${BAM_PREFIX}.bam

        secphase_index -i ${BAM_PREFIX}.bam
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones : zones
    }
    output {
        File bamSecphaseIndex = glob("*.bam.secphase.index")[0]
    }
}


task secphase {
    input {
        File bam
        File bamSecphaseIndex
        File diploidAssemblyFastaGz
        File? phasedVcf
        File? variantBed
        String options = "--hifi"
        String prefix = "secphase"
        # runtime configurations
        Int memSize=32
        Int threadCount=32
        Int diskSize=128
        String dockerImage="mobinasri/secphase:v0.4.2"
        Int preemptible=2
        String zones="us-west2-a"
    }
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        ##set -o pipefail
        # cause a bash script to exit immediately when a command fails
        ##set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        ##set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        ##set -o xtrace

        BAM_FILENAME=$(basename ~{bam})
        BAM_PREFIX=${BAM_FILENAME%.bam}

        ln ~{bam} ${BAM_PREFIX}.bam
        ln ~{bamSecphaseIndex} ${BAM_PREFIX}.bam.secphase.index

        ln ~{diploidAssemblyFastaGz} asm.fa.gz
        gunzip -c asm.fa.gz > asm.fa
        samtools faidx asm.fa

        mkdir output
        if [[ -n "~{phasedVcf}" ]];then
            ln ~{phasedVcf} phased.vcf
            ln ~{variantBed} variant.bed
            echo "Running variant/marker dual mode"
            secphase ~{options} -@~{threadCount}  -i ${BAM_PREFIX}.bam -f asm.fa --outDir output --prefix ~{prefix} -v phased.vcf -B variant.bed
        else
            echo "Running marker mode"
            secphase ~{options} -@~{threadCount}  -i ${BAM_PREFIX}.bam -f asm.fa --outDir output --prefix ~{prefix}
        fi
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones : zones
    }
    output {
        File outLog = "output/~{prefix}.out.log"
        File initalVariantBlocksBed = "output/~{prefix}.initial_variant_blocks.bed"
        File modifiedReadBlocksVariantsBed = "output/~{prefix}.modified_read_blocks.variants.bed"
        File modifiedReadBlocksMarkersBed = "output/~{prefix}.modified_read_blocks.markers.bed"
        File markerBlocksBed = "output/~{prefix}.marker_blocks.bed"
        File variantBlocksBed = "output/~{prefix}.variant_blocks.bed"
    }
}

task concatLogs {
    input {
        Array[File] logs
        String filename
        # runtime configurations
        Int memSize=2
        Int threadCount=1
        Int diskSize=32
        String dockerImage="mobinasri/secphase:v0.4.2"
        Int preemptible=2
        String zones="us-west2-a"
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
        cat ~{sep=" " logs} > output/~{filename}.txt
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones : zones
    }
    output {
        File log = glob("output/*.txt")[0]
    }
}

task mergeBeds {
    input {
        Array[File] beds
        String filename
        # runtime configurations
        Int memSize=2
        Int threadCount=1
        Int diskSize=32
        String dockerImage="mobinasri/secphase:v0.4.2"
        Int preemptible=2
        String zones="us-west2-a"
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
        cat ~{sep=" " beds} | bedtools sort -i - | bedtools merge -i - > output/~{filename}
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones : zones
    }
    output {
        File mergedBed = glob("output/*.bed")[0]
    }
}

task splitByName {
    input {
        File bamFile
        Int NReadsPerBam = 400000
        # runtime configurations
        Int memSize=16
        Int threadCount=8
        Int diskSize=512
        String dockerImage="mobinasri/secphase:v0.4.2"
        Int preemptible=2
        String zones="us-west2-a"
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

        BAM_FILENAME=$(basename ~{bamFile})
        BAM_PREFIX=${BAM_FILENAME%.bam}

        mkdir output
        split_bam_by_readname -i ~{bamFile} -o output -n ~{NReadsPerBam} -t~{threadCount}
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones : zones
    }
    output {
        Array[File] splitBams = glob("output/*.bam")
    }
}

task sortByName {
    input {
        File bamFile
        String excludeSingleAlignment="yes"
        # runtime configurations
        Int memSize=16
        Int threadCount=8
        Int diskSize=1024
        String dockerImage="mobinasri/secphase:v0.4.2"
        Int preemptible=2
        String zones="us-west2-a"
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

        BAM_FILENAME=$(basename ~{bamFile})
        BAM_PREFIX=${BAM_FILENAME%.bam}

        mkdir output
        if [ ~{excludeSingleAlignment} == "yes" ]; then
        samtools view ~{bamFile} | cut -f1 | sort | uniq -c > readnames.txt
        cat readnames.txt | awk '$1 > 1' | tr -s ' ' | cut -d ' ' -f3 > selected_readnames.txt
        extract_reads -i ~{bamFile} -o output/${BAM_PREFIX}.bam -r selected_readnames.txt
        else
        ln ~{bamFile} output/${BAM_PREFIX}.bam
        fi
        samtools sort -n -@8 output/${BAM_PREFIX}.bam > output/${BAM_PREFIX}.sorted_by_qname.bam
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones : zones
    }
    output {
        File outputBam = glob("output/*.sorted_by_qname.bam")[0]
    }
}

task sortByContig {
    input {
        File bamFile
        # runtime configurations
        Int memSize=8
        Int threadCount=4
        Int diskSize=128
        String dockerImage="mobinasri/secphase:v0.4.2"
        Int preemptible=2
        String zones="us-west2-a"
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

        BAM_FILENAME=$(basename ~{bamFile})
        BAM_PREFIX=${BAM_FILENAME%.bam}

        mkdir output
        samtools sort -@~{threadCount} ~{bamFile} > output/${BAM_PREFIX}.bam
        samtools index output/${BAM_PREFIX}.bam
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones : zones
    }
    output {
        File outputBam = glob("output/*.bam")[0]
        File outputBai = glob("output/*.bai")[0]
    }
}

