version 1.1 

import "../../../ext/hpp_production_workflows/QC/wdl/tasks/extract_reads.wdl" as extractReads_t
import "../../../ext/hpp_production_workflows/QC/wdl/tasks/arithmetic.wdl" as arithmetic_t
import "merge_bams.wdl" as mergeBams_t

task verityMap {
    input {
       	File assemblyFastaGz
        File readsFastq
        String suffix
        String mode # either hifi or ont
        # runtime configurations
        Int memSize=16
        Int threadCount=8
        Int diskSize=128
        String dockerImage="mobinasri/veritymap:latest"
        Int preemptible=2
        String zones = "us-west2-a"
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

        # extract HOR assembly
        FILENAME=$(basename ~{assemblyFastaGz})
        PREFIX=${FILENAME%%.fa*.gz}
        gunzip -c ~{assemblyFastaGz} > ${PREFIX}.fa
        samtools faidx ${PREFIX}.fa

        # run VerityMap
        python3 ${VERITY_MAP_PY} --reads ~{readsFastq} -o output -t~{threadCount} -d ~{mode} --careful ${PREFIX}.fa

        # convert sam to bam, sort and index bam file
        SAM_PATH=$(ls output/*.sam)
        SAM_FILENAME=$(basename ${SAM_PATH})
        SAM_PREFIX=${SAM_FILENAME%%.sam} 
        samtools view -hb output/${SAM_PREFIX}.sam > output/${SAM_PREFIX}.bam
        samtools sort output/${SAM_PREFIX}.bam > output/${SAM_PREFIX}.sorted.bam
        samtools index output/${SAM_PREFIX}.sorted.bam

        # Copy bam and bai
        mv output/${SAM_PREFIX}.sorted.bam .
        mv output/${SAM_PREFIX}.sorted.bam.bai .
        
        # Rename output folder and separate fasta file
        mv output ${PREFIX}.~{suffix}
        mkdir output
        tar cvzf ${PREFIX}.~{suffix}.tar.gz ${PREFIX}.~{suffix}

        du -s -BG output/${SAM_PREFIX}.sorted.bam | sed 's/G.*//' > outputsize.txt
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
        File outputTarGz = glob("*.tar.gz")[0]
        File sortedBamFile = glob("*.sorted.bam")[0]
        File sortedBai = glob("*.sorted.bai")[0]
        Int fileSizeGB = read_int("outputsize.txt")
    }
}

