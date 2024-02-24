version 1.0

import "../../../ext/hpp_production_workflows/QC/wdl/tasks/extract_reads.wdl" as extractReads_t
import "../../../ext/hpp_production_workflows/QC/wdl/tasks/arithmetic.wdl" as arithmetic_t
import "merge_bams.wdl" as mergeBams_t

workflow longReadAlignment {
    input {
        String aligner="winnowmap"
        String preset
        String sampleName
        String sampleSuffix
        Array[File] readFiles
        File assembly
        File? referenceFasta
        Int preemptible=2
        Int extractReadsDiskSize=256
        String zones
    }

    scatter (readFile in readFiles) {
        call extractReads_t.extractReads as extractReads {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=extractReadsDiskSize,
                dockerImage="tpesout/hpp_base:latest"
        }

         ## align reads to the assembly
         call alignmentBam as alignment{
             input:
                 aligner =  aligner,
                 preset = preset,
                 refAssembly=assembly,
                 readFastq_or_queryAssembly = extractReads.extractedRead,
                 diskSize = extractReads.fileSizeGB * 3,
                 preemptible = preemptible,
                 zones = zones
        }
    }

    call arithmetic_t.sum as bamSize {
        input:
            integers=alignment.fileSizeGB
    }

    ## merge the bam files
    call mergeBams_t.merge as mergeBams{
        input:
            sampleName = "${sampleName}.${sampleSuffix}",
            sortedBamFiles = alignment.sortedBamFile,
            # runtime configurations
            diskSize = floor(bamSize.value * 2.5),
            preemptible = preemptible,
            zones = zones
    }
    output {
        File sortedBamFile = mergeBams.mergedBam
        File baiFile = mergeBams.mergedBai
    }

}

task indexBam{
    input{
        String bam 
        # runtime configurations
        Int memSize=16
        Int threadCount=8
        Int diskSize=ceil(size(bam, "GB")) + 32
        String dockerImage="mobinasri/long_read_aligner:v0.3.0"
        Int preemptible=2
        String zones="us-west2-a"
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        BAM_FILENAME=$(basename ~{bam})
        BAM_PREFIX=${BAM_FILENAME%%.bam}

        ln -s ~{bam} ${BAM_PREFIX}.bam
        samtools index ${BAM_PREFIX}.bam
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones: zones
    }
    output {
        File bamIndex = glob("*.bai")[0]
    }
}

task alignmentBam{
    input{
        String aligner
        # For winnowmap-v2.03 -> map-pb/map-ont/asm5/asm10
        # For minimap2-v2.24 -> map-pb/map-hifi/map-ont/asm5/asm10 
        # For veritymap v2.1.2-alpha -> hifi-haploid/hifi-haploid-complete/hifi-diploid/ont-haploid-complete
        String preset
        String suffix=""
        String options=""
        File readFastq_or_queryAssembly
        File refAssembly
        Int kmerSize=15
        # runtime configurations    
        Int memSize=64
        Int threadCount=32
        Int diskSize
        String dockerImage="mobinasri/long_read_aligner:v0.3.0"
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

       

        gunzip -c ~{refAssembly} > asm.fa
        # Sort fasta based on contig names
        seqkit sort -nN asm.fa > asm.sorted.fa	

        REF_FILENAME=$(basename ~{refAssembly})
        REF_FILENAME_NO_GZ=${QUERY_FILENAME%.gz}
        REF_PREFIX=${REF_FILENAME_NO_GZ%.*}

        QUERY_FILENAME=$(basename ~{readFastq_or_queryAssembly})
        QUERY_FILENAME_NO_GZ=${QUERY_FILENAME%.gz}
        QUERY_PREFIX=${QUERY_FILENAME_NO_GZ%.*}
        
        ALIGNMENT_PREFIX="${QUERY_PREFIX}.to_${REF_PREFIX}"
        
        if [[ ~{aligner} == "winnowmap" ]]; then
            # Run meryl for winnowmap
            meryl count k=~{kmerSize} output merylDB asm.sorted.fa
            meryl print greater-than distinct=0.9998 merylDB > repetitive_k~{kmerSize}.txt

            # run winnowmap
            winnowmap -W repetitive_k~{kmerSize}.txt -a -x ~{preset} ~{options} -t~{threadCount} asm.sorted.fa ~{readFastq_or_queryAssembly} | samtools view -h -b > ${ALIGNMENT_PREFIX}.bam
        elif [[ ~{aligner} == "minimap2" ]] ; then
            # Run minimap2
            minimap2 -k ~{kmerSize} -a -x ~{preset} ~{options} -t~{threadCount} asm.sorted.fa ~{readFastq_or_queryAssembly} | samtools view -h -b > ${ALIGNMENT_PREFIX}.bam
        elif [[ ~{aligner} == "veritymap" ]] ; then
            # Run veritymap
            python3 ${VERITY_MAP_PY} --reads ~{readFastq_or_queryAssembly} -o output -t~{threadCount} -d ~{preset} ~{options} asm.sorted.fa

            # Convert sam to bam, sort and index bam file
            SAM_PATH=$(ls output/*.sam)
            SAM_FILENAME=$(basename ${SAM_PATH})
            SAM_PREFIX=${SAM_FILENAME%%.sam}
            samtools view -hb output/${SAM_PREFIX}.sam > ${ALIGNMENT_PREFIX}.bam
            # To save some space
            rm -rf output/${SAM_PREFIX}.sam

        else
             echo "UNSUPPORTED ALIGNER (expect minimap2 or winnowmap): ~{aligner}"
             exit 1
        fi
        
        if [ -z ~{suffix} ]; then
            OUTPUT_FILE=${ALIGNMENT_PREFIX}.sorted.bam
        else
            OUTPUT_FILE=${ALIGNMENT_PREFIX}.~{suffix}.sorted.bam  
        fi
        samtools sort -@~{threadCount} -o ${OUTPUT_FILE} ${ALIGNMENT_PREFIX}.bam
        du -s -BG ${OUTPUT_FILE} | sed 's/G.*//' > outputsize.txt
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones: zones
        cpuPlatform: "Intel Cascade Lake"
    }
    output {
        File sortedBamFile = glob("*.sorted.bam")[0]
        Int fileSizeGB = read_int("outputsize.txt")
    }
}

task alignmentPaf{
    input{
        String aligner
        String preset
        String suffix=""
        String options=""
        File readFastq_or_queryAssembly
        File refAssembly
        Int kmerSize=15
        # runtime configurations    
        Int memSize=64
        Int threadCount=32
        Int diskSize
        String dockerImage="mobinasri/long_read_aligner:v0.1"
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

        
        if [[ ~{aligner} == "winnowmap" ]]; then
            meryl count k=~{kmerSize} output merylDB ~{refAssembly}
            meryl print greater-than distinct=0.9998 merylDB > repetitive_k~{kmerSize}.txt
            ALIGNER_CMD="winnowmap -W repetitive_k~{kmerSize}.txt"
        elif [[ ~{aligner} == "minimap2" ]] ; then
            ALIGNER_CMD="minimap2"
        else
             echo "UNSUPPORTED ALIGNER (expect minimap2 or winnowmap): ~{aligner}"
             exit 1
        fi
        
        fileBasename=$(basename ~{readFastq_or_queryAssembly})

        if [ -z ~{suffix} ]; then
            OUTPUT_FILE=${fileBasename%.*.*}.paf
        else
            OUTPUT_FILE=${fileBasename%.*.*}.~{suffix}.paf
        fi

        if [ -z ~{preset} ]; then
            ${ALIGNER_CMD} ~{options} -t~{threadCount} ~{refAssembly} ~{readFastq_or_queryAssembly} > ${OUTPUT_FILE}
        else
            ${ALIGNER_CMD} ~{options} -x ~{preset} -t~{threadCount} ~{refAssembly} ~{readFastq_or_queryAssembly} > ${OUTPUT_FILE}
        fi

        du -s -BG ${OUTPUT_FILE} | sed 's/G.*//' > outputsize.txt
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones: zones
    }
    output {
        File pafFile = glob("*.paf")[0]
        Int fileSizeGB = read_int("outputsize.txt")
    }
}

