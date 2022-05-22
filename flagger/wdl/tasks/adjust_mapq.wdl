version 1.0 

#import "../../../QC/wdl/tasks/extract_reads.wdl" as extractReads_t
#import "../../../QC/wdl/tasks/arithmetic.wdl" as arithmetic_t
import "read_set_splitter.wdl" as readSetSplitter_t
import "long_read_aligner.wdl" as longReadAligner_t

workflow runAdjustMapq{
    input{
        String sampleName
        String suffix
        File bam
        Int maxMapq
        File matAssemblyFastaGz
        File patAssemblyFastaGz
        Int readSetSplitNumber = 8
        String aligner
        String preset
    }
    call getLowMapq{
        input:
            bam = bam,
            maxMapq = maxMapq
    }
    call phaseBam{
        input:
            bam = getLowMapq.lowMapqBam
    }

    call extractReads as extractReadsPat {
        input:
            readFile=phaseBam.patBam,
            memSizeGB=4,
            threadCount=4,
            diskSizeGB=512
    }
    call readSetSplitter_t.readSetSplitter as readSetSplitterPat {
        input:
            readFastqs = [extractReadsPat.extractedRead],
            splitNumber = readSetSplitNumber,
            diskSize = floor(extractReadsPat.fileSizeGB * 2.5)
    }
    scatter (readFastqAndSize in zip(readSetSplitterPat.splitReadFastqs, readSetSplitterPat.readSizes)) {
         ## align reads to the assembly
         call longReadAligner_t.alignmentBam as alignmentPat{
             input:
                 aligner =  aligner,
                 preset = preset,
                 refAssembly=patAssemblyFastaGz,
                 readFastq_or_queryAssembly = readFastqAndSize.left,
                 diskSize = 8 + floor(readFastqAndSize.right) * 6,
                 preemptible = 2,
                 zones = "us-west2-a"
        }
    }
    call extractReads as extractReadsMat {
        input:
            readFile=phaseBam.matBam,
            memSizeGB=4,
            threadCount=4,
            diskSizeGB=512
    }
    call readSetSplitter_t.readSetSplitter as readSetSplitterMat {
        input:
            readFastqs = [extractReadsMat.extractedRead],
            splitNumber = readSetSplitNumber,
            diskSize = floor(extractReadsMat.fileSizeGB * 2.5)
    }
    scatter (readFastqAndSize in zip(readSetSplitterMat.splitReadFastqs, readSetSplitterMat.readSizes)) {
         ## align reads to the assembly
         call longReadAligner_t.alignmentBam as alignmentMat{
             input:
                 aligner =  aligner,
                 preset = preset,
                 refAssembly=matAssemblyFastaGz,
                 readFastq_or_queryAssembly = readFastqAndSize.left,
                 diskSize = 8 + floor(readFastqAndSize.right) * 6,
                 preemptible = 2,
                 zones = "us-west2-a"
        }
    }
    call getMapqTable {
        input:
            bams = flatten([alignmentMat.sortedBamFile, alignmentPat.sortedBamFile]),
            prefix = "${sampleName}.${suffix}"
    }
    output {
        File mapqTable = getMapqTable.mapqTable
    }
}

task getLowMapq {
    input {
        File bam
        Int maxMapq
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=512
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

        FILENAME=$(basename ~{bam})
        PREFIX=${FILENAME%.bam}
        mkdir output
        samtools view -h ~{bam} | awk 'substr($1,1,1) == "@" || $5 < ~{maxMapq}' > output/${PREFIX}.low_mapq.bam
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File lowMapqBam = glob("output/*.bam")[0]
    }
}

task phaseBam {
    input {
        File bam
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=512
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

        FILENAME=$(basename ~{bam})
        PREFIX=${FILENAME%.bam}
        mkdir output
        samtools view -h ~{bam} | awk 'substr($1,1,1) == "@" || $3 ~ /.*#1#.*/' > output/${PREFIX}.pat.bam
        samtools view -h ~{bam} | awk 'substr($1,1,1) == "@" || $3 ~ /.*#2#.*/' > output/${PREFIX}.mat.bam
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File patBam = glob("output/*.pat.bam")[0]
        File matBam = glob("output/*.mat.bam")[0]
    }
}

task getMapqTable {
    input {
        Array[File] bams
        String prefix
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=512
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
        for bam in ~{sep=" " bams}
        do
            samtools view ${bam} | cut -f1,3,4,5 >> output/~{prefix}.mapq_table.unsorted.txt
        done
        sort -k2,2 -k3nr,3 output/~{prefix}.mapq_table.unsorted.txt > output/~{prefix}.mapq_table.txt
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File mapqTable = glob("output/*.mapq_table.txt")[0]
    }
}


task extractReads {
    input {
        File readFile
        File? referenceFasta
        Int memSizeGB = 4
        Int threadCount = 8
        Int diskSizeGB = 128
        String dockerImage = "quay.io/masri2019/hpp_base:latest"
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
        

        FILENAME=$(basename -- "~{readFile}")
        PREFIX="${FILENAME%.*}"

        mkdir output
        # -F 0x000 is because we may have only the secondary alignment of a read
        # avoid writing a read multiple times (that's what the awk statemen does) because we may have both sec/pri alignments of a single read
        samtools fastq -F 0x000 -@~{threadCount} ~{readFile} | awk '(substr($1,1,1) == "@"){a[$1] += 1; read=$1}(a[read] == 1){print $0}'> output/${PREFIX}.fq

        OUTPUTSIZE=`du -s -BG output/ | sed 's/G.*//'`
        if [[ "0" == $OUTPUTSIZE ]] ; then
            OUTPUTSIZE=`du -s -BG ~{readFile} | sed 's/G.*//'`
        fi
        echo $OUTPUTSIZE >outputsize
    >>>

    output {
        File extractedRead = flatten([glob("output/*"), [readFile]])[0]
        Int fileSizeGB = read_int("outputsize")
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }

    parameter_meta {
        readFile: {description: "Reads file in BAM"}
    }
}


