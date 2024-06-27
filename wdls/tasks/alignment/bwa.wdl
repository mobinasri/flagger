version 1.0


import "../../../ext/hpp_production_workflows/QC/wdl/tasks/extract_reads.wdl" as extractReads_t


workflow bwaAlignment{
    input {
        String sampleName
        String suffix
        File cramFile
        File assembly
        File referenceFasta
    }

    ## build bwa index files for the assembly
    call buildBwaIndex{
        input:
            assembly = assembly
    }
    call extractReads_t.extractReads as extractReads {
        input:
            readFile=cramFile,
            referenceFasta=referenceFasta,
            memSizeGB=4,
            threadCount=4,
            diskSizeGB=ceil(size(cramFile, "GB") * 3) + 64
    }
    call BwaAlignment{
        input:
            referenceFasta = referenceFasta,
            readFastq = extractReads.extractedRead,
            indexTar = buildBwaIndex.indexTar,
            outputName = "${sampleName}_${suffix}",
            diskSize = ceil(size(extractReads.extractedRead, "GB") * 2) + 64
    }
    output{
        File sortedBamFile = BwaAlignment.sortedBamFile
        File sortedBaiFile = BwaAlignment.sortedBaiFile
    }
 
}

task buildBwaIndex{
    input{
        File assembly
        # runtime configurations
        Int memSize=16
        Int threadCount=4
        Int diskSize=64
        String dockerImage="quay.io/masri2019/hpp_bwa:latest"
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
        
        if [[ ~{assembly} =~ .*f(ast)?a\.gz$ ]] ; then    
            zcat ~{assembly} > asm.fa
        elif [[ ~{assembly} =~ .*f(ast)?a$ ]] ; then
            cp ~{assembly} asm.fa
        else
             echo "UNSUPPORTED READ FORMAT (expect .fa .fasta .fa.gz .fasta.gz): $(basename ~{assembly})"
             exit 1
        fi

        # build bwa index for the given assembly
        bwa index asm.fa
        mkdir index
        mv asm.* index/
        tar -cf index.tar index
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }

    output {
        File indexTar = "index.tar"
    }
}

task BwaAlignment{
    input{
        File readFastq
        File indexTar
        String outputName
        File? referenceFasta
        String bwaParams = ""
        # runtime configurations
        Int memSize=64
        Int threadCount=32
        Int diskSize=512
        String dockerImage="quay.io/masri2019/hpp_bwa:latest"
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

        # extract the previously generated bwa index
        tar -xf ~{indexTar} --strip-components 1

        # bwa alignment
        bwa mem ~{bwaParams} -t~{threadCount} asm.fa ~{readFastq} | samtools view -b -h > ~{outputName}.bam
        samtools sort -@~{threadCount} -o ~{outputName}.sorted.bam ~{outputName}.bam
        samtools index ~{outputName}.sorted.bam
    >>>

    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File sortedBamFile = glob("*.sorted.bam")[0]
        File sortedBaiFile = glob("*.sorted.bam.bai")[0]
    }
}


