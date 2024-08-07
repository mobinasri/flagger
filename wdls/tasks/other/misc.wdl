version 1.0 


task getCanonicalBasesBed {
    input {
        File assemblyFasta
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=32
        String dockerImage="mobinasri/flagger:v0.4.0--98d66028a969211773077f005511f3d78afdc21c"
        Int preemptible=2
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        FILENAME=$(basename ~{assemblyFasta})

        EXTENSION=${FILENAME##*.}

        if [[ ${EXTENSION} == "gz" ]];then
            PREFIX=${FILENAME%%.fa*(sta).gz}
            gunzip -c ~{assemblyFasta} > ${PREFIX}.fa
        else
            PREFIX=${FILENAME%%.fa*(sta)}
            ln -s ~{assemblyFasta} ${PREFIX}.fa 
        fi
        

        # ignore Ns
        python3 /home/scripts/get_contig_coords.py --inputFasta ${PREFIX}.fa | bedtools sort -i - > ${PREFIX}.canonical_only.bed

    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File canonicalBasesBed = glob("*.bed")[0]
    }
}

task gzipCompress {
    input {
        File inputFile
        # runtime configurations
        Int memSize=8
        Int threadCount=4
        Int diskSize=32
        String dockerImage="mobinasri/bio_base:v0.4.0"
        Int preemptible=2
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace


        FILENAME=$(basename ~{inputFile})
        EXTENSION=${FILENAME##*.}

        mkdir output
        if [[ ${EXTENSION} != "gz" ]];then
            pigz -p~{threadCount} -c ~{inputFile} > output/${FILENAME}.gz
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
        File gzCompressedFile = flatten([glob("output/*"), [inputFile]])[0]
    }
}


task createDipAsm {
    input {
        File hap1AssemblyFasta
        File hap2AssemblyFasta
        String outputName
        # runtime configurations
        Int memSize=8
        Int threadCount=4
        Int diskSize=32
        String dockerImage="mobinasri/bio_base:v0.4.0"
        Int preemptible=2
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace


        FILENAME=$(basename ~{hap1AssemblyFasta})
        EXTENSION=${FILENAME##*.}
        if [[ ${EXTENSION} == "gz" ]];then
            CAT_COMMAND="zcat"
        else
            CAT_COMMAND="cat"
        fi

        ${CAT_COMMAND} ~{hap1AssemblyFasta} ~{hap2AssemblyFasta} > ~{outputName}.fa
        samtools faidx ~{outputName}.fa
        pigz -p~{threadCount} ~{outputName}.fa

    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File diploidAssemblyFastaGz = glob("*.fa.gz")[0]
        File diploidAssemblyFastaIndex = glob("*.fai")[0]
    }
}

task createFile{
     input {
        String content = ""
        String filename = "mock.txt"
        # runtime configurations
        Int memSize=2
        Int threadCount=2
        Int diskSize=2
        String dockerImage="mobinasri/bio_base:v0.4.0"
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        mkdir output
        echo ~{content} > output/~{filename}
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
    }
    output {
        File outFile = glob("output/*")[0] 
    }
} 
