version 1.0 


task getCanonicalBasesBed {
    input {
        File assemblyFastaGz
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=32
        String dockerImage="mobinasri/flagger:v0.4.0"
        Int preemptible=2
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        FILENAME=$(basename ~{assemblyFastaGz})
        PREFIX=${FILENAME%%.fa*(sta).gz}
        ln -s ~{assemblyFastaGz} ${PREFIX}.fa.gz
        gunzip -c ${PREFIX}.fa.gz > ${PREFIX}.fa
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
        if [[ ${EXTENSION} == "gz" ]];then
            ln -s ~{inputFile} output/${FILENAME}
        else
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
        File gzCompressedFile = glob("output/*")[0]
    }
}


task createDipAsm {
    input {
        File hap1AssemblyFastaGz
        File hap2AssemblyFastaGz
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

        zcat ~{hap1AssemblyFastaGz} ~{hap2AssemblyFastaGz} > ~{outputName}.fa
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
