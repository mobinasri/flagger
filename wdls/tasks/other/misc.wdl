version 1.0 

task getIndexLabeledBed{
    input{
        File bed
        String suffix=""
        File? canonicalBasesDiploidBed
        Boolean addHapCoordinates = false
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=8
        String dockerImage="mobinasri/flagger:v1.0.0"
        Int preemptible=2
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        BED_FILE="~{bed}"
        PREFIX=$(basename ${BED_FILE%%.bed})

        touch additional_hap.bed

        if [[ -n "~{canonicalBasesDiploidBed}" && -n "${true="ADD_HAP" false="" addHapCoordinates}" ]]; then
            bedtools subtract -a ~{canonicalBasesDiploidBed} -b ~{bed} | \
                awk '{print $1"\t"$2"\t"$3"\tHap\t0\t.\t"$2"\t"$3"\t0,138,0"}' | \
                bedtools sort -i - > additional_hap.bed
        fi


        mkdir output
        if [ -n "~{suffix}" ]
        then
            OUTPUT="output/${PREFIX}.~{suffix}.index_labeled.bed"
        else
            OUTPUT="output/${PREFIX}.index_labeled.bed"
        fi

        cat ~{bed} additional_hap.bed | \
            cut -f1-4 | \
            sed -e 's|Col_Err|3|g' -e 's|Col_Del|3|g' | sed -e 's|Err|0|g' -e 's|Hap|2|g' -e 's|Dup|1|g' -e 's|Col|3|g' | \
            bedtools sort -i - > ${OUTPUT}

    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output{
        File labeledBed = glob("output/*.bed")[0]
    }
}


task getCanonicalBasesBed {
    input {
        File assemblyFasta
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=32
        String dockerImage="mobinasri/flagger:v1.0.0"
        Int preemptible=2
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace


        FA_PREFIX=$(echo $(basename ~{assemblyFasta}) | sed -e 's/\.fa$//' -e 's/\.fa.gz$//' -e 's/\.fasta$//' -e 's/\.fasta.gz$//')

        if [[ "~{assemblyFasta}" == *.gz ]]; then
            gunzip -c ~{assemblyFasta} > ${FA_PREFIX}.fa
        else
            cp ~{assemblyFasta} ${FA_PREFIX}.fa
        fi

        # ignore Ns
        python3 /home/scripts/get_contig_coords.py --inputFasta ${FA_PREFIX}.fa | bedtools sort -i - > ${FA_PREFIX}.canonical_only.bed

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
