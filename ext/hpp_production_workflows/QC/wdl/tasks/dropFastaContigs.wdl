version 1.0

workflow dropFastaContigs {

    call dropContigs 

    output {
        File FinalAssembly = dropContigs.FinalAssembly
    }
}



task dropContigs {

    input {
        String sampleName
        String outputFileTag
        File inputFasta
        
        File? contigsToDrop 

        Int memSizeGB = 4
        Int diskSizeGB = 64
        String dockerImage = "biocontainers/samtools:v1.9-4-deb_cv1"
    }

    String outputFasta = "${sampleName}.${outputFileTag}.fa.gz"

    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace

        
        inputFastaFN=$(basename -- "~{inputFasta}")
        
        ## if contigsToDrop exists and is not empty: drop from inputFasta
        if [ -s ~{contigsToDrop} ]; then

            ## first check if inputFasta needs to be unzipped
            if [[ $inputFastaFN =~ \.gz$ ]]; then
                cp ~{inputFasta} .
                gunzip -f $inputFastaFN
                inputFastaFN="${inputFastaFN%.gz}"
            else
                ln -s ~{inputFasta}
            fi 


            ## Index inputFasta to get file with contigs names + sizes
            samtools faidx $inputFastaFN
            inputFastaFai="${inputFastaFN}.fai"

            ## Print all contigs in inputFasta. Write to file without the contigs that we are dropping
            cat ${inputFastaFai} | cut -f1 | grep -v -f ~{contigsToDrop} > contigsToKeep.txt

            ## Write contigs to keep to outputFasta
            samtools faidx $inputFastaFN `cat contigsToKeep.txt` | gzip > ~{outputFasta}


        ## contigsToDrop is empty or doesn't exist: copy existing inputFasta to outputFasta (ensure is gzipped)
        else
            ## If already gzipped, just copy over
            if [[ $inputFastaFN =~ \.gz$ ]]; then
                cp ~{inputFasta} ~{outputFasta}

            ## Otherwise copy over and gzip
            else
                cat ~{inputFasta} | gzip > ~{outputFasta}
            fi
        fi 


    >>>

    output {

        File FinalAssembly  = outputFasta
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}