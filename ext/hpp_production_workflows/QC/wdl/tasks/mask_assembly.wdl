version 1.0

workflow runMaskAssembly {

    call maskAssembly 

    output {
        File FinalAssembly = maskAssembly.FinalAssembly
    }
}



task maskAssembly {

    input {
        String sampleName
        String outputFileTag
        File inputFasta
        
        File? adapterBed 

        Int memSize = 4
        Int diskSizeGB = 64
        String dockerImage = "biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1"
    }

    String outputFasta   = "${sampleName}.${outputFileTag}.fa"
    String outputFastaGz = "${sampleName}.${outputFileTag}.fa.gz"

    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace

        
        inputFastaFN=$(basename -- "~{inputFasta}")
        
        ## if adapterBed exists and is not empty: mask regions in inputFasta
        if [ -s ~{adapterBed} ]; then

            ## first check if inputFasta needs to be unzipped
            if [[ $inputFastaFN =~ \.gz$ ]]; then
                cp ~{inputFasta} .
                gunzip -f $inputFastaFN
                inputFastaFN="${inputFastaFN%.gz}"
            else
                ln -s ~{inputFasta}
            fi 


            ## mask fasta in adapterBed regions
            bedtools maskfasta -fi ${inputFastaFN} -bed ~{adapterBed} -fo ~{outputFasta}

            gzip ~{outputFasta}


        ## adapterBed is empty or doesn't exist: copy existing inputFasta to outputFasta (ensure is gzipped)
        else
            ## If already gzipped, just copy over
            if [[ $inputFastaFN =~ \.gz$ ]]; then
                cp ~{inputFasta} ~{outputFastaGz}

            ## Otherwise copy over and gzip
            else
                cat ~{inputFasta} | gzip > ~{outputFastaGz}
            fi
        fi 


    >>>

    output {

        File FinalAssembly  = outputFastaGz
    }

    runtime {
        memory: memSize
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}