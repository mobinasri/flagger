version 1.0

workflow runMinigraphMap {

    call minigraphMap 

    output {
        File outputPafGZ = minigraphMap.outputPafGZ
    }
}



task minigraphMap {

    input {
        File inputFasta
        File inputGraph
        String sampleName
        String outputFileTag = "minigraph"
        String args = ""
        
        Int threadCount = 8
        Int memSizeGB   = 64
        Int diskSizeGB  = 64
        String dockerImage = "humanpangenomics/hpp_minigraph:latest"
    }

    String outputPafGz = "${sampleName}.${outputFileTag}.paf.gz"

    command <<<

        set -eux -o pipefail
        

        inputFastaFN=$(basename -- "~{inputFasta}")

        ## first check if inputFasta needs to be unzipped
        if [[ $inputFastaFN =~ \.gz$ ]]; then
            cp ~{inputFasta} .
            gunzip -f $inputFastaFN
            inputFastaFN="${inputFastaFN%.gz}"
        else
            ln -s ~{inputFasta}
        fi 


        ## Run minigraph
        minigraph -t ~{threadCount} ~{args} ~{inputGraph} ${inputFastaFN} | gzip > ~{outputPafGz}

    >>>

    output {

        File outputPafGZ  = outputPafGz
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}