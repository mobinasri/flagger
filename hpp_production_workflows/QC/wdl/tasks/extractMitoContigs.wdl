version 1.0

workflow extractMitoContigs {

    call extractContigs 

    output {
        File OutputContigFasta = extractContigs.OutputContigFasta
    }
}



task extractContigs {

    input {
        String sampleName
        String outputFileTag
        File matInputFasta
        File patInputFasta
        File matContigsToExtract
        File patContigsToExtract 

        Int memSizeGB = 4
        Int diskSizeGB = 64
        String dockerImage = "biocontainers/samtools:v1.9-4-deb_cv1"
    }

    String outputFasta = "${sampleName}.${outputFileTag}.fa"

    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace

        
        matInputFastaFN=$(basename -- "~{matInputFasta}")
        patInputFastaFN=$(basename -- "~{patInputFasta}")        

        ## first check if inputFastas needs to be unzipped

        ## MATERNAL
        if [[ $matInputFastaFN =~ \.gz$ ]]; then
            cp ~{matInputFasta} .
            gunzip -f $matInputFastaFN
            matInputFastaFN="${matInputFastaFN%.gz}"
        else
            ln -s ~{matInputFasta}
        fi 

        # PATERNAL
        if [[ $patInputFastaFN =~ \.gz$ ]]; then
            cp ~{patInputFasta} .
            gunzip -f $patInputFastaFN
            patInputFastaFN="${patInputFastaFN%.gz}"
        else
            ln -s ~{patInputFasta}
        fi 



        ## Extract mito fasta sequences and write to file
        samtools faidx $matInputFastaFN `cat ~{matContigsToExtract}` > mat_contigs.fa
        samtools faidx $patInputFastaFN `cat ~{patContigsToExtract}` > pat_contigs.fa
        
        
        ## Combine contigs from maternal and paternal assemblies
        cat mat_contigs.fa pat_contigs.fa > ~{outputFasta}

    >>>

    output {

        File OutputContigFasta  = outputFasta
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}