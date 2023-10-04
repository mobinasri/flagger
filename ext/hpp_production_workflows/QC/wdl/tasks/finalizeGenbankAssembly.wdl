version 1.0

workflow finalizeGenbankAssembly {

    call renameAndUnmask 

    output {
        File FinalAssembly  = renameAndUnmask.FinalAssembly
    }
}

task renameAndUnmask {

    input {
        String sampleName
        String outputFileTag
        Int mat_pat_int
        File inputFasta
        

        Int memSize = 4
        Int diskSizeGB = 64
        String dockerImage = "biocontainers/samtools:v1.9-4-deb_cv1"
    }

    String outputFasta    = "${sampleName}.${outputFileTag}.fa.gz"

    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace

        if [[ ~{mat_pat_int} -eq 2 ]]; then
            ## Label MT contig
            sed '/#2#MT/s/\.1 /\.1#MT /' ~{inputFasta} > mt_renamed.fa
        else
            ## Just copy over file
            cp ~{inputFasta} mt_renamed.fa
        fi

        ## Remove everything after first space
        sed '/^>/s/ .*$//' \
            mt_renamed.fa \
            > mt_renamed_space_removed.fa

        ## Add the sample name and haplotype integer to the headers
        sed "s/^>/>~{sampleName}\#~{mat_pat_int}\#/" \
            mt_renamed_space_removed.fa \
            > header_corrected.fa    

        ## remove soft masking
        awk '{if(/^[^>]/)$0=toupper($0);print $0}' \
            header_corrected.fa \
            | gzip > ~{outputFasta}

    >>>

    output {

        File FinalAssembly  = outputFasta
    }

    runtime {
        memory: memSize
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}