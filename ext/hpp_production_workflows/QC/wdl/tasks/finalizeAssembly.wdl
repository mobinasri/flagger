version 1.0

workflow finalizeAssembly {

    call renameContigsAddMT 

    output {
        File FinalAssembly  = renameContigsAddMT.FinalAssembly
    }
}

task renameContigsAddMT {

    input {
        String sampleName
        String outputFileTag
        Int mat_pat_int
        File inputFastaGZ
        
        File? mitoAssembly

        Int memSizeGB = 4
        Int diskSizeGB = 64
        String dockerImage = "biocontainers/samtools:v1.9-4-deb_cv1"
    }

    String unzippedOrigFa = basename(inputFastaGZ, ".gz")
    String renamedFasta   = "${sampleName}.renamed.fa"
    String outputFasta    = "${sampleName}.${outputFileTag}.fa.gz"

    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace

        ## gunzip input fasta 
        gunzip -c ~{inputFastaGZ} > ~{unzippedOrigFa}

        ## Rename contig names to sampleName#1/2#contigName format (1 = paternal, 2 = maternal)
        sed "s/^>/>~{sampleName}\#~{mat_pat_int}\#/" ~{unzippedOrigFa} > ~{renamedFasta}

        ## Check if we have a mito assembly (not all samples do). If so, add it to maternal haplotype
        if [ -s ~{mitoAssembly} ]; then

            ## Only add to maternal haplotype
            if [[ ~{mat_pat_int} == 2 ]]; then
                cat ~{renamedFasta} ~{mitoAssembly} | gzip > ~{outputFasta}
            else
                cat ~{renamedFasta} | gzip > ~{outputFasta}
            fi
        else
            ## There is no mito assembly, just copy over the file
            cat ~{renamedFasta} | gzip > ~{outputFasta}
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