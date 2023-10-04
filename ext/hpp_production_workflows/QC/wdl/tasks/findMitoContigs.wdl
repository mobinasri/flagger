version 1.0

workflow findMitoContigs {

    input {
        File chrM_fa
        String sampleName
        String parent
        File inputFastaGZ
    }
    
    call blastFasta {
        input:
            chrM_fa=chrM_fa,
            sampleName=sampleName,
            parent=parent,
            inputFastaGZ=inputFastaGZ,
    }

    call parseBlastOutput {
        input:
            sampleName=sampleName,
            parent=parent,
            blastOutput=blastFasta.blastOutput
    }

    output {
        File blastOutput       = blastFasta.blastOutput
        File parsedBlastOutput = parseBlastOutput.parsedBlastOutput
    }
}


task blastFasta {

    input {
        File chrM_fa
        String sampleName
        String parent
        File inputFastaGZ
        
        Int memSize = 4
        Int diskSizeGB = 64
        String dockerImage = "ncbi/blast:latest"
    }

    String inputFasta = basename(inputFastaGZ, ".gz")
    String blastOutputName = "${sampleName}.${parent}.BlastOutput.txt"

    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace

        ## Create blastn database
        makeblastdb \
            -in ~{chrM_fa} \
            -title chrM_Reference \
            -dbtype nucl \
            -out chrM_ref/chrM_ref


        ## gunzip input fasta 
        gunzip -c ~{inputFastaGZ} > ~{inputFasta}
        

        ## Query fasta against MT database
        blastn \
            -query ~{inputFasta} \
            -db chrM_ref/chrM_ref \
            -num_threads 2 \
            -dust no \
            -soft_masking false \
            -out ~{blastOutputName} \
            -outfmt '6 std qlen slen'
    >>>

    output {
        File blastOutput = blastOutputName
    }

    runtime {
        memory: memSize
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}


task parseBlastOutput {

    input {
        String sampleName
        String parent
        File blastOutput

        Int memSize = 4
        Int diskSizeGB = 64
        String dockerImage = "juklucas/parse_mito_blast:latest"
    }

    String parsedBlastOutputName = "${sampleName}.${parent}.ParsedBlastOutput.txt"

    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace

        ## parse blast output (aggregate blast hits into contig-level summary; 
        ## only outputs MT-only contigs). Slightly modified version of parse_blast.py 
        ## script found in MitoHiFi (https://github.com/marcelauliano/MitoHiFi)
        parse_mito_blast_results.py ~{blastOutput} ~{parsedBlastOutputName}
        
    >>>

    output {
        File parsedBlastOutput = parsedBlastOutputName
    }

    runtime {
        memory: memSize
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}