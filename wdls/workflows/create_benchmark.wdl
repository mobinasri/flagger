version 1.1 

workflow runCreateBenchmark{
    input {
        File refFastaGz
        File asmFastaGz
    }
    ## align query assembly to the ref assembly
    call aligner_t.alignmentBam as asm2ref{
        input:
            aligner = aligner,
            preset = preset,
            suffix = suffix,
            refAssembly = refFastaGz,
            readFastq_or_queryAssembly = asmFastaGz,
            kmerSize = 19,
            diskSize = 64,
            preemptible = 2,
            zones = zones
    }
    call bam2paf {
        File 
    } 
    call getDuplicatedAndErroneous {
        input:
            alignmentBam = alignment
    }
    ## align unmapped ref blocks to the assembly
    call aligner_t.alignmentBam as asm2ref{
        input:
            aligner = aligner,
            preset = preset,
            suffix = suffix,
            refAssembly = asmFastaGz,
            readFastq_or_queryAssembly = get,
            kmerSize = 19,
            diskSize = 64,
            preemptible = 2,
            zones = zones
    }
    call bam2paf {
    }
    call getCollapsed {
        input:
            alignmentBam = alignment
    }
    output {
        File collapsedBed = createBenchmark.collapsedBed
        File duplicate
    }
}

task createBenchmark {
    input {
        File refFastaGz
        File asmFastaGz
        # runtime configurations
        Int memSize=32
        Int threadCount=16
        Int diskSize=32
        String dockerImage="mobinasri/long_read_aligner:v0.4.0"
        Int preemptible=2
    }
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace
      
        # Get file prefixes
        ASM_FILENAME=$(basename ~{asmFastaGz})
        ASM_PREFIX=${ASM_FILENAME%%.fa?(sta).gz}
        
        REF_FILENAME=$(basename ~{refFastaGz})
        REF_PREFIX=${REF_FILENAME%%.fa?(sta).gz}

        # Align assembly to reference
        mkdir output
        minimap2 -ax asm5 --cs --eqx -Y -L ~{refFastaGz} ~{asmFastaGz} | samtools view -hb > output/${ASM_PREFIX}.${REF_PREFIX}.bam
        samtools sort output/${ASM_PREFIX}.${REF_PREFIX}.bam > output/${ASM_PREFIX}.${REF_PREFIX}.sorted.bam
        samtools index output/${ASM_PREFIX}.${REF_PREFIX}.sorted.bam

        # Extract duplicated blocks (blocks with multiple alignments)
        bedtools bamtobed -i output/${ASM_PREFIX}.${REF_PREFIX}.sorted.bam | bedtools merge -i - > output/${ASM_PREFIX}.${REF_PREFIX}.duplicated.bed
        # Project duplicated blocks onto the assembly coords
        

    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File  = glob("output/*.flagger_final.bed")[0]
        File simplifiedFinalBed = glob("output/*.flagger_final.simplified.bed")[0]
    }
}

