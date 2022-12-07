version 1.0

import "long_read_aligner.wdl" as aligner_t
workflow asm2asmAlignment {
    input {
        String aligner="winnowmap"
        String preset
        File queryAssemblyFastaGz
        File refAssemblyFastaGz
        String suffix=""
        # splitAssembly is recommended to be true if the assembly is having almost
        # the same quality as reference (including complete centromeres) since minimap2/winnowmap tests showed
        # that it takes forever to do the alignment of such assemblies to a reference like chm13v2.0
        Boolean splitAssembly = false # split contigs before alignment
        Int splitSize = 20000000 # maximum size (in bases) of split contigs
        String zones = "us-west2-a"
    }
    if (splitAssembly){
        call runSplitAssembly{
            input:
                assemblyFastaGz = queryAssemblyFastaGz,
                splitSize = splitSize
        }
    }

    File queryAssemblyFastaGzProcessed = select_first([runSplitAssembly.splitAssemblyFastaGz, queryAssemblyFastaGz])
    ## align query assembly to the ref assembly
    call aligner_t.alignmentBam{
        input:
            aligner =  aligner,
            preset = preset,
            suffix = suffix,
            refAssembly = refAssemblyFastaGz,
            readFastq_or_queryAssembly = queryAssemblyFastaGzProcessed,
            kmerSize = 19,
            zones = zones
    }
    output {
        File sortedBamFile = alignmentBam.sortedBamFile
    }
}

task runSplitAssembly {
    input {
        File assemblyFastaGz
        Int splitSize
        # runtime configurations
        Int memSize=8
        Int threadCount=8
        Int diskSize=128
        String dockerImage="mobinasri/flagger:v0.2"
        Int preemptible=2
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        FILENAME=$(basename ~{assemblyFastaGz})
        PREFIX=${FILENAME%%.f(ast)?a.gz}

        gunzip -c ~{assemblyFastaGz} > ${PREFIX}.fa
        samtools faidx ${PREFIX}.fa
        python3 /home/programs/src/split_fai_by_length.py --fai asm.fa.fai --splitSize ~{splitSize} > ${PREFIX}.bed
        bedtools getfasta -fi ${PREFIX}.fa -bed ${PREFIX}.bed -fo ${PREFIX}.split.fa
        pigz -p8 ${PREFIX}.split.fa
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File splitAssemblyFastaGz = glob("*.split.fa.gz")
    }
}
