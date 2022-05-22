version 1.0

import "long_read_aligner.wdl" as aligner_t
workflow asm2asmAlignment {
    input {
        String aligner="winnowmap"
        String preset
        File queryAssemblyFastaGz
        File refAssemblyFastaGz
        String suffix=""
        String zones = "us-west2-a"
    }
    ## align query assembly to the ref assembly
    call aligner_t.alignmentBam{
        input:
            aligner =  aligner,
            preset = preset,
            suffix = suffix,
            refAssembly = refAssemblyFastaGz,
            readFastq_or_queryAssembly = queryAssemblyFastaGz,
            kmerSize = 19,
            diskSize = 64,
            preemptible = 2,
            zones = zones
    }
    output {
        File sortedBamFile = alignmentBam.sortedBamFile
    }
}
