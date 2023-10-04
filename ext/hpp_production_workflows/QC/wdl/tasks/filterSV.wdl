version 1.0

# This is a task level wdl workflow to run a python script to filter SV calls

workflow runFilterSV {

    call Filter
    
    output {
        File outputFile = Filter.vcfOut
    }
}

task Filter{
    input {
        File inputVcf
        String SampleName
        String? outputFileTag
        
        String dockerImage = "kishwars/t2t_polishing@sha256:418486a1e88c48555ad4f7158c0a9923762182e7c9cd883342ffe0a161d89de6" # 0.1
        Int memSize = 128
        Int threadCount = 64
        Int diskSizeGB = 128
    }
    parameter_meta{
     inputVcf: "Reads aligned to assembly. Must be in BAM format."
     SampleName: "Sample name. Will be used in output VCF file."
     outputFileTag: "Output file tag to tag files by type of read data (HiFi/Ont)."
    }
    
    String outputName = "~{SampleName}.~{outputFileTag}_filtered.vcf"
    
    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        ln -s /opt/filter.py .
        python3 filter.py ~{inputVcf} > ~{outputName}

    >>>
    output {
        File vcfOut = glob("*filtered.vcf")[0]
    }
    runtime {
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
