version 1.0

# This is a task level wdl workflow to run IRIS to filter HIFI SV calls

workflow runIris {

    call Iris
    
    output{
        File outputFile = Iris.vcfOut
    }
}

task Iris{
    input{
        File genomeIn
        File readsIn
        File vcfIn
        String SampleName
        String? IrisOut = "IRISOutputDir"
        
        String dockerImage = "quay.io/biocontainers/irissv@sha256:e854b554b11377b9b47f32d1c33b13d84b9fde2c5b99045d946f1e01568ec6a1" # 1.0.4--hdfd78af_2
        Int memSizeGB = 128
        Int threadCount = 64
        Int diskSizeGB = 128

    }
    parameter_meta{
        genomeIn: "FASTA file containing the reference genome."
        readsIn: "BAM file containing the reads."
        vcfIn: "VCF file with variant calls/supporting reads determined by Sniffles."
        SampleName: "Sample name. Will be used in output VCF file."
        IrisOut: "Name of the output directory to store output VCF file produced. If not provided, default value will be used."
    }
    
    String vcfOut = "~{SampleName}_iris.vcf"
    
    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        iris --hifi --keep_long_variants --keep_files genome_in=~{genomeIn} reads_in=~{readsIn} \
        vcf_in=~{vcfIn} vcf_out=~{vcfOut} out_dir=~{IrisOut}

    >>>
    output{
        File vcfOut = glob("*iris.vcf")[0]
    }
    runtime{
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
