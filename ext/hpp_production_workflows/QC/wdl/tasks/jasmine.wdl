version 1.0

# This is a task level wdl workflow to run Jasmine

workflow runJasmine {

    call Jasmine
    
    output{
        File SV_filelist = Jasmine.SV_filelist
        File outputFile = Jasmine.vcfOut
    }
}

task Jasmine{
    input{
        Array[File] InputVCFs
        String SampleName
        
        Int? maxDist = 500
        Float? minSeqID = 0.3
        Int? specReads = 3
        
        String dockerImage = "quay.io/biocontainers/jasminesv@sha256:1b591512db1dbe32b34a09cebb285f0293683e3effb7f299a68cbad2ae14c955" # 1.1.4
        Int memSize = 128
        Int threadCount = 64
        Int diskSizeGB = 128

    }
    
    parameter_meta{
        InputVCFs: "Array of VCF files of each sample to be used."
        SampleName: "Sample name. Will be used in output VCF file."
        maxDist: "A constant integer value such that the distance threshold for every variant will be equal to this value."
        minSeqID: "The sequence identity threshold required for two insertions to be merged."
        specReads: "The minimum number of reads required to support a variant for it to be considered a specific call."
    }
    
    String SV_like_errors = "~{SampleName}_SV_like_errors.vcf"
    
    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail

        echo "${InputVCFs[*]}"

        VCFS=(~{sep=' ' InputVCFs})

        printf "%s\n" "${VCFS[@]}" > SV_filelist.txt

        jasmine max_dist=~{maxDist} min_seq_id=~{minSeqID} spec_reads=~{specReads} --output_genotypes \
        file_list=SV_filelist.txt out_file=~{SV_like_errors}

    >>>
    output{
        File SV_filelist = glob("*.txt")[0]
        File vcfOut = glob("*.vcf")[0]
    }
    runtime{
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
