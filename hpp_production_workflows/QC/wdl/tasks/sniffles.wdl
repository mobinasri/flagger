version 1.0

# This is a task level wdl workflow to run the program SNIFFLES

workflow runSniffles {

  call Sniffles
  
  output{
    File outputFile = Sniffles.vcfOut
  }
}

task Sniffles{
  input{
    File inputBam
    String SampleName
    String? outputFileTag
    
    String dockerImage = "quay.io/biocontainers/sniffles@sha256:a403144dc9aad093a6aca476ec3eea40e27d82efaba682b753e944264f5e512d" # 1.0.12--h8b12597_1
    Int memSizeGB = 128
    Int threadCount = 64
    Int diskSizeGB = 128

  }
  
 parameter_meta{
     inputBam: "PacBio and Oxford Nanopore read data. Must be in BAM format."
     SampleName: "Sample name. Will be used in output VCF file."
     outputFileTag: "Output file tag to tag files by type of read data (HiFi/Ont)."
 }
 
 String outputName = "~{SampleName}.~{outputFileTag}_sniffles.vcf"
 
  command <<<
      # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        sniffles -m ~{inputBam} -v ~{outputName} -d 500 -n -1 -s 3

  >>>
  output{
    File vcfOut = glob("*.vcf")[0]
  }
  runtime{
    memory: memSizeGB + " GB"
    cpu: threadCount
    disks: "local-disk " + diskSizeGB + " SSD"
    docker: dockerImage
    preemptible: 1
  }
}
