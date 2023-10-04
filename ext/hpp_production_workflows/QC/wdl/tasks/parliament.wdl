version 1.0

# This is a task level wdl workflow to run the program PARLIAMENT

workflow runParliament{

    call Parliament
  
  output{
    File ParliamentVCF = Parliament.vcfOut
  }
}

task Parliament{
  input{
    File inputBam
    File indexBam
    File refGenome
    File indexGenome

    String SampleName = SampleName
    Boolean? filterShortContigs = true
    String? otherArgs = "--breakdancer --breakseq --manta --cnvnator --lumpy --delly_deletion --genotype --svviz_only_validated_candidates"

    String dockerImage = "dnanexus/parliament2@sha256:9076e0cb48f1b0703178778865a6f95df48a165fbea8d107517d36b23970a3d3" # latest
    Int memSize = 128
    Int threadCount = 64
    Int diskSizeGB = 128
  }

  parameter_meta{
    inputBam: "Illumina BAM file for which to call structural variants containing mapped reads."
    indexBam: "Corresponding index for the Illumina BAM file."
    refGenome: "Genome reference file that matches the reference used to map the Illumina inputs."
    indexGenome: "Corresponding index for the reference genome file."
    filterShortContigs: "If true, contigs shorter than 1 MB will be filtered out. Default is true. Enter false to keep short contigs."
    otherArgs: "Other optional arguments can be defined here. Refer to https://github.com/dnanexus/parliament2#help for more details."
    }
    
    String prefix = "~{SampleName}.Parl"

  command <<<
      # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        # copy input files to the /in folder to make them accessible to the parliament2.py script
        cp ~{inputBam} /home/dnanexus/in
        cp ~{refGenome} /home/dnanexus/in
        cp ~{indexBam} /home/dnanexus/in
        cp ~{indexGenome} /home/dnanexus/in

        # initialize command
        
        cmd=(python /home/dnanexus/parliament2.py --bam ~{basename(inputBam)} -r ~{basename(refGenome)})
        cmd+=( --prefix ~{prefix} --bai ~{basename(indexBam)} --fai ~{basename(indexGenome)} )
        
        # pass filter_short_contigs argument based on user input, default being true and Run PARLIAMENT
        if ["~{filterShortContigs}" == "false"];
        then
            cmd+=(~{otherArgs})
        else
            cmd+=(--filter_short_contigs ~{otherArgs})
        fi
        
        # run command

        "${cmd[@]}"


        # copy output files to output folder
        cp /home/dnanexus/out/~{prefix}.combined.genotyped.vcf .
  >>>
  output{
    File vcfOut = "~{prefix}.combined.genotyped.vcf"
  }
  runtime{
    memory: memSize
    cpu: threadCount
    disks: "local-disk " + diskSizeGB + " SSD"
    docker: dockerImage
    preemptible: 1
  }
}
