version 1.0

import "../../../QC/wdl/tasks/extract_reads.wdl" as extractReads_t

workflow RunDownSampling{
    input {
        Array[File] readFiles
        Float downsampledCoverage
        Float refLength = 3100000000.0
        File? referenceFasta
        Int fileExtractionDiskSizeGB = 512
    }

    scatter (readFile in readFiles){
        call extractReads_t.extractReads as extractReads {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSize=4,
                threadCount=4,
                diskSizeGB=fileExtractionDiskSizeGB,
                dockerImage="tpesout/hpp_base:latest"
        }
        call getCoverage {
            input:
                readFastq = extractReads.extractedRead,
                refLength = refLength,
                diskSize = fileExtractionDiskSizeGB
        }
    }

    call sum {
         input:
             numbers = getCoverage.coverage
    }

    scatter (readFastq in extractReads.extractedRead){ 
        call downsample {
            input:
                readFastq=readFastq,
                samplingRate = downsampledCoverage / sum.sum,
                suffix = "${downsampledCoverage}X",
                refLength = refLength,
                memSize=8,
                threadCount=8,
                diskSizeGB= fileExtractionDiskSizeGB
        }
    }

    output {
        Array[File] downsampledFastqGzs = downsample.downSampledFastqGz
        Array[Float] initialCoverages = getCoverage.coverage
        Array[Float] downsampledCoverages = downsample.downSampledCoverage
        
    }
}


task sum {
  input {
    Array[Float] numbers
    Int memSize=2
    Int threadCount=2
    Int diskSize=16
    Int preemptible=2
    String dockerImage="quay.io/masri2019/hpp_seqtk:latest"
    String zones="us-west2-a"
  }
  command <<<
  python -c "print(~{sep="+" numbers})"
  >>>
  runtime {
        docker: dockerImage
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones: zones
  }
  output {
    Float sum = read_float(stdout())
  }
}


task getCoverage {
    input{
        File readFastq
        Float refLength
        # runtime configurations
        Int memSize=4
        Int threadCount=4
        Int diskSize=256
        Int preemptible=2
        String dockerImage="quay.io/masri2019/hpp_seqtk:latest"
        String zones="us-west2-a"
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

        mkdir downsampled
        FILE_NAME=$(basename ~{readFastq})
        PREFIX=${FILE_NAME%.*}

        ln ~{readFastq} ${PREFIX}.fq
        samtools faidx ${PREFIX}.fq
        cat ${PREFIX}.fq.fai | awk -vrefLength=~{refLength} '{totalBases+=$2}END{printf "%.2f\n",totalBases/refLength}' > coverage.txt

    >>>

    runtime {
        docker: dockerImage
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones: zones
    }

    output {
        Float coverage = read_float("coverage.txt")
    }
}

task downsample {
    input{
        File readFastq
        Float samplingRate
        String suffix
        Float refLength
        # runtime configurations
        Int memSize=4
        Int threadCount=8
        Int diskSizeGB=256
        Int preemptible=2
        String dockerImage="quay.io/masri2019/hpp_seqtk:latest"
        String zones="us-west2-a"
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

        FILE_NAME=$(basename ~{readFastq})
        PREFIX=${FILE_NAME%.*}
    
        mkdir downsampled
        seqtk sample ~{readFastq} ~{samplingRate} > downsampled/${PREFIX}.~{suffix}.fq
        samtools faidx downsampled/${PREFIX}.~{suffix}.fq
        cat downsampled/${PREFIX}.~{suffix}.fq.fai | awk -v refLength=~{refLength} '{totalBases += $2}END{printf "%.2f\n", totalBases/refLength}' > cov.txt
        pigz -p8 downsampled/${PREFIX}.~{suffix}.fq > downsampled/${PREFIX}.~{suffix}.fq.gz
        
    >>>

    runtime {
        docker: dockerImage
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible : preemptible
        zones: zones
    }

    output {
        Float downSampledCoverage = read_float("cov.txt")
        File downSampledFastqGz = glob("downsampled/*.fq.gz")[0]
    }
}
