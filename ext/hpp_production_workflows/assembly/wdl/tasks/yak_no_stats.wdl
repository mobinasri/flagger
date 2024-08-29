version 1.0

import "../../../QC/wdl/tasks/extract_reads.wdl" as extractReads_t
import "../../../QC/wdl/tasks/arithmetic.wdl" as arithmetic_t
import "../../../QC/wdl/tasks/yak.wdl" as yak_t

workflow runYak {

    input {
        Array[File] sampleReadsILM
        String sampleName
        File? referenceFasta
        Int fileExtractionDiskSizeGB = 256
        String dockerImage = "tpesout/hpp_yak:latest"
    }

    scatter (readFile in sampleReadsILM) {
        call extractReads_t.extractReads as sampleReadsExtracted {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=fileExtractionDiskSizeGB
        }
    }

    call arithmetic_t.sum as sampleReadSize {
        input:
            integers=sampleReadsExtracted.fileSizeGB
    }

    call yak_t.yakCount as yakCountSample {
        input:
            readFiles=sampleReadsExtracted.extractedRead,
            sampleName=sampleName,
            diskSizeGB=sampleReadSize.value * 2,
            dockerImage=dockerImage
    }

    output {
        File outputYak = yakCountSample.outputYak
    }

}

