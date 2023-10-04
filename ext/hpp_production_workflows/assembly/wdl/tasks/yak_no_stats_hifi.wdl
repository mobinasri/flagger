version 1.0

import "../../../QC/wdl/tasks/extract_reads.wdl" as extractReads_t
import "../../../QC/wdl/tasks/arithmetic.wdl" as arithmetic_t
import "../../../QC/wdl/tasks/yak.wdl" as yak_t
import "filter_hifi_adapter.wdl" as filterHiFiAdapter_t

workflow runYak {

    input {
        Array[File] sampleReadsHiFi
        String sampleName
        File? referenceFasta
        Int fileExtractionDiskSizeGB = 256
        String dockerImage = "tpesout/hpp_yak:latest"
    }

    scatter (readFile in sampleReadsHiFi) {
        call extractReads_t.extractReads as sampleReadsExtracted {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSize=4,
                threadCount=4,
                diskSizeGB=fileExtractionDiskSizeGB,
                dockerImage=dockerImage
        }
        call filterHiFiAdapter_t.cutadapt {
            input:
                readFastq = sampleReadsExtracted.extractedRead
        }
    }

    call arithmetic_t.sum as sampleReadSize {
        input:
            integers=cutadapt.fileSizeGB
    }

    call yak_t.yakCount as yakCountSample {
        input:
            readFiles=cutadapt.filteredReadFastq,
            sampleName=sampleName,
            diskSizeGB=sampleReadSize.value * 2,
            dockerImage=dockerImage
    }

    output {
        File outputYak = yakCountSample.outputYak
    }

}

