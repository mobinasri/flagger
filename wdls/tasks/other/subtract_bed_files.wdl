version 1.1

import "bedtools.wdl" as bedtools_t

workflow unionBedFiles {
    input {
        File firstBed
        File secondBed
        String suffix
    }
    call bedtools_t.subtract {
        input:
            firstBed = firstBed,
            secondBed = secondBed,
            outputPrefix = basename("${firstBed}", ".bed") + "${suffix}"
    }
    output {
       File subtractBed = subtract.subtractBed
    }
}
