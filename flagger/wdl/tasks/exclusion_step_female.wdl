version 1.0

import "bedtools.wdl" as bedtools_t

workflow exclusionStepFemale {
    input {
        File autosome_nonCntr_bed
        File autosome_cntr_bed
    }
    call bedtools_t.subtract {
       input:
           firstBed = autosome_nonCntr_bed,
           secondBed = autosome_cntr_bed,
           outputPrefix = basename("${autosome_nonCntr_bed}", ".bed") + ".excluded"
    }
    output {
       File autosome_nonCntr_excluded_bed = subtract.subtractBed
    }
}
