version 1.0

import "bedtools.wdl" as bedtools_t

workflow exclusionStepMale {
    input {
        File autosome_nonCntr_bed
        File sex_nonCntr_bed
        File autosome_cntr_bed
        File sex_cntr_bed
    }
    call bedtools_t.union as union_1 {
       input:
           bedFiles = [sex_cntr_bed, autosome_cntr_bed, sex_nonCntr_bed]
    }
    call bedtools_t.union as union_2 {
       input:
           bedFiles = [sex_cntr_bed, autosome_cntr_bed]
    }
    call bedtools_t.subtract as subtract_1 {
       input:
           firstBed = autosome_nonCntr_bed,
           secondBed = union_1.unionBed,
           outputPrefix = basename("${autosome_nonCntr_bed}", ".bed") + ".excluded"
    }
    call bedtools_t.subtract as subtract_2 {
       input:
           firstBed = sex_nonCntr_bed,
           secondBed = union_2.unionBed,
           outputPrefix = basename("${sex_nonCntr_bed}", ".bed") + ".excluded"
    }
    call bedtools_t.subtract as subtract_3 {
       input:
           firstBed = autosome_cntr_bed,
           secondBed = sex_cntr_bed,
           outputPrefix = basename("${autosome_cntr_bed}", ".bed") + ".excluded"
    }
    output {
       File autosome_nonCntr_excluded_bed = subtract_1.subtractBed
       File sex_nonCntr_excluded_bed = subtract_2.subtractBed
       File autosome_cntr_excluded_bed = subtract_3.subtractBed
    }
}
