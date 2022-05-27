version 1.0

import "bedtools.wdl" as bedtools_t

workflow unionBedFiles {
    call bedtools_t.union
    output {
       File unionBed = union.unionBed
    }
}
