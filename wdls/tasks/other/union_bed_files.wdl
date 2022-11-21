version 1.0

import "bedtools.wdl" as bedtools_t

workflow unionBedFiles {
    input{
        Array[File] bedFiles
        String sample 
        String prefix
    }
    call bedtools_t.union{
        input:
            bedFiles = bedFiles,
            outputPrefix = "${sample}.${prefix}"
    }
    output {
       File unionBed = union.unionBed
    }
}
