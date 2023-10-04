version 1.0

workflow run_primrose {
    meta {
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
        description: "Outputs hifi bam file with 5-Methylcytosine (5mC) predictions of each CpG in PacBio HiFi reads with [primrose](https://github.com/PacificBiosciences/primrose)"
    }

    input {
        Array[File] input_bams
    }


    scatter (input_bam in input_bams) {
        call primrose {
            input:
                input_bam  = input_bam
        }
    }


   output {
        Array[File] output_5mc_bam = primrose.output_bam
    }  
}


task primrose {
    input {
        File input_bam

        Int memSize   = 64
        Int threadCount = 8
        Int addlDisk    = 50        
        Int preempts    = 3
    }

    parameter_meta {
        input_bam: "Input for primrose are PacBio HiFi reads with kinetics. File should be named *.hifi_reads.bam"
    }

    # Estimate disk size required
    Int input_bam_size   = ceil(size(input_bam, "GB"))    
    Int final_disk_dize  = input_bam_size * 2 + addlDisk


    command <<<
        set -eux -o pipefail

        FILENAME=$(basename -- "~{input_bam}" .hifi_reads.bam)

        primrose ~{input_bam} ${FILENAME}.5mc.hifi_reads.bam

    >>>

    output {
        File output_bam = glob("*.5mc.hifi_reads.bam")[0]
    }

    runtime {
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + final_disk_dize + " SSD"
        docker: "humanpangenomics/primrose@sha256:93ed60851f6a43e46e20b4075b0fea146e8ecc0b86d00b6800b987f0253fcd2e"
        preemptible: preempts
    }
}
