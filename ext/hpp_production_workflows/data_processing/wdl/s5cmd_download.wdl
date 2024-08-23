version 1.0

workflow s3_download_s5cmd_wf {
    meta {
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
        description: "Download s3 objects that are stored in open buckets that don't charge egress with s5cmd"
    }

    parameter_meta {
        s3_uris: "Files should be named s3://blah/blah.blah"
    }

    input {
        Array[String] s3_uris
    }


    scatter (s3_uri in s3_uris) {
        call s3_download_s5cmd {
            input:
                s3_uri  = s3_uri
        }
    }


   output {
        Array[File] files = s3_download_s5cmd.file
    }  
}


task s3_download_s5cmd {
    input {
        String s3_uri

        Int memSizeGB   = 16
        Int threadCount = 4
        Int diskSize    = 250        
        Int preempts    = 2

    }


    command <<<
        set -eux -o pipefail

        mkdir download 

        s5cmd --no-sign-request cp  \
            ~{s3_uri} \
            download/

    >>>

    output {
        File file = glob("download/*")[0]
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        docker: "humanpangenomics/s5cmd@sha256:5148765f287827dca1cbcbce92075b1ea66213aaeda6113f681190d9e75796ad"
        preemptible: preempts
    }
}
