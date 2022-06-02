version 1.0

workflow misjoinCheck {

    call paftoolsMisjoinCheck 

    output {
        File misjoinSummary = paftoolsMisjoinCheck.misjoinSummary
    }
}



task paftoolsMisjoinCheck {

    input {
        File minigraphPaf
        File centromereBed
        String sampleName
        String outputFileTag = "misjoin_summary"
        String args = "-ec"
        
        Int memSizeGB   = 4
        Int diskSizeGB  = 32
        String dockerImage = "humanpangenomics/hpp_paftools:latest"
    }

    String misjoinSummaryFn = "${sampleName}.${outputFileTag}.txt"

    command <<<

        set -eux -o pipefail
        
        ## paftools call to identify misjoins
        paftools.js misjoin ~{args} ~{centromereBed} ~{minigraphPaf} > ~{misjoinSummaryFn}
        
    >>>

    output {

        File misjoinSummary  = misjoinSummaryFn
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}