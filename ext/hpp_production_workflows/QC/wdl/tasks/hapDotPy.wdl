version 1.0

workflow runhapDotPy {

    call hapDotPy

    output{
        File hapDotPyVCF    = hapDotPy.vcfOut
        File hapDotPyVCFIdx = hapDotPy.vcfIdxOut
        File hapDotPyTar    = hapDotPy.happyTar
    }
}

task hapDotPy{
    input{
        File truthVCF
        File queryVCF
        File assembly
        File assemblyIndex
        String sample

        Boolean passOnly = true

        Int memSize = 24
        Int threadCount = 8
        Int diskSizeGB = 128
        String dockerImage = "jmcdani20/hap.py@sha256:0812d37e7210011e407914deb3da2094ec8258077b21d3211e694e5ca303b489" # v0.3.12


    }


    String outputPrefix = "~{sample}_happy"

    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        ## Soft link fasta and index so they are in the same directory
        REF=$(basename ~{assembly})
        REF_IDX=$(basename ~{assemblyIndex}) 

        ln -s ~{assembly} ./$REF
        ln -s ~{assemblyIndex} ./$REF_IDX

        
        ## Pass argument if callRegions is set, if not just pass empty string
        if [[ ~{passOnly} ]]
        then
            PASS_ONLY_TOKEN="--pass-only"
        else
            PASS_ONLY_TOKEN=""
        fi


        ## make directory to put output into
        mkdir happy_out

        ## Run hapDotPy
        /opt/hap.py/bin/hap.py \
            ~{truthVCF} \
            ~{queryVCF} \
            -r $REF \
            -o happy_out/~{outputPrefix} \
            $PASS_ONLY_TOKEN \
            --engine=vcfeval \
            --threads=~{threadCount}

        tar czvf ~{sample}_happy.tar.gz happy_out/
    >>>
    output{
        File vcfOut     = "happy_out/~{outputPrefix}.vcf.gz"
        File vcfIdxOut  = "happy_out/~{outputPrefix}.vcf.gz.tbi"
        File happyTar   = "~{sample}_happy.tar.gz"
    }

    runtime{
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}