version 1.0

workflow hic_group_xy_wf {
    
    meta {
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
        description: "Partitions sex chromosomes in HiC assemblies based on chrX and chrY kmers. See yak's documentation."
    }

    call group_xy

    output {
        File outputHap1FastaGz = group_xy.outputHap1FastaGz
        File outputHap2FastaGz = group_xy.outputHap2FastaGz
    }
}


task group_xy {

    input {
        File hap1_gz
        File hap2_gz
        String childID
        Boolean isMaleSample

        File chrY_no_par_yak
        File chrX_no_par_yak
        File par_yak

        Int threadCount = 16
        Int memSizeGB  = 32
        Int diskSizeGB = 96
        String dockerImage = "humanpangenomics/hifiasm@sha256:d6851d4686b7183b0364274c79311284fb5ab91c957addcc309bb599d42c495f"
    }

    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace

        ## Male samples have inconsistent sex chrom assignment. Use Yak to partition properly
        ## based on chrX and chrY yaks...
        if [[ ~{isMaleSample} == true ]]; then
        
            yak sexchr \
                -K2g \
                -t16 \
                ~{chrY_no_par_yak} \
                ~{chrX_no_par_yak} \
                ~{par_yak} \
                ~{hap1_gz} \
                ~{hap2_gz} \
                > cnt.txt

            groupxy.pl \
                cnt.txt \
                | awk '$4==1' | cut -f2 \
                    | seqtk subseq -l60 \
                    <(zcat ~{hap1_gz} ~{hap2_gz}) - \
                    | pigz \
                    > ~{childID}.hap1.corrected.fa.gz

            groupxy.pl \
                cnt.txt | \
                awk '$4==2' | cut -f2 \
                    | seqtk subseq -l60 \
                    <(zcat ~{hap1_gz} ~{hap2_gz}) - \
                    | pigz \
                    > ~{childID}.hap2.corrected.fa.gz

            ## rename (in case of input/output name conflicts)
            mv ~{childID}.hap1.corrected.fa.gz ~{childID}.hap1.xygrouped.fa.gz
            mv ~{childID}.hap2.corrected.fa.gz ~{childID}.hap2.xygrouped.fa.gz


        ## No need to do anything for female samples. Just rename and exit.
        else 
            cp ~{hap1_gz} ~{childID}.hap1.xygrouped.fa.gz
            cp ~{hap2_gz} ~{childID}.hap2.xygrouped.fa.gz
        fi
    >>>

    output {

        File outputHap1FastaGz = "~{childID}.hap1.xygrouped.fa.gz"
        File outputHap2FastaGz = "~{childID}.hap2.xygrouped.fa.gz"
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        cpu: threadCount
        docker: dockerImage
        preemptible: 1
    }
}