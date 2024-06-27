version 1.0 

workflow RunShortQC {
    input {
        File hap1Fasta   # paternal for trios
        File hap2Fasta   # maternal for trios
        File childYak
        File genesFasta
        File hs38Paf
        String childID

        File? patYak
        File? matYak

        # runtime configurations
        Int memSize
        Int threadCount
        String dockerImage = "humanpangenomics/hpp_qc_stats@sha256:6a64ac0be88ce9ca760eb7713922f65e66f9d09076b9d17c2c416e5d558bf1d0"
        Int preemptible=1
        Int diskSize
        String zones="genomics.default-zones"
    }

    parameter_meta {
        hap1Fasta: "(Gzipped) assembly to run through QC; paternal is passed as hap1 when trio yaks are available"
        hap2Fasta: "(Gzipped) assembly to run through QC; maternal is passed as hap2 when trio yaks are available"
        childYak: "yak file for child illumina or HiFi reads"
        genesFasta: "from gs://hifiasm/Homo_sapiens.GRCh38.cdna.all.fa"
        hs38Paf: "from gs://hifiasm/hs38.paf"
        childID: "ID of assembled sample; used to name output file"
        patYak: "(optional) yak file for paternal illumina or HiFi reads"
        matYak: "(optional) yak file for maternal illumina or HiFi reads"
    }

    call shortQC{
        input:
            hap1Fasta = hap1Fasta ,
            hap2Fasta = hap2Fasta ,
            childYak = childYak,
            patYak = patYak ,
            matYak = matYak ,
            genesFasta = genesFasta ,
            hs38Paf = hs38Paf ,
            childID = childID ,
            # runtime configurations ,
            memSize = memSize ,
            threadCount = threadCount ,
            dockerImage = dockerImage ,
            preemptible = preemptible,
            diskSize = diskSize,
            zones = zones
    }

    output {
        File qcStatsText = shortQC.qcStatsText
    }
}

task shortQC {
    input{
        File hap1Fasta   # paternal for trios
        File hap2Fasta   # maternal for trios
        File childYak
        File genesFasta
        File hs38Paf
        String childID

        File? patYak
        File? matYak
        
        # runtime configurations    
        Int memSize
        Int threadCount
        String dockerImage
        Int preemptible
        Int diskSize
        String zones
    }
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        # Computing length statistics
        k8 ${CAL_N50_PATH} ~{hap1Fasta} > hap1.len.stats.txt
        k8 ${CAL_N50_PATH} ~{hap2Fasta} > hap2.len.stats.txt

        # Aligning genes to assembly
        minimap2 -cx splice:hq -t ~{threadCount} ~{hap1Fasta} ~{genesFasta} > hap1.paf
        minimap2 -cx splice:hq -t ~{threadCount} ~{hap2Fasta} ~{genesFasta} > hap2.paf

        # hs38.ensembl.v99.cdna is already aligned to hg38 and its paf file should be given as an input
        # link to hs38.ensembl.v99.cdna : ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
        ln -s ~{hs38Paf} hs38.paf
        # Computing statistics for gene completeness
        paftools.js asmgene -a hs38.paf hap1.paf > hap1.gene.stats.txt
        paftools.js asmgene -a hs38.paf hap2.paf > hap2.gene.stats.txt

        # Computing QV stats
        yak qv -t 32 -p -K 3.2g -l 100k ~{childYak} <(zcat ~{hap1Fasta}) > pat.yak_qv.txt
        yak qv -t 32 -p -K 3.2g -l 100k ~{childYak} <(zcat ~{hap2Fasta}) > mat.yak_qv.txt
        

        if [[ -f "~{patYak}" ]] ; then
            # Compute error rates (only if we have trio info)
            yak trioeval -t ~{threadCount} ~{patYak} ~{matYak} ~{hap1Fasta} > pat.error.stats.txt
            yak trioeval -t ~{threadCount} ~{patYak} ~{matYak} ~{hap2Fasta} > mat.error.stats.txt

            # Merge all statistics in a single well-organized text file
            python3 ${QC_STATS_GENERATOR_PATH} \
                --childID ~{childID} \
                --hap1LenStats hap1.len.stats.txt \
                --hap2LenStats hap2.len.stats.txt \
                --hap1GeneStats hap1.gene.stats.txt \
                --hap2GeneStats hap2.gene.stats.txt \
                --patErrorStats pat.error.stats.txt \
                --matErrorStats mat.error.stats.txt \
                --hap1QVStats pat.yak_qv.txt \
                --hap2QVStats mat.yak_qv.txt \
                --output ~{childID}.qc.stats.txt

        else
            # Merge all statistics in a single well-organized text file
            # Don't include hamming/switch error statistics in call for non-trios!
            python3 ${QC_STATS_GENERATOR_PATH} \
                --childID ~{childID} \
                --hap1LenStats hap1.len.stats.txt \
                --hap2LenStats hap2.len.stats.txt \
                --hap1GeneStats hap1.gene.stats.txt \
                --hap2GeneStats hap2.gene.stats.txt \
                --hap1QVStats pat.yak_qv.txt \
                --hap2QVStats mat.yak_qv.txt \
                --output ~{childID}.qc.stats.txt
        fi 

    >>>

    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible: preemptible
        zones: zones
    }

    output {
        File qcStatsText = "~{childID}.qc.stats.txt"
    }
}
