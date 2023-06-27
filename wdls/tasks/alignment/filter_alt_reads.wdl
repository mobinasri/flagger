version 1.0 

workflow runFilterAltReads {
    call filterAltReads
    output {
        File filteredBam = filterAltReads.filteredBam
        File altBam = filterAltReads.altBam
        File filteredBamIdex = filterAltReads.filteredBamIndex
        File altBamIndex = filterAltReads.altBamIndex
    }
}

task filterAltReads {
    input {
        File vcf
        File bam
        String moreOptions
        Float vafCutoff = 0.5
        Int qCutoff = 10
        # runtime configurations
        Int memSize=8
        Int threadCount=8
        Int diskSize=512
        String dockerImage="mobinasri/flagger:v0.3.1"
        Int preemptible=2
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
        
        FILENAME=$(basename ~{bam})
        PREFIX=${FILENAME%.bam}
        bcftools view -Ov -f PASS -m2 -M2 -v snps ~{vcf} -e 'FORMAT/VAF<~{vafCutoff} | FORMAT/GQ<~{qCutoff}' > bi_snps.passed.vcf
        mkdir output
        filter_alt_reads -i ~{bam} -o output/$PREFIX.alt_filtered.bam -f output/$PREFIX.alt.bam -v bi_snps.passed.vcf -t~{threadCount} ~{moreOptions}
        samtools index -@{threadCount} output/$PREFIX.alt_filtered.bam
        samtools index -@{threadCount} output/$PREFIX.alt.bam
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File filteredBam = glob("output/*.alt_filtered.bam")[0]
        File altBam = glob("output/*.alt.bam")[0]
        File filteredBamIndex = glob("output/*.alt_filtered.bam.bai")[0]
        File altBamIndex = glob("output/*.alt.bam.bai")[0]
    }
}

