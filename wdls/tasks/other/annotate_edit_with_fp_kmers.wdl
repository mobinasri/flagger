version 1.0 

workflow runAnnotateVcfWithFPKmers{
    call annotateVCFwithFPKmers
    output{
        File inducedFPKmerBed = annotateVCFwithFPKmers.inducedFPKmerBed
        File fixedFPKmerBed = annotateVCFwithFPKmers.fixedFPKmerBed
        File editsAnnotatedInducedFPKmerVcf = annotateVCFwithFPKmers.editsAnnotatedInducedFPKmerVcf
        File editsAnnotatedFixedFPKmerVcf = annotateVCFwithFPKmers.editsAnnotatedFixedFPKmerVcf
    }
}
task annotateVCFwithFPKmers {
    input {
        File hap1PolishedToRawPaf
        File hap2PolishedToRawPaf
        File polishedMerquryTarGz
        File rawMerquryTarGz
        File deeppolisherVcfGz
        String sample
        # runtime configurations
        Int memSize=16
        Int threadCount=8
        Int diskSize=64
        String dockerImage="mobinasri/flagger@sha256:5d738412b56bac5a64227569c1d6e57e7920e3d3e5724c17ab233f92279bcff6"
        Int preemptible=2
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace
        
        # extract merqury bed files for the polished assembly
        mkdir polished && cd polished
        cp ~{polishedMerquryTarGz} .
        tar -xvzf *.tar.gz
        cat *.asm_only.bed *.altHap_only.bed | bedtools sort -i - | bedtools merge -i - -c 1 -o count > ../polished_merged_fp_blocks_with_counts.bed
        
        # extract merqury bed files for the raw assembly
        cd ../ && mkdir raw && cd raw
        cp ~{rawMerquryTarGz} .
        tar -xvzf *.tar.gz
        cat *.asm_only.bed *.altHap_only.bed | bedtools sort -i - | bedtools merge -i - -c 1 -o count > ../raw_merged_fp_blocks_with_counts.bed

        cd ../

        # concat hap1 and hap2 paf files
        cat ~{hap1PolishedToRawPaf} ~{hap2PolishedToRawPaf} > dip_polished_to_raw.paf
       
        # project polished blocks to raw coordinates using the merged paf file
        python3 /home/programs/src/project_blocks_multi_thread.py \
            --mode asm2ref \
            --paf dip_polished_to_raw.paf \
            --blocks polished_merged_fp_blocks_with_counts.bed \
            --outputProjectable polished_merged_fp_blocks_with_counts.projectable.bed \
            --outputProjection polished_merged_fp_blocks_with_counts.projection_to_raw.bed \
            --threads ~{threadCount}

        bedtools sort -i polished_merged_fp_blocks_with_counts.projection_to_raw.bed > polished_merged_fp_blocks_with_counts.projection_to_raw.sorted.bed
        
        # get FP kmer blocks induced by polishing
        bedtools subtract -a polished_merged_fp_blocks_with_counts.projection_to_raw.sorted.bed \
                          -b raw_merged_fp_blocks_with_counts.bed -A > ~{sample}.induced_fp_kmer_blocks.bed
        
        # get FP kmer blocks fixed by polishing
        bedtools subtract -a raw_merged_fp_blocks_with_counts.bed \
                          -b polished_merged_fp_blocks_with_counts.projection_to_raw.sorted.bed -A > ~{sample}.fixed_fp_kmer_blocks.bed

        # get annotated vcf files, last 4 columns would be the FP kmer block coordinates along with the count of kmers
        bedtools intersect \
                   -header \
                   -a ~{deeppolisherVcfGz} \
                   -b ~{sample}.induced_fp_kmer_blocks.bed \
                   -wb > ~{sample}.deeppolisher.annotated_induced_fp_kmers.vcf
        bedtools intersect \
                   -header \
                   -a ~{deeppolisherVcfGz} \
                   -b ~{sample}.fixed_fp_kmer_blocks.bed \
                   -wb > ~{sample}.deeppolisher.annotated_fixed_fp_kmers.vcf
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File inducedFPKmerBed = glob("*.induced_fp_kmer_blocks.bed")[0]
        File fixedFPKmerBed = glob("*.fixed_fp_kmer_blocks.bed")[0]
        File editsAnnotatedInducedFPKmerVcf = glob("*_fixed_fp_kmers.vcf")[0]
        File editsAnnotatedFixedFPKmerVcf = glob("*_induced_fp_kmers.vcf")[0]
    }
}

