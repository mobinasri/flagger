version 1.0

workflow VariantCalling{
    input {
        File assemblyFastaGz
        File bam
        File bamIndex
        Int minMAPQ
        Int ploidy = 2
        Int numberOfCallerNodes=16
        Int nodeThreadCount=16
        String sampleName
        String sampleSuffix
        String platform
    }
    call splitBam {
        input:
            assemblyFastaGz = assemblyFastaGz,
            bam = bam,
            bamIndex = bamIndex,
            splitNumber = numberOfCallerNodes,
            threadCount = numberOfCallerNodes
    }
    scatter (part in zip(splitBam.splitBams, splitBam.splitBeds)) {
        call callVariant{
            input:
                assemblyFastaGz = assemblyFastaGz,
                bam = part.left,
                bed = part.right,
                ploidy = ploidy,
                minMAPQ = minMAPQ
        }
    }
    call mergeVcf{
        input:
            vcfGzFiles = callVariant.vcfGz,
            outputName = "${sampleName}.${sampleSuffix}.${platform}"
    }
    output{
        File vcfGz = mergeVcf.vcfGz 
    }
}

task splitBam{
    input{
        File assemblyFastaGz
        File bam
        File bamIndex
        Int splitNumber
        Int memSize=32
        Int threadCount
        Int diskSize=512
        String dockerImage="mobinasri/bio_base:v0.1"
        Int preemptible=2
        String zones="us-west2-a"
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

        ## unzip fasta file and produce its index file        
        ASSEMBLY_NAME=$(basename ~{assemblyFastaGz})
        ASSEMBLY_PREFIX=${ASSEMBLY_NAME%%.fa.gz}
        gunzip -c ~{assemblyFastaGz} > ${ASSEMBLY_PREFIX}.fa
        samtools faidx ${ASSEMBLY_PREFIX}.fa

        ## hard link the bam and bai files to the working directory
        BAM_NAME=$(basename ~{bam})
        BAM_PREFIX=${BAM_NAME%%.bam}
        ln -f ~{bam} > ${BAM_PREFIX}.bam
        ln -f ~{bamIndex} > ${BAM_PREFIX}.bam.bai

        ## make a bed file that covers the whole assembly
        cat ${ASSEMBLY_PREFIX}.fa.fai | awk '{print $1"\t"0"\t"$2}' > ${ASSEMBLY_PREFIX}.bed

        ## split the bed file of the whole assembly into multiple bed files
        mkdir split_beds split_bams
        python3 ${SPLIT_BED_PY} --bed ${ASSEMBLY_PREFIX}.bed --n 16 --dir split_beds --prefix ${ASSEMBLY_PREFIX}

        ## make a separate bam for each bed file
        seq 1 ~{splitNumber} | xargs -I {} -n 1 -P ~{threadCount} sh -c "samtools view -h -b -L split_beds/${ASSEMBLY_PREFIX}_{}.bed ${BAM_PREFIX}.bam > split_bams/${BAM_PREFIX}_{}.bam"
    >>>
    runtime {
        docker: dockerImage
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones: zones
    }
    output {
        Array[File] splitBams = glob("split_bams/*.bam")
        Array[File] splitBeds = glob("split_beds/*.bed")
    }
}

task callVariant{
    input{
        File bam
        File assemblyFastaGz
        File bed
        Int ploidy
        Int minMAPQ
        String extraArgs=""
        # runtime configurations
        Int memSize=32
        Int threadCount=16
        Int diskSize=256
        String dockerImage="mobinasri/bio_base:v0.1"
        Int preemptible=2
        String zones="us-west2-a"
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
        
        ## hard link the bam file to the working directory and produce its index file
        BAM_NAME=$(basename ~{bam})
        BAM_PREFIX=${BAM_NAME%%.bam}
        ln -f ~{bam} > ${BAM_PREFIX}.bam
        samtools index ${BAM_PREFIX}.bam

        ## unzip the fasta file and produce its index
        gunzip -c ~{assemblyFastaGz} > asm.fa
        samtools faidx asm.fa

        ## split the previously split bed file again into smaller chunks for the sake of parallelism
        mkdir split_beds
        python3 ${SPLIT_BED_PY} --bed ~{bed} --n ~{threadCount} --dir split_beds --prefix tmp
        
        ## run the variant caller for each small temporary bed file
        mkdir vcf_files
        seq 1 ~{threadCount} | xargs -I {} -n 1 -P ~{threadCount} sh -c "bcftools mpileup -a FORMAT/AD -a INFO/AD ~{extraArgs} -q~{minMAPQ} -B -d 100000 -f asm.fa -R split_beds/tmp_{}.bed ${BAM_PREFIX}.bam | bcftools call -cv --ploidy ~{ploidy} -Oz -o vcf_files/tmp.{}.vcf.gz"
        
        ## make a header for the merged vcf file
        zcat vcf_files/tmp.1.vcf.gz | awk 'substr($0,1,1) == "#"' | awk -v colname="${BAM_PREFIX}" '{if ($1 == "#CHROM"){ for(i =1; i < 10; i++){printf("%s\t",$i)}; printf("%s\n",colname)} else {print $0}}' > merged.vcf
        zcat vcf_files/*.vcf.gz | awk 'substr($0,1,1) != "#"' >> merged.vcf
        
        ## sort the merged vcf file and produce the final gzipped vcf file
        mkdir final_vcf
        bcftools sort -o final_vcf/${BAM_PREFIX}.vcf.gz -Oz merged.vcf
    >>>
    runtime {
        docker: dockerImage
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones: zones
    }
    output {
        File vcfGz = glob("final_vcf/*.vcf.gz")[0]
    }
}


task mergeVcf{
    input{
        Array[File] vcfGzFiles
        String outputName
        # runtime configurations
        Int memSize=32
        Int threadCount=16
        Int diskSize=512
        String dockerImage="mobinasri/bio_base:v0.1"
        Int preemptible=2
        String zones="us-west2-a"
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

        mkdir vcf_files
        ln ~{sep=" " vcfGzFiles} vcf_files
        files=(vcf_files/*) 

        ## make a header for the merged vcf file
        zcat ${files[0]} | awk 'substr($0,1,1) == "#"' | awk -v colname="~{outputName}" '{if ($1 == "#CHROM"){ for(i =1; i < 10; i++){printf("%s\t",$i)}; printf("%s\n",colname)} else {print $0}}' > merged.vcf
        zcat vcf_files/*.vcf.gz | awk 'substr($0,1,1) != "#"' >> merged.vcf

        ## sort the merged vcf file and produce the final gzipped vcf file
        mkdir final_vcf
        bcftools sort -o final_vcf/~{outputName}.vcf.gz -Oz merged.vcf
    >>>
    runtime {
        docker: dockerImage
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones: zones
    }
    output {
        File vcfGz = "final_vcf/~{outputName}.vcf.gz"
    }

}

