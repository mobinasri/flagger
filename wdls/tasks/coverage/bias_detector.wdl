version 1.0


workflow RunBiasDetector{
    call biasDetector
    output {
        Array[File] biasedRegionBedArray = biasDetector.biasedRegionBedArray
        Array[Float] biasedRegionFactorArray = biasDetector.biasedRegionFactorArray
    }
}



task biasDetector {
    input {
        Array[File] bedArray
        File inputBam
        File inputBai
        File assemblyFastaGz
        # runtime configurations
        Int memSize=64
        Int threadCount=16
        Int diskSize=floor(size(inputBam, "GB")) + 64
        String dockerImage="mobinasri/flagger:v0.3.4_bias_detector"
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

        # make a bed file spanning the whole genome
        gunzip -c ~{assemblyFastaGz} > asm.fa
        samtools faidx asm.fa
        cat asm.fa.fai | awk '{print $1"\t"0"\t"$2}' | bedtools sort -i - > asm.bed

        # make a bed file to be used a the baseline bed file
        cat ~{sep=" " bedArray} | bedtools sort -i - | bedtools merge -i - > all_given_beds.bed
        bedtools subtract -a asm.bed -b all_given_beds.bed > baseline.bed
       
        # make a json file pointing to all bed files
        echo "{" > bed_files.json
        for BED_FILE in ~{sep=" " bedArray};do
             BED_NAME="$(basename ${BED_FILE%%.bed})"
             echo \"${BED_NAME}\" : \"${BED_FILE}\", >> bed_files.json
        done
        echo \"baseline\":\"${PWD}/baseline.bed\" >> bed_files.json
        echo "}" >> bed_files.json
        
        ln -s ~{inputBam} alignment.bam
        ln -s ~{inputBai} alignment.bam.bai
        bias_detector -i alignment.bam -j bed_files.json -b "baseline" -t~{threadCount} > bias_table.tsv

        cat bias_table.tsv
        mkdir -p output
        
        tail -n +2 bias_table.tsv | while read line; do \
            BIAS_STATUS=$(echo "${line}" | awk '{print $2}')
            if [ ${BIAS_STATUS} == "biased" ]; then
                # save the factor values in a text file; one factor per line
                echo "${line}" | awk '{print $4+1}' >> output/factors.txt

                # save paths to the selected bed files in a text file; one path per line
                BED_FILE=$(echo ${line} | awk '{print $5}')
                BED_FILENAME="$(basename ${BED_FILE})"
                BED_PREFIX=${BED_FILENAME%.bed}
                cp ${BED_FILE} output/${BED_FILENAME}
                echo "output/${BED_FILENAME}" >> output/bed_files.txt
                echo "${BED_PREFIX}" >> output/names.txt
            fi
        done
          
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        Array[File] biasedRegionBedArray  = read_lines("output/bed_files.txt")
        Array[Float] biasedRegionFactorArray = read_lines("output/factors.txt")
        Array[String] biasedRegionNameArray = read_lines("output/names.txt")
    }
}
