version 1.0 

workflow runCorrectBam{
    call correctBam
    output {
        File correctedBam = correctBam.correctedBam
        File correctedBamIndex = correctBam.correctedBamIndex
        File excludedReadIdsText = correctBam.excludedReadIdsText
    }
}

task correctBam {
    input {
        File? phasingLogText
        File? mapqTableText
        File bam
        String? options
        String suffix
        Boolean flagRemoveMultiplePrimary = false
        Boolean flagRemoveSupplementary = false
        # runtime configurations
        Int memSize=8
        Int threadCount=8
        Int diskSize=512
        String dockerImage="mobinasri/secphase:v0.2.0"
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

        mkdir output
        OPTIONS="~{options}"

        if [ -n "~{phasingLogText}" ]
        then
            OPTIONS="${OPTIONS} --phasingLog ~{phasingLogText}"
        fi

        if [ -n "~{mapqTableText}" ]
        then
            OPTIONS="${OPTIONS} --mapqTable ~{mapqTableText}"
        fi

        # Create an empty file in case flagRemoveMultiplePrimary and flagRemoveSupplementary both are false
        touch ${PREFIX}.excluded_read_ids.txt

        if [ -n "~{true="REMOVE" false="" flagRemoveMultiplePrimary}" ]
        then
            samtools view -F 0x904 ~{bam} | cut -f 1 | sort | uniq -c | awk '$1 > 1{print $2}' > ${PREFIX}.excluded_read_ids.txt 
        fi
        
        if [ -n "~{true="REMOVE" false="" flagRemoveSupplementary}" ]
        then
            samtools view -f 0x800 ~{bam} | cut -f 1 | sort -u >> ${PREFIX}.excluded_read_ids.txt
        fi
 
        if [ -n "~{true="REMOVE" false="" flagRemoveSupplementary || flagRemoveMultiplePrimary}" ]
        then
            OPTIONS="${OPTIONS} --exclude ${PREFIX}.excluded_read_ids.txt"
        fi

 
        correct_bam ${OPTIONS} -i ~{bam} -o output/$PREFIX.~{suffix}.bam -n~{threadCount} 
        samtools index -@~{threadCount} output/$PREFIX.~{suffix}.bam
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File correctedBam = glob("output/*.bam")[0]
        File correctedBamIndex = glob("output/*.bam.bai")[0]
        File excludedReadIdsText = glob("*.excluded_read_ids.txt")[0]
    }
}

