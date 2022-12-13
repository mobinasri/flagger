version 1.0 

workflow runSecPhase{
    input {
        File inputBam
        File diploidAssemblyFastaGz
        String secphaseOptions = "--hifi"
        String secphaseDockerImage = "mobinasri/secphase:dev-v0.1"
        Boolean debugMode = false
    }
    call sortByName{
         input:
             bamFile = inputBam,
             diskSize = 4 * ceil(size(inputBam, "GB")) + 64
    }
    call splitByName{
         input:
             bamFile = sortByName.outputBam,
             diskSize = 2 * ceil(size(sortByName.outputBam, "GB")) + 64
    }
    scatter (splitBam in splitByName.splitBams) { 
        call secphase{
            input:
                bamFile = splitBam,
                diploidAssemblyFastaGz = diploidAssemblyFastaGz,
                debugMode = debugMode,
                options = secphaseOptions,
                dockerImage = secphaseDockerImage,
                diskSize = ceil(size(splitBam, "GB")) + 64
        }
    }
    call concatLogs as concatErrLogs{
        input:
            logs = secphase.errLog,
            filename = basename("${inputBam}", ".bam") + ".phasing_err",
            diskSize = 256
    }
    call concatLogs as concatOutLogs{
        input:
            logs = secphase.outLog,
            filename = basename("${inputBam}", ".bam") + ".phasing_out"
    }

    output{
        File errLog = concatErrLogs.log
        File outLog = concatOutLogs.log
    }
}

task secphase {
    input {
        File bamFile
        File diploidAssemblyFastaGz
        String options = "--hifi"
        Boolean debugMode
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=128
        String dockerImage="mobinasri/secphase:dev-v0.1"
        Int preemptible=2
        String zones="us-west2-a"
    }
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        ##set -o pipefail
        # cause a bash script to exit immediately when a command fails
        ##set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        ##set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        ##set -o xtrace

        BAM_FILENAME=$(basename ~{bamFile})
        BAM_PREFIX=${BAM_FILENAME%.bam}

        gunzip -c ~{diploidAssemblyFastaGz} > asm.fa
        samtools faidx asm.fa

        mkdir output
        COMMAND=~{true="secphase_debug" false="secphase" debugMode}
        ${COMMAND} ~{options} -i ~{bamFile} -f asm.fa 2> err.log > out.log || true
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones : zones
    }
    output {
        File errLog = "err.log"
        File outLog = "out.log"
    }
}

task concatLogs {
    input {
        Array[File] logs
        String filename
        # runtime configurations
        Int memSize=2
        Int threadCount=1
        Int diskSize=32
        String dockerImage="mobinasri/secphase:dev-v0.1"
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

        mkdir output
        cat ~{sep=" " logs} > output/~{filename}.txt
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones : zones
    }
    output {
        File log = glob("output/*.txt")[0]
    }
}


task splitByName {
    input {
        File bamFile
        Int NReadsPerBam = 400000
        # runtime configurations
        Int memSize=16
        Int threadCount=8
        Int diskSize=512
        String dockerImage="mobinasri/secphase:dev-v0.1"
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

        BAM_FILENAME=$(basename ~{bamFile})
        BAM_PREFIX=${BAM_FILENAME%.bam}

        mkdir output
        split_bam_by_readname -i ~{bamFile} -o output -n ~{NReadsPerBam} -t~{threadCount}
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones : zones
    }
    output {
        Array[File] splitBams = glob("output/*.bam")
    }
}


task sortByName {
    input {
        File bamFile
        String excludeSingleAlignment="yes"
        # runtime configurations
        Int memSize=16
        Int threadCount=8
        Int diskSize=1024
        String dockerImage="mobinasri/secphase:dev-v0.1"
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
        
        BAM_FILENAME=$(basename ~{bamFile})
        BAM_PREFIX=${BAM_FILENAME%.bam}
       
        mkdir output
        if [ ~{excludeSingleAlignment} == "yes" ]; then
            samtools view ~{bamFile} | cut -f1 | sort | uniq -c > readnames.txt
            cat readnames.txt | awk '$1 > 1' | tr -s ' ' | cut -d ' ' -f3 > selected_readnames.txt
            extract_reads -i ~{bamFile} -o output/${BAM_PREFIX}.bam -r selected_readnames.txt
        else
            ln ~{bamFile} output/${BAM_PREFIX}.bam
        fi
        samtools sort -n -@8 output/${BAM_PREFIX}.bam > output/${BAM_PREFIX}.sorted_by_qname.bam
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones : zones
    }
    output {
        File outputBam = glob("output/*.sorted_by_qname.bam")[0]
    }
}

task sortByContig {
    input {
        File bamFile
        # runtime configurations
        Int memSize=8
        Int threadCount=4
        Int diskSize=128
        String dockerImage="mobinasri/secphase:dev-v0.1"
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

        BAM_FILENAME=$(basename ~{bamFile})
        BAM_PREFIX=${BAM_FILENAME%.bam}

        mkdir output
        samtools sort -@~{threadCount} ~{bamFile} > output/${BAM_PREFIX}.bam
        samtools index output/${BAM_PREFIX}.bam
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones : zones
    }
    output {
        File outputBam = glob("output/*.bam")[0]
        File outputBai = glob("output/*.bai")[0]
    }
}

