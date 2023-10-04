version 1.0 

workflow runBreakIntoContigs{
    call breakIntoContigs
    output {
        File assemblyContigsFaGz = breakIntoContigs.assemblyContigsFaGz
    }
}

task breakIntoContigs {
    input {
        File assemblyFaGz
        # runtime configurations
        Int memSize=8
        Int threadCount=4
        Int diskSize=128
        String dockerImage="quay.io/masri2019/hpp_hifiasm:0.18.2"
        Int preemptible=2
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace
        
        FILENAME=$(basename ~{assemblyFaGz})
        PREFIX=${FILENAME%%.fa.gz}

        mkdir output
        python3 /home/scripts/break_into_contigs.py --inputFasta <(zcat ~{assemblyFaGz}) | pigz -p8 - > output/${PREFIX}.contigs.fa.gz
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File assemblyContigsFaGz = glob("output/*.fa.gz")[0]
    }
}

