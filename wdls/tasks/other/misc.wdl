version 1.1 

task createDipAsm {
    input {
        File hap1AssemblyFastaGz
        File hap2AssemblyFastaGz
        String outputName
        # runtime configurations
        Int memSize=8
        Int threadCount=4
        Int diskSize=32
        String dockerImage="mobinasri/bio_base:v0.1"
        Int preemptible=2
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace
        
        zcat ~{hap1AssemblyFastaGz} ~{hap2AssemblyFastaGz} > ~{outputName}.fa
        samtools faidx ~{outputName}.fa
        pigz -p~{threadCount} ~{outputName}.fa

    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File diploidAssemblyFastaGz = glob("*.fa.gz")[0]
        File diploidAssemblyFastaIndex = glob("*.fai")[0]
    }
}

task createFile{
     input {
        String content = ""
        File filename = "mock.txt"
        # runtime configurations
        Int memSize=2
        Int threadCount=2
        Int diskSize=2
        String dockerImage="mobinasri/bio_base:v0.1"
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        mkdir output
        echo ~{content} > output/~{filename}
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
    }
    output {
        File outFile = glob("output/*")[0] 
    }
} 

