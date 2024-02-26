version 1.1 

workflow runTest{
    input{
        Array[File] files = []
    }
    scatter (file in files){
        call test {
           input :
              file1 = file 
        }
    }
    output {
        Array[File] out = test.out
    }    
}
task test {
    input {
        File file1
        # runtime configurations
        Int memSize=2
        Int threadCount=2
        Int diskSize=4
        String dockerImage="mobinasri/bio_base:v0.1"
        Int preemptible=2
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace
        
        cat ~{file1} > test.txt

    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File out = glob("test.txt")[0]
    }
}

