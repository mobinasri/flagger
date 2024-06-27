version 1.0

workflow runLiftoff {
	call liftoff
        output{
            File outputGff3 = liftoff.outputGff3
        }
}

task liftoff {
    input {
        String sample
        String suffix
        File assemblyFastaGz
        File referenceFastaGz
        File geneGff3
        Int memSizeGB = 32
        Int threadCount = 8
        Int diskSizeGB = 128
        String dockerImage = "mobinasri/liftoff:1.6.2"
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
        cd output

        gunzip -c ~{assemblyFastaGz} > asm.fa
        gunzip -c ~{referenceFastaGz} > ref.fa

        # parameteres taken from https://www.science.org/doi/10.1126/science.abj6987 supplementary file
        liftoff -p ~{threadCount} -sc 0.95 -copies -polish -g ~{geneGff3} -o ~{sample}.~{suffix}.gff3 asm.fa ref.fa
	        
	>>>
	output {
		File outputGff3 = glob("output/*.gff3")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 2
        cpuPlatform: "Intel Cascade Lake" 
    }
}
