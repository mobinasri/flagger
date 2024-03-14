version 1.0

import "../tasks/dipcall.wdl" as dipcall_t
import "../tasks/asmgene.wdl" as asmgene_t
import "../tasks/quast.wdl" as quast_t
import "../tasks/yak_non_trio.wdl" as yak_non_trio_t

workflow standardQualityControlNonTrio {

    input {
        String sampleName
        Array[File] sampleReadsILM
        File assemblyFastaHap1
        File assemblyFastaHap2
        Boolean isMaleSample
        File referenceFasta
        File geneAnnotationFile
    }

    ### Dipcall: ###
    call dipcall_t.dipcall as dipcall {
        input:
            assemblyFastaMat=assemblyFastaHap1,
            assemblyFastaPat=assemblyFastaHap2,
            referenceFasta=referenceFasta,
            isMaleSample=isMaleSample
    }

    ### Minimap2 Gene Stats ###
    call asmgene_t.asmgene as asmgeneHap1 {
        input:
            assemblyFasta=assemblyFastaHap1,
            genesFasta=geneAnnotationFile,
            referenceFasta=referenceFasta
    }
    call asmgene_t.asmgene as asmgeneHap2 {
        input:
            assemblyFasta=assemblyFastaHap2,
            genesFasta=geneAnnotationFile,
            referenceFasta=referenceFasta
    }

    ### Quast ###
    call quast_t.quast as quastHap1 {
        input:
            assemblyFasta=assemblyFastaHap1,
            extraArguments="--large --est-ref-size 3100000000 --no-icarus"
    }
    call quast_t.quast as quastHap2 {
        input:
            assemblyFasta=assemblyFastaHap2,
            extraArguments="--large --est-ref-size 3100000000 --no-icarus"
    }


    ### Yak ###
    call yak_non_trio_t.runNonTrioYakAssemblyStats as yak {
        input:
            sampleReadsILM=sampleReadsILM,
            referenceFasta=referenceFasta,
            assemblyFastaHap1=assemblyFastaHap1,
            assemblyFastaHap2=assemblyFastaHap2
    }

    ### Consolidation ###
    call consolidate {
        input:
            sampleName = sampleName,
            dipcallFullOutput = dipcall.outputTarball,
            hap1GeneStats = asmgeneHap1.geneStats,
            hap2GeneStats = asmgeneHap2.geneStats,
            hap1QuastResults = quastHap1.outputTarball,
            hap2QuastResults = quastHap2.outputTarball,
            yakResults = yak.outputTarball
    }

	output {
        File dipcallVCF = dipcall.outputVCF
        File dipcallBED = dipcall.outputBED
        File asmgeneHap1Summary = asmgeneHap1.geneStats
        File asmgeneHap2Summary = asmgeneHap2.geneStats
        File quastHap1Summary = quastHap1.outputSummary
        File quastHap2Summary = quastHap2.outputSummary
        File yakSummary = yak.outputSummary
        File allResults = consolidate.allResults
	}
}


task consolidate {
    input{
        String sampleName
        File dipcallFullOutput
        File hap1GeneStats
        File hap2GeneStats
        File hap1QuastResults
        File hap2QuastResults
        File yakResults

        # runtime configurations
        Int memSizeGB=8
        Int threadCount=8
        Int diskSizeGB=256
        String dockerImage="tpesout/hpp_base:latest"
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

        # init
        OUT="~{sampleName}_StandardQC"
        mkdir $OUT

        # asmgene
        mkdir -p $OUT/asmgene/hap1
        mkdir -p $OUT/asmgene/hap2
        cp ~{hap1GeneStats} $OUT/asmgene/hap1/
        cp ~{hap2GeneStats} $OUT/asmgene/hap2/

        # dipcall
        mkdir $OUT/dipcall
        cd $OUT/dipcall
        tar xvf ~{dipcallFullOutput}
        mv *dipcall/* . ; rmdir *dipcall
        cd ../..

        # quast
        mkdir -p $OUT/quast/hap1
        mkdir -p $OUT/quast/hap2
        cd $OUT/quast/hap1
        tar xvf ~{hap1QuastResults}
        mv *quast/* . ; rmdir *quast
        cd ../hap2
        tar xvf ~{hap2QuastResults}
        mv *quast/* . ; rmdir *quast
        cd ../../..

        # yak
        mkdir -p $OUT/yak
        cd $OUT/yak
        tar xvf ~{yakResults}
        cd ../..

        # finalize
        tar czvf ~{sampleName}_StandardQC.tar.gz $OUT

    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
    }

    output {
        File allResults = glob("*_StandardQC.tar.gz")[0]
    }
}

