version 1.0

import "../tasks/dipcall.wdl" as dipcall_t
import "../tasks/asmgene.wdl" as asmgene_t
import "../tasks/quast.wdl" as quast_t
import "../tasks/yak_no_qv.wdl" as yak_t

workflow standardQualityControl {

    input {
        String sampleName
        Array[File] maternalReadsILM
        Array[File] paternalReadsILM
        File paternalAssembly
        File maternalAssembly
        Boolean isMaleSample
        File referenceFasta
        File geneAnnotationFile
    }
    
    ### Dipcall (v0.3) ###
    call dipcall_t.dipcall as dipcall {
        input:
            assemblyFastaPat=paternalAssembly,
            assemblyFastaMat=maternalAssembly,
            referenceFasta=referenceFasta,
            isMaleSample=isMaleSample
    }

    ### Minimap2 Gene Stats ###
    call asmgene_t.asmgene as asmgenePaternal {
        input:
            assemblyFasta=paternalAssembly,
            genesFasta=geneAnnotationFile,
            referenceFasta=referenceFasta
    }
    call asmgene_t.asmgene as asmgeneMaternal {
        input:
            assemblyFasta=maternalAssembly,
            genesFasta=geneAnnotationFile,
            referenceFasta=referenceFasta
    }

    ### Quast ###
    call quast_t.quast as quastPaternal {
        input:
            assemblyFasta=paternalAssembly,
            extraArguments="--large --est-ref-size 3100000000 --no-icarus"
    }
    call quast_t.quast as quastMaternal {
        input:
            assemblyFasta=maternalAssembly,
            extraArguments="--large --est-ref-size 3100000000 --no-icarus"
    }

    ### Yak ###
    call yak_t.runYakAssemblyStats as yak {
        input:
            maternalReadsILM=maternalReadsILM,
            paternalReadsILM=paternalReadsILM,
            referenceFasta=referenceFasta,
            assemblyFastaPat=paternalAssembly,
            assemblyFastaMat=maternalAssembly
    }

    ### Consolidation ###
    call consolidate {
        input:
            sampleName = sampleName,
            dipcallFullOutput = dipcall.outputTarball,
            paternalGeneStats = asmgenePaternal.geneStats,
            maternalGeneStats = asmgeneMaternal.geneStats,
            paternalQuastResults = quastPaternal.outputTarball,
            maternalQuastResults = quastMaternal.outputTarball,
            yakResults = yak.outputTarball
    }

	output {
        File dipcallVCF = dipcall.outputVCF
        File dipcallBED = dipcall.outputBED
        File asmgenePaternalSummary = asmgenePaternal.geneStats
        File asmgeneMaternalSummary = asmgeneMaternal.geneStats
        File quastPaternalSummary = quastPaternal.outputSummary
        File quastMaternalSummary = quastMaternal.outputSummary
        File yakSummary = yak.outputSummary
        File allResults = consolidate.allResults
	}
}




task consolidate {
    input{
        String sampleName
        File dipcallFullOutput
        File paternalGeneStats
        File maternalGeneStats
        File paternalQuastResults
        File maternalQuastResults
        File yakResults
        # runtime configurations
        Int memSize=8
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
        mkdir -p $OUT/asmgene/mat
        mkdir -p $OUT/asmgene/pat
        cp ~{paternalGeneStats} $OUT/asmgene/pat/
        cp ~{maternalGeneStats} $OUT/asmgene/mat/

        # dipcall
        mkdir $OUT/dipcall
        cd $OUT/dipcall
        tar xvf ~{dipcallFullOutput}
        mv *dipcall/* . ; rmdir *dipcall
        cd ../..

        # quast
        mkdir -p $OUT/quast/mat
        mkdir -p $OUT/quast/pat
        cd $OUT/quast/mat
        tar xvf ~{maternalQuastResults}
        mv *quast/* . ; rmdir *quast
        cd ../pat
        tar xvf ~{paternalQuastResults}
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
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
    }

    output {
        File allResults = glob("*_StandardQC.tar.gz")[0]
    }
}

