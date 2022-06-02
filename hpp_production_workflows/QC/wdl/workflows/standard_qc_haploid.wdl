version 1.0

import "../tasks/merqury.wdl" as merqury_t
import "../tasks/meryl.wdl" as meryl_t
import "../tasks/asmgene.wdl" as asmgene_t
import "../tasks/quast.wdl" as quast_t
import "../tasks/yak.wdl" as yak_t

workflow standardQualityControlHaploid {

    input {
        String sampleName
        Array[File] sampleReadsILM
        Array[File] maternalReadsILM
        Array[File] paternalReadsILM
        File assembly
        Boolean isMaleSample
        File referenceFasta
        File geneAnnotationFile
    }

    ### Minimap2 Gene Stats ###
    call asmgene_t.asmgene as asmgene {
        input:
            assemblyFasta=assembly,
            genesFasta=geneAnnotationFile,
            referenceFasta=referenceFasta
    }

    ### Quast ###
    call quast_t.quast as quast {
        input:
            assemblyFasta=assembly,
            extraArguments="--large --est-ref-size 3100000000 --no-icarus"
    }

    ### Merqury ###
    call meryl_t.runMeryl as meryl {
        input:
            sampleReadsILM=sampleReadsILM,
            maternalReadsILM=maternalReadsILM,
            paternalReadsILM=paternalReadsILM,
            referenceFasta=referenceFasta
    }
    call merqury_t.merqury as merqury {
        input:
            assemblyFasta=assembly,
            kmerTarball=meryl.sampleMerylDB,
            matKmerTarball=meryl.maternalHapmer,
            patKmerTarball=meryl.paternalHapmer
    }

    ### Yak ###
    call yak_t.runYakAssemblyStats as yak {
        input:
            sampleReadsILM=sampleReadsILM,
            maternalReadsILM=maternalReadsILM,
            paternalReadsILM=paternalReadsILM,
            referenceFasta=referenceFasta,
            assemblyFastaPat=assembly,
            assemblyFastaMat=assembly
    }

    ### Consolidation ###
    call consolidate {
        input:
            sampleName = sampleName,
            geneStats = asmgene.geneStats,
            quastResults = quast.outputTarball,
            merylHapmerImages = meryl.hapmerImages,
            merquryResults = merqury.outputTarball,
            yakResults = yak.outputTarball
    }

	output {
        File asmgeneSummary = asmgene.geneStats
        File quastSummary = quast.outputSummary
        File merylSampleDB = meryl.sampleMerylDB
        File merylMaternalHapmer = meryl.maternalHapmer
        File merylPaternalHapmer = meryl.paternalHapmer
        File merquryQV = merqury.QV
        File yakSummary = yak.outputSummary
        File allResults = consolidate.allResults
	}
}




task consolidate {
    input{
        String sampleName
        File geneStats
        File quastResults
        File merylHapmerImages
        File merquryResults
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
        OUT="~{sampleName}_StandardQC_Haploid"
        mkdir $OUT

        # asmgene
        mkdir $OUT/asmgene/
        cp ~{geneStats} $OUT/asmgene/

        # meryl/merqury
        mkdir $OUT/merqury
        cd $OUT/merqury
        tar xvf ~{merylHapmerImages}
        tar xvf ~{merquryResults}
        cd ../..

        # quast
        mkdir -p $OUT/quast/
        cd $OUT/quast/
        tar xvf ~{quastResults}
        mv *quast/* . ; rmdir *quast
        cd ../..

        # yak
        mkdir -p $OUT/yak
        cd $OUT/yak
        tar xvf ~{yakResults}
        cd ../..

        # finalize
        tar czvf ~{sampleName}_StandardQC_Haploid.tar.gz $OUT/
        ls -lah

    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
    }

    output {
        File allResults = glob("*_StandardQC_Haploid.tar.gz")[0]
    }
}

