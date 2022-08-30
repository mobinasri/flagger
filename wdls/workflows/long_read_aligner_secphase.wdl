version 1.0

import "long_read_aligner_scattered.wdl" as aligner_t
import "../../ext/secphase/wdls/workflows/correct_bam.wdl" as correct_bam_t
import "../../ext/secphase/wdls/workflows/secphase.wdl" as secphase_t

workflow longReadFullAlignment {
    input {
        String aligner="winnowmap"
        String preset
        String sampleName
        String sampleSuffix
        Array[File] readFiles
        Int splitNumber = 16
        File assembly
        Int kmerSize = 15
        String alignerOptions="--eqx -I8g"
        String fastqOptions=""
        File? referenceFasta
        Int preemptible=2
        Int extractReadsDiskSize=512
        String zones="us-west2-a"
    }

    # Aligning reads using either minimap2 or winnowmap
    call aligner_t.longReadAlignmentScattered as readAligner{
        input:
            aligner = aligner,
            preset = preset,
            sampleName = sampleName,
            sampleSuffix = sampleSuffix,
            readFiles = readFiles,
            splitNumber = splitNumber,
            assembly = assembly,
            referenceFasta = referenceFasta,
            preemptible = 2,
            extractReadsDiskSize = extractReadsDiskSize,
            alignerOptions = alignerOptions,
            fastqOptions = fastqOptions,
            kmerSize = kmerSize,
            zones = "us-west2-a"
    }

    # Running Secphase to find the alignments that need correction
    call secphase_t.runSecPhase as secphase {
        input:
            inputBam = readAligner.bamFile,
            diploidAssemblyFastaGz = assembly,
            debugMode = false
    }

    # Running correctBam to swap primary/secondary tag for detected alignments
    call correct_bam_t.correctBam {
        input:
            phasingLogText = secphase.outLog,
            bam = readAligner.bamFile,
            options = "--minReadLen 100 --minAlignment 100 --maxDiv 0.5",
            suffix = "secphase",
            flagRemoveMultiplePrimary = false,
            flagRemoveSupplementary = false,
            # runtime configurations
            memSize=8,
            threadCount=8,
            diskSize= ceil(2 * size(readAligner.bamFile, "GB")) + 64,
            dockerImage="mobinasri/secphase:dev-v0.1",
            preemptible=2
    }

    output {
        File finalBam = correctBam.correctedBam
        File finalBamIndex = correctBam.correctedBamIndex
        File excludedReadIdsText = correctBam.excludedReadIdsText
    }

}
