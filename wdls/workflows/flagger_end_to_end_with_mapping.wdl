version 1.0

import "../tasks/other/misc.wdl" as misc_t
import "flagger_end_to_end.wdl" as flagger_t
import "long_read_aligner_scattered.wdl" as aligner_t

workflow FlaggerEndToEndWithMapping{
    meta {
        author: "Mobin Asri"
        email: "masri@ucsc.edu"
        description: "Mapping reads and running Flagger for evaluating a diploid assembly using long read alignments.More information at https://github.com/mobinasri/flagger"
    }
    parameter_meta {
        sampleName: "Sample name for example 'HG002'"
        suffixForFlagger: "Suffix string that contains information about this analysis for example 'hifi_winnowmap_flagger_for_hprc'"
        suffixForMapping: "Suffix string that contains information about this alignment. It will be appended to the name of the final alignment. For example 'hifi_winnowmap_v2.03_hprc_y2'"
        hap1AssemblyFasta: "Path to uncompressed or gzip-compressed fasta file of the 1st haplotype."
        hap2AssemblyFasta: "Path to uncompressed or gzip-compressed fasta file of the 2nd haplotype."
        readfiles: "An array of read files. Their format can be either fastq, fq, fastq.gz, fq.gz, bam or cram. For cram format referenceFastaForReadExtraction should also be passed."
        aligner: "Name of the aligner. It can be either minimap2, winnowmap or veritymap. (Default = winnowmap)"
        preset: "Paremeter preset should be selected based on aligner and sequencing platform. Common presets are map-pb/map-hifi/map-ont for minimap2, map-pb/map-ont for winnowmap and hifi-haploid/hifi-haploid-complete/hifi-diploid/ont-haploid-complete for veritymap"
        kmerSize: "The kmer size for using minimap2 or winnowmap. With winnowmap kmer size should be 15 and with minimap2 kmer size should be 17 and 19 for using the presets map-ont and map-hifi/map-pb respectively."
        alignerOptions: "Aligner options. It can be something like '--eqx --cs -Y -L -y' for minimap2/winnowmap. Note that if assembly is diploid and aligner is either minimap2 or winnowmap '-I8g' is necessary. If the reads contain modification tags and these tags are supposed to be present in the final alignment file, alignerOptions should contain '-y' and the aligner should be either minimap2 or winnowmap. If running secphase is enabled it is recommended to add '-p0.5' to alignerOptions; it will keep more secondary alignments so secphase will have more candidates per read. For veritymap '--careful' can be used but not recommended for whole-genome assembly since it increases the runtime dramatically."
        readExtractionOptions: "The options to be used while converting bam to fastq with samtools fastq. If the reads contain epigentic modification tags it is necessary to use '-TMm,Ml,MM,ML'"
        referenceFastaForReadExtraction: "If reads are in CRAM format then the related reference should be passed to this parameter for read extraction."
        enableAddingMDTag: "If true it will call samtools calmd to add MD tags to the final bam file. [Default = true]"
        splitNumber: "The number of chunks which the input reads should be equally split into. Note that enableSplittingReadsEqually should be set to true if user wants to split reads into equally sized chunks. [Default = 16]"
        enableSplittingReadsEqually: "If true it will merge all reads together and then split them into multiple chunks of roughly equal size. Each chunk will then be aligned via a separate task. This feature is useful for running alignment on cloud/slurm systems where there are  multiple nodes available with enough computing power and having alignments tasks distributed among small nodes is more efficient or cheaper than running a single alignment task in a large node. If the  whole workflow is being on a single node it is not recommened to use this feature since mergin and splitting reads takes its own time. [Default = false]"
        minReadLengthForMapping: "If it is greater than zero, a task will be executed for filtering reads shorter than this value before alignment. [Default = 0]"
        alignerThreadCount : "The number of threads for mapping in each alignment task [Default = 16]"
        alignerMemSize : "The size of the memory in Gb for mapping in each alignment task [Default = 48]"
        alignerDockerImage : "The mapping docker image [Default = 'mobinasri/long_read_aligner:v0.4.0']"
        correctBamOptions: "Options for the correct_bam program that can filters short/highly divergent alignments [Default = '--primaryOnly --minReadLen 5000 --minAlignment 5000 --maxDiv 0.1']"
        preemptible: "Number of retries to use preemptible nodes on Terra/GCP. [Default = 2]"
        zones: "Name of the zone for taking nodes on Terra/GCP. (Default = us-west2-a)"
        maxReadDivergenceForFlagger: "Alignments with gap-compressed ratio higher than this will be filtered in the pre-process step of flagger. (Default: 0.1)"
        potentialBiasesBedArray: "Array of bed files each of which contains regions with potential coverage bias for example one bed file can contain HSat2 regions in haplotype 1. (Default: [])"
        sexBed: "Optional bed file containing regions assigned to X/Y chromosomes. (can be either in ref or asm coordinates)"
        SDBed: "Optional Bed file containing Segmental Duplications. (can be either in ref or asm coordinates)"
        cntrBed: "Optional Bed file containing peri/centromeric satellites (ASat, HSat, bSat, gSat) without 'ct' blocks."
        cntrCtBed: "Optional Bed file containing centromere transition 'ct' blocks."
        additionalStratificationBedArray: "Array of additional stratification bed files for final stats tsv file. (Default: [])"
        additionalStratificationNameArray: "Array of names for the stratifications provided in the argument additionalStratificationBedArray. (Default: [])"
        enableProjectingBedsFromRef2Asm: "If True it means that the given bed files are in ref coors (e.g. chm13v2) and they have to be projected to asm coors. (Default: false)"
        projectionReferenceFastaGz: "The given bed files are in the coordinates of this reference. A reference should be passed if enableProjectingBedsFromRef2Asm is true. (Default: '')"
        enableRunningSecphase : "If true it will run secphase in the marker mode using the wdl parameters starting with 'secphase' otherwise skip it. (Default: false)"
        secphaseDockerImage: "Docker image for running Secphase (Default: mobinasri/secphase:v0.4.3)"
        secphaseOptions: "String containing secphase options (can be either --hifi or --ont). (Default --hifi)"
        secphaseVersion: "Secphase version. (Default: v0.4.3)"
        enableOutputtingWig: "If true it will make wig files from cov files and output them. wig files can be easily imported into IGV sessions (Default: true)"
        enableOutputtingBam: "If true it will output read alignment bam file and its related index (Default: false)"
        windowSize: "The size of the window flagger uses for finding coverage distrubutions (Default: 5000000)"
        sortPdfPagesByHaplotype: "Sort the coverage distributions plotted in the output pdf by haplotype (Default: false)"
        hap1ContigPattern: "The pattern that will be used for finding the names of the contigs belonging to haplotype1. It will be skipped if sortPdfPagesByHaplotype is false. (Default: hap1)"
        hap2ContigPattern: "The pattern that will be used for finding the names of the contigs belonging to haplotype2. It will be skipped if sortPdfPagesByHaplotype is false. (Default: hap2)"
    }
    input{
        String sampleName
        String suffixForMapping
        String suffixForFlagger
        File hap1AssemblyFasta
        File hap2AssemblyFasta
        Array[File] readFiles
        String preset 

        String aligner="winnowmap"
        Int kmerSize = 15
        String alignerOptions="--eqx -Y -L -y"
        String readExtractionOptions="-TMM,ML,Mm,Ml"
        File? referenceFastaForReadExtraction
        Boolean enableAddingMDTag=true
        Int splitNumber = 16
        Boolean enableSplittingReadsEqually=false
        
        Int minReadLengthForMapping = 0
        Int alignerThreadCount = 16
        Int alignerMemSize = 48
        String alignerDockerImage = "mobinasri/long_read_aligner:v0.4.0"

        Float maxReadDivergenceForFlagger = 0.1
        String correctBamOptions="--primaryOnly --minReadLen 5000 --minAlignment 5000 --maxDiv 0.1"

        Array[File] potentialBiasesBedArray = []
        File? sexBed
        File? SDBed
        File? cntrBed # censat annotation with no "ct"
        File? cntrCtBed
        Array[File] additionalStratificationBedArray=[]
        Array[String] additionalStratificationNameArray=[]
        
        Boolean enableProjectingBedsFromRef2Asm = true
        File projectionReferenceFasta = ""

        # Secphase will be executed in the alignment task
        # and the input bam file created for flagger task
        # has the secphase log file applied on the alignments
        Boolean enableRunningSecphase = false
        String secphaseDockerImage = "mobinasri/secphase:v0.4.3--c99e0e9f3561192e127b2d913c932c3e68aa21bf"
        String secphaseOptions = "--hifi"
        String secphaseVersion = "v0.4.3"

        Boolean enableOutputtingWig = true
        Int preemptible=2
        String zones="us-west2-a"
        Boolean enableOutputtingBam=false

        Int windowSize=5000000
        Boolean sortPdfPagesByHaplotype=false
        String hap1ContigPattern="hap1"
        String hap2ContigPattern="hap2"
    }

    # Create a diploid assembly 
    # from the given haplotypes
    call misc_t.createDipAsm {
        input:
            hap1AssemblyFasta = hap1AssemblyFasta,
            hap2AssemblyFasta = hap2AssemblyFasta,
            outputName = "${sampleName}.dip.asm"
    }

    call aligner_t.longReadAlignmentScattered as alignmentTask{
        input:
            readFiles = readFiles,
            assemblyFasta = createDipAsm.diploidAssemblyFastaGz,
            aligner = aligner,
            preset = preset,
            kmerSize = kmerSize,
            alignerOptions = alignerOptions,
            readExtractionOptions = readExtractionOptions,
            sampleName = sampleName,
            suffix = suffixForMapping,
            referenceFastaForReadExtraction = referenceFastaForReadExtraction,
            enableAddingMDTag = enableAddingMDTag,
            splitNumber = splitNumber,
            enableSplittingReadsEqually = enableSplittingReadsEqually,

            # minReadLength will be used for filtering
            # short reads before mapping reads
            minReadLength = minReadLengthForMapping,

            alignerThreadCount = alignerThreadCount,
            alignerMemSize = alignerMemSize,
            alignerDockerImage = alignerDockerImage,

            enableRunningSecphase = enableRunningSecphase,
            secphaseOptions = secphaseOptions,
            secphaseDockerImage = secphaseDockerImage,
            secphaseVersion = secphaseVersion,
            correctBamOptions = correctBamOptions,
            preemptible = preemptible,
            zones = zones
    }

    call flagger_t.FlaggerEndToEnd as flaggerTask{
        input: 
            sampleName = sampleName,
            suffix = suffixForFlagger,
            hap1AssemblyFasta = hap1AssemblyFasta,
            hap2AssemblyFasta = hap2AssemblyFasta,
            readAlignmentBam = alignmentTask.bamFile,
            readAlignmentBai = alignmentTask.baiFile,
            maxReadDivergence = maxReadDivergenceForFlagger,
            potentialBiasesBedArray = potentialBiasesBedArray,
            sexBed = sexBed,
            SDBed = SDBed,
            cntrBed = cntrBed,
            cntrCtBed = cntrCtBed,
            additionalStratificationBedArray = additionalStratificationBedArray,
            additionalStratificationNameArray = additionalStratificationNameArray,
            enableProjectingBedsFromRef2Asm = enableProjectingBedsFromRef2Asm,
            projectionReferenceFasta = projectionReferenceFasta,
            
            # running secphase is disabled here since if enableRunningSecphase
            # was true, secphase must have been executed in the alignment task
            # so it would be redundant if we enable it here again
            enableRunningSecphase = false,
            enableOutputtingWig = enableOutputtingWig,

            windowSize = windowSize,
            sortPdfPagesByHaplotype = sortPdfPagesByHaplotype,
            hap1ContigPattern = hap1ContigPattern,
            hap2ContigPattern = hap2ContigPattern
    }

    if (enableOutputtingBam){
        File readAlignmentBamToOutput = alignmentTask.bamFile
        File readAlignmentBaiToOutput = alignmentTask.baiFile
    }

    output {
        # output bam file if enableOutputtingBam = true
        File? readAlignmentBam = readAlignmentBamToOutput
        File? readAlignmentBai = readAlignmentBaiToOutput          
   
        # bed files and factors for biased regions
        Array[File] detectedBiasedRegionBedArray = flaggerTask.detectedBiasedRegionBedArray
        Array[Float] detectedBiasedRegionFactorArray = flaggerTask.detectedBiasedRegionFactorArray

        # get projected bed files if there is any
        File? projectionSexBed = flaggerTask.projectionSexBed
        File? projectionSDBed = flaggerTask.projectionSDBed
        File? projectionCntrBed = flaggerTask.projectionCntrBed
        Array[File]? projectionAdditionalStratificationBedArray = flaggerTask.projectionAdditionalStratificationBedArray

        # flagger preprocess files
        Float modeCoverageFloat = flaggerTask.modeCoverageFloat
        File covGz = flaggerTask.covGz
        File highMapqCovGz = flaggerTask.highMapqCovGz
        File excludedReadIdsText = flaggerTask.excludedReadIdsText

        # wig files if user enabled outputting them
        File? wig = flaggerTask.wig
        File? highMapqWig = flaggerTask.highMapqWig

        # secphase output
        File? secphaseOutputLog = alignmentTask.secphaseOutputLog
        File? secphaseModifiedReadBlocksMarkersBed = alignmentTask.secphaseModifiedReadBlocksMarkersBed
        File? secphaseMarkerBlocksBed = alignmentTask.secphaseMarkerBlocksBed

        # flagger outputs for all alignments
        File finalBed = flaggerTask.finalBed
        File finalBedNoHap = flaggerTask.finalBedNoHap
        File miscFilesTarGz = flaggerTask.miscFilesTarGz
        File pdf = flaggerTask.pdf

        # flagger statistics stratified by long contigs, centeromere, SD and sex
        # all alignments
        File statsTsv = flaggerTask.statsTsv
        File statsPercOnlyTsv = flaggerTask.statsPercOnlyTsv
    }
}

