version 1.0

import "../tasks/coverage/make_summary_table.wdl" as make_summary_table_t
import "../../ext/secphase/wdls/workflows/correct_bam.wdl" as correct_bam_t
import "../tasks/coverage/bam_coverage.wdl" as cov_t
import "../tasks/other/misc.wdl" as misc_t
import "../tasks/hmm_flagger/hmm_flagger.wdl" as hmm_flagger_t
import "../tasks/other/project_blocks_for_flagger.wdl" as project_t
import "../../ext/secphase/wdls/workflows/secphase.wdl" as secphase_t
import "../tasks/alignment/asm2asm_aligner.wdl" as asm2asm_t 
import "../tasks/alignment/produce_fai.wdl" as fai_t
import "../tasks/coverage/cov2wig.wdl" as cov2wig_t
import "../tasks/coverage/augment_coverage_by_labels.wdl" as augment_cov_t

workflow HMMFlaggerEndToEnd{
    meta {
        author: "Mobin Asri"
        email: "masri@ucsc.edu"
        description: "Running HMM-Flagger for evaluating a diploid assembly using long read alignments. More information at https://github.com/mobinasri/flagger"
    }
    parameter_meta {
        sampleName: "(Required) Sample name for example 'HG002'"
        suffix: "(Required) Suffix string that contains information about this analysis for example 'hifi_winnowmap_flagger_for_hprc'"
        hap1AssemblyFasta: "(Required) Path to uncompressed or gzip-compressed fasta file of the 1st haplotype."
        hap2AssemblyFasta: "(Required) Path to uncompressed or gzip-compressed fasta file of the 2nd haplotype."
        readAlignmentBam: "(Required) Path to sorted read alignment bam."
        readAlignmentBai: "(Required) Path to bam index for read alignment"
        alphaTsv : "(Required) The dependency factors for adjusting emission parameters with previous emission. This parameter is a tsv file with 4 rows and 4 columns with no header line. All numbers should be between 0 and 1."
        maxReadDivergence: "Alignments with gap-compressed ratio higher than this will be filtered in the pre-process step. (Default: 0.1)"
        minReadLength: "Minimum read size (Default: 5000)"
        minAlignmentLength: "Minimum alignment size (Default: 5000)"
        downSamplingRate: "Rate of downsampling (Default: 1.0 which means no down-sampling)"
        sexBed: "(Optional) bed file containing regions assigned to X/Y chromosomes. (in asm coordinates)"
        sexBedToBeProjected: "(Optional) bed file containing regions assigned to X/Y chromosomes. (in ref coordinates)"
        SDBed: "(Optional) Bed file containing Segmental Duplications. (in asm coordinates)"
        SDBedToBeProjected: "(Optional) Bed file containing Segmental Duplications. (in ref coordinates)"
        cntrBed: "(Optional) Bed file containing peri/centromeric satellites (ASat, HSat, bSat, gSat) without 'ct' blocks. (in asm coordinates)"
        cntrBedToBeProjected: "(Optional) Bed file containing peri/centromeric satellites (ASat, HSat, bSat, gSat) without 'ct' blocks. (in ref coordinates)"
        cntrCtBed: "(Optional) Bed file containing centromere transition 'ct' blocks. (in asm coordinates)"
        cntrCtBedToBeProjected: "(Optional) Bed file containing centromere transition 'ct' blocks. (in ref coordinates)"
        annotationsBedArray: "(Optional) list of annotations to be used for augmenting coverage files and stratifying final HMM-Flagger results. These annotations should be in the coordinates of the assembly under evaluation."
        annotationsBedArrayToBeProjected: "(Optional) list of annotations to be used for augmenting coverage files and stratifying final HMM-Flagger results. These annotations are not in the coordinates of the assembly and they need to be projected from the given projection reference to the assembly coordinates."
        biasAnnotationsBedArray: "(Optional) list similar to annotationsBedArray but these annotations potentially have read coverage biases like HSat2. (in the assembly coordinates)"
        biasAnnotationsBedArrayToBeProjected: "(Optional) list similar to annotationsBedArrayToBeProjected but these annotations potentially have read coverage biases like HSat2. (in the reference coordinates)"
        projectionReferenceFasta: "(Optional) If any of the parameters ending with 'ToBeProjected' is not empty a reference fasta should be passed for performing the necessary projections (Default: '')"
        enableRunningSecphase : "If true it will run secphase in the marker mode using the parameters starting with 'secphase' otherwise skip it. (Default: false)"
        secphaseDockerImage: "Docker image for running Secphase (Default: mobinasri/secphase:v0.4.3)"
        secphaseOptions: "String containing secphase options (can be either --hifi or --ont). (Default --hifi)"
        secphaseVersion: "Secphase version. (Default: v0.4.3)"
        enableOutputtingBigWig: "If true it will make bigwig files from cov files and output them. bigwig files can be easily imported into IGV sessions (Default: true)"
        includeContigListText : "(Optional) Create coverage file and run HMM-Flagger only on these contigs (listed in a text file with one contig name per line). (Default: all contigs)"
        binArrayTsv : "(Optional)  A tsv file (tab-delimited) that contains bin arrays for stratifying results by event size. Bin intervals can have overlap. It should contain three columns. 1st column is the closed start of the bin and the 2nd column is the open end. The 3rd column has a name for each bin. (Default: all sizes in a single bin named ALL_SIZES)"
        chunkLen : "The length of chunks for running HMM-Flagger. Each chunk will be processed in a separate thread before merging results together. (Default: 20000000)"
        windowLen : "The length of windows for running HMM-Flagger. The coverage values will be averaged over the bases in each window and then the average value will be considered as an emission. (Default: 4000)"
        labelNames : "The names of the labels/states for reporting in the final summary tsv files (Default: 'Err,Dup,Hap,Col')"
        trackName : "The track name in the final BED file (Default: hmm_flagger_v1.0)"
        numberOfIterations : "Number of EM iterations for estimating HMM parameters (Default:100)"
        convergenceTolerance : "Convergence tolerance. The EM iteration will stop once the difference between all model parameter values in two consecutive iterations is less than this value. (Default = 0.001)"
        maxHighMapqRatio : "Maximum ratio of high mapq coverage for duplicated state (Default = 0.25)"
        minHighMapqRatio : "Minimum ratio of high mapq coverage for collapsed state (Default = 0.5)"
        flaggerMoreOptions : "(Optional) More options for HMM-Flagger provided in a single string (Default = '')"
        modelType : "Model type can be either 'gaussian', 'negative_binomial', or 'trunc_exp_gaussian' (Default = 'trunc_exp_gaussian')"
        flaggerMinimumBlockLenArray : "Array of minimum lengths for converting short non-Hap blocks into Hap blocks. Given numbers should be related to the states Err, Dup and Col respectively. (Default: [0,0,0])"
        flaggerMemSize : "Memory size in GB for running HMM-Flagger (Default : 32)"
        flaggerThreadCount : "Number of threads for running HMM-Flagger (Default : 16)"
        flaggerDockerImage : "Docker image for HMM-Flagger (Default : mobinasri/flagger:v1.1.0)"
        truthBedForMisassemblies : "(Optional) A BED file containing the coordinates and labels of the truth misassemblies. It can be useful when the misassemblies are simulated (e.g. with Falsifier) (Default: None)"
    }
    input{
        String sampleName
        String suffix

        File hap1AssemblyFasta
        File hap2AssemblyFasta
        File readAlignmentBam
	File readAlignmentBai
	File alphaTsv
        Float maxReadDivergence = 0.1
        Int minReadLength = 5000
        Int minAlignmentLength = 5000
        Float downSamplingRate = 1.0
              
        File? includeContigListText
        File? binArrayTsv
        Int chunkLen = 20000000
        Int windowLen = 4000
        String labelNames = "Err,Dup,Hap,Col"
        String trackName = "hmm_flagger_v1.1.0"
        Int numberOfIterations = 100
        Float convergenceTolerance = 0.001
        Float maxHighMapqRatio=0.25
        Float minHighMapqRatio=0.5
        String? flaggerMoreOptions
        String modelType = "trunc_exp_gaussian"
        Array[Int] flaggerMinimumBlockLenArray = []
        Int flaggerMemSize=32
        Int flaggerThreadCount=16
        String flaggerDockerImage="mobinasri/flagger:v1.1.0"

        File? sexBed
        File? SDBed
        File? cntrBed # censat annotation with no "ct"
        File? cntrCtBed

        File? sexBedToBeProjected
        File? SDBedToBeProjected
        File? cntrBedToBeProjected # censat annotation with no "ct"
        File? cntrCtBedToBeProjected
        
        Array[File] annotationsBedArray = []
        Array[File] annotationsBedArrayToBeProjected = []
        Array[File] biasAnnotationsBedArray = []
        Array[File] biasAnnotationsBedArrayToBeProjected = []

        File? projectionReferenceFasta

        Boolean enableRunningSecphase = false
        String secphaseDockerImage = "mobinasri/secphase:v0.4.3--c99e0e9f3561192e127b2d913c932c3e68aa21bf"
        String secphaseOptions = "--hifi"
        String secphaseVersion = "v0.4.3"

        Boolean enableOutputtingBigWig = true
        File? truthBedForMisassemblies

    }

    # Create a diploid assembly 
    # from the given haplotypes
    call misc_t.createDipAsm {
        input:
            hap1AssemblyFasta = hap1AssemblyFasta,
            hap2AssemblyFasta = hap2AssemblyFasta,
            outputName = "${sampleName}.dip.asm"
    }

    # Index diploid assembly
    call fai_t.produceFai {
        input:
            fasta = createDipAsm.diploidAssemblyFastaGz
    }

    # Run Secphase if it is enabled by user
    # Secphase is for fixing reads that were
    # not mapped to the correct haplotype
    # more info https://github.com/mobinasri/secphase
    if (enableRunningSecphase) {
        call secphase_t.runSecPhase as secphase{
            input:
                inputBam = readAlignmentBam,
                diploidAssemblyFastaGz = createDipAsm.diploidAssemblyFastaGz,
                secphaseOptions = secphaseOptions,
                secphaseDockerImage = secphaseDockerImage,
                version = secphaseVersion,
        }
    }

    # Preprocess read alignment file by applying
    # the changes suggested by secphase and
    # filtering reads/alignments based on length
    # and divergence rate
    if (enableRunningSecphase || (0 < minReadLength) || (0 < minAlignmentLength) || (maxReadDivergence < 1.0)){
        call correct_bam_t.correctBam {
            input:
                bam = readAlignmentBam,
                phasingLogText = secphase.outLog,
                suffix = "corrected",
                options = "--primaryOnly --minReadLen ${minReadLength} --minAlignment ${minAlignmentLength} --maxDiv ${maxReadDivergence}",
                flagRemoveSupplementary = false,
                flagRemoveMultiplePrimary = false,
                dockerImage = secphaseDockerImage,
        }
    }


    Boolean needsAtLeastOneProjection = (0 < length(annotationsBedArrayToBeProjected)) || (0 < length(biasAnnotationsBedArrayToBeProjected) || defined(sexBedToBeProjected) || defined(SDBedToBeProjected) || defined(cntrBedToBeProjected))

    # Map each haplotype to the given reference
    # and then project bed files to the assembly 
    # coordinates. This projection will be performed
    # only for the annotations provided in the parameters
    # ending with ToBeProjected
    #
    #   1. annotationsBedArrayToBeProjected
    #   2. biasAnnotationsBedArrayToBeProjected
    #   3. sexBedToBeProjected
    #   4. SDBedToBeProjected
    #   5. cntrBedToBeProjected
    #   6. cntrCtBedToBeProjected (can be ignored since it is only for better projection of cntr blocks)
    #
    # if both arrays are empty then no projection will
    # be executed
    if (needsAtLeastOneProjection && defined(projectionReferenceFasta)) {
        call asm2asm_t.asm2asmAlignment as hap1ToRef {
            input :
                aligner="minimap2",
                preset="asm5",
                queryAssemblyFasta=hap1AssemblyFasta,
                refAssemblyFasta=select_first([projectionReferenceFasta]),
        }
        call asm2asm_t.asm2asmAlignment as hap2ToRef {
            input :
                aligner="minimap2",
                preset="asm5",
                queryAssemblyFasta=hap2AssemblyFasta,
                refAssemblyFasta=select_first([projectionReferenceFasta]),
        }
        call project_t.runProjectBlocksForFlagger as project{
            input:
                sampleName = sampleName,
                hap1AssemblyBam = hap1ToRef.sortedBamFile,
                hap2AssemblyBam = hap2ToRef.sortedBamFile, 
                refBiasedBlocksBedArray = biasAnnotationsBedArrayToBeProjected,
                additionalBedArray = annotationsBedArrayToBeProjected,
                refSexBed = sexBedToBeProjected,
                refSDBed = SDBedToBeProjected,
                refCntrBed = cntrBedToBeProjected,
                refCntrCtBed = cntrCtBedToBeProjected,
        }
    }

    # Get bed files containing potentially biased 
    # regions in asm coordinates, which can be 
    # projections from reference
    Array[File] potentialBiasesBedArrayInAsmCoor = flatten([select_first([project.projectionBiasedBedArray, []]), select_first([biasAnnotationsBedArray, []])])
    
    # Get bed files containing additional stratifications
    # regions in asm coordinates, which can be
    # projections from reference
    Array[File] additionalStratificationBedArrayInAsmCoor = flatten([select_first([project.projectionAdditionalBedArray, []]), select_first([annotationsBedArray, []])])
     
    # Get bed files for SD, centromere and 
    # sex which can be projections from reference

    if (defined(project.projectionSDBed) || defined(SDBed)){
        File SDBedInAsmCoor = select_first([project.projectionSDBed, SDBed])
    }

    if (defined(project.projectionCntrBed) || defined(cntrBed)){
        File cntrBedInAsmCoor = select_first([project.projectionCntrBed, cntrBed])
    }

    if (defined(project.projectionSexBed) || defined(sexBed)){
        File sexBedInAsmCoor = select_first([project.projectionSexBed, sexBed])
    }

    call collectAnnotations{
        input:
            fasta = createDipAsm.diploidAssemblyFastaGz,
            minContigLen = 1000000,
            minContigLenName = "1Mb",
            includeContigListText = includeContigListText,
            biasAnnotationsBedArray = potentialBiasesBedArrayInAsmCoor,
            otherAnnotationsBedArray = additionalStratificationBedArrayInAsmCoor,
            difficultBed_1 = cntrBedInAsmCoor,
            difficultString_1 = "Cntr",
            difficultBed_2 = SDBedInAsmCoor,
            difficultString_2 = "SD",
            sexBed = sexBedInAsmCoor,
            dockerImage = flaggerDockerImage
    }

    # convert bam to cov and add annotation indices to the cov file
    call cov_t.bam2cov {
        input:
            bam = select_first([correctBam.correctedBam, readAlignmentBam]),
            bai = select_first([correctBam.correctedBamIndex, readAlignmentBai]),
            fasta = createDipAsm.diploidAssemblyFastaGz,
            suffix = "",
            mapqThreshold = 20,
            clipRatioThreshold = 0.1,
            downsampleRate = downSamplingRate,
            annotationBedArray = collectAnnotations.annotationBedArray,
            biasAnnotationNameArray = collectAnnotations.biasAnnotationNameArray,
            baselineAnnotationName = "WHOLE_GENOME_DEFAULT",
            includeContigListText = includeContigListText,
            runBiasDetection = (length(collectAnnotations.biasAnnotationNameArray) > 0),
            format = "all",
            memSize = 32,
            threadCount = 8,
            dockerImage = flaggerDockerImage,
    }


    call hmm_flagger_t.hmmFlagger {
        input:
            coverage = bam2cov.coverageGz,
            binArrayTsv = binArrayTsv,
            chunkLen = chunkLen,
            windowLen = windowLen,
            labelNames = labelNames,
            trackName = trackName,
            numberOfIterations = numberOfIterations,
            convergenceTolerance = convergenceTolerance,
            maxHighMapqRatio = maxHighMapqRatio,
            minHighMapqRatio = minHighMapqRatio,
            moreOptions = flaggerMoreOptions,
            minimumBlockLenArray = flaggerMinimumBlockLenArray, 
            alphaTsv = alphaTsv,
            modelType = modelType,
            memSize = flaggerMemSize,
            threadCount = flaggerThreadCount,
            dockerImage = flaggerDockerImage,
    }

    # Get coordinates of canonical bases only (no "N" which may come from scaffolding)
    call misc_t.getCanonicalBasesBed as dipCanonical{
        input:
            assemblyFasta = createDipAsm.diploidAssemblyFastaGz,
            dockerImage = flaggerDockerImage
    }


    # make a truth bed file that contains
    # state indices instead of the names of the states
    if (defined(truthBedForMisassemblies)){
        call misc_t.getIndexLabeledBed as labelTruth{
            input:
                bed = select_first([truthBedForMisassemblies]),
                canonicalBasesDiploidBed = dipCanonical.canonicalBasesBed,
                addHapCoordinates = true
        }
    }

    # make a prediction bed file that contains
    # state indices instead of the names of the states
    call misc_t.getIndexLabeledBed as labelPrediction{
        input:
            bed = hmmFlagger.predictionBed
    }

    # Augment coverage file with prediction labels
    # it will add truth labels if available
    call augment_cov_t.augmentCoverageByLabels {
        input:
            coverage = bam2cov.coverageGz,
            fai = produceFai.fai, 
	    numberOfLabels = 4,
            truthBed = labelTruth.labeledBed,
            predictionBed = labelPrediction.labeledBed,
            includeContigListText = includeContigListText,
            suffix="augmented",
            dockerImage = flaggerDockerImage,
    }

    call make_summary_table_t.makeSummaryTable {
        input:
           coverage = augmentCoverageByLabels.augmentedCoverageGz,
           binArrayTsv = binArrayTsv,
           dockerImage = flaggerDockerImage
    }

    # make bigwig files from cov files
    # bigwig files can be easily imported into IGV sessions
    if (enableOutputtingBigWig) {
        call cov2wig_t.cov2bigwig as cov2bigwig{
            input:
                coverage = bam2cov.coverageGz,
                windowLen = 1000,
                trackName = "${sampleName}.${suffix}",
                fai = produceFai.fai,
                dockerImage = flaggerDockerImage,
        }
    }

    # add labels for N bases and remove Hap labels from bed file
    call getFinalBed {
        input:
            predictionBed = hmmFlagger.predictionBed,
            canonicalBasesDiploidBed = dipCanonical.canonicalBasesBed, 
            sampleName = sampleName,
            suffix = suffix,
            dockerImage = flaggerDockerImage
    }

    if (defined(truthBedForMisassemblies)){
        File benchmarkingSummaryTsvOutput = makeSummaryTable.benchmarkingSummaryTsv
        File contiguitySummaryTsvOutput = makeSummaryTable.contiguitySummaryTsv
    }


    output {
        File coverageGz = augmentCoverageByLabels.augmentedCoverageGz
        File biasTableTsv = bam2cov.biasTableTsv

        File? benchmarkingSummaryTsv = benchmarkingSummaryTsvOutput
        File? contiguitySummaryTsv = contiguitySummaryTsvOutput
        File fullStatsTsv = makeSummaryTable.fullStatsTsv

        File finalBed = getFinalBed.finalBed
        File predictionBed = hmmFlagger.predictionBed
        File loglikelihoodTsv = hmmFlagger.loglikelihoodTsv
        File miscFilesTarGz = hmmFlagger.outputTarGz

        # get projected bed files if there is any
        File? projectionSexBed = project.projectionSexBed
        File? projectionSDBed = project.projectionSDBed
        File? projectionCntrBed = project.projectionCntrBed
        Array[File]? projectionAnnotationsBedArray = project.projectionAdditionalBedArray
        Array[File]? projectionBiasAnnotationsBedArray = project.projectionBiasedBedArray

        # bigwig files if user enabled outputting them
        Array[File]? bigwigArray = cov2bigwig.bigwigArray
        
        # secphase output
        File? secphaseOutputLog = secphase.outLog
        File? secphaseModifiedReadBlocksMarkersBed = secphase.modifiedReadBlocksMarkersBed
        File? secphaseMarkerBlocksBed = secphase.markerBlocksBed
    }
}

task getFinalBed {
    input {
        File predictionBed
        File canonicalBasesDiploidBed
        String sampleName
        String suffix
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=32
        String dockerImage="mobinasri/flagger:v0.4.0"
        Int preemptible=2
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
        
        # make a bed file for Ns only
        # the color is black with an rgb of "0,0,0"
        # the label is "NNN"
        bedtools subtract \
            -a ~{predictionBed} \
            -b ~{canonicalBasesDiploidBed} | \
            cut -f1-3 | \
            bedtools merge -i - | \
            awk '{print $1"\t"$2"\t"$3"\tNNN\t0\t.\t"$2"\t"$3"\t0,0,0"}' > non_canonical.bed

        # make a BED with no Hap tracks
        cat ~{predictionBed} | \
            grep -v "Hap" | \
            grep -v "^track" | \
            bedtools intersect -a - -b ~{canonicalBasesDiploidBed}  > canonical.no_Hap.bed
        
        
        # add track name
        echo "track name=\"~{sampleName}.~{suffix}\" visibility=2 itemRgb=\"On\"" > output/~{sampleName}.~{suffix}.hmm_flagger.no_Hap.bed

        # merge canonical and non-canonical tracks in the final bed
        cat non_canonical.bed canonical.no_Hap.bed | bedtools sort -i - >> output/~{sampleName}.~{suffix}.hmm_flagger.no_Hap.bed
        
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File finalBed = glob("output/*.no_Hap.bed")[0]
    }
}


task collectAnnotations{
    input{
        File fasta
        Int minContigLen = 1000000
        String minContigLenName = "1Mb"
        File? includeContigListText
        Array[File] biasAnnotationsBedArray = []
        Array[File] otherAnnotationsBedArray = []
        File? difficultBed_1
        File? difficultBed_2
        String? difficultString_1
        String? difficultString_2  
        File? sexBed
        # runtime configurations
        Int memSize=8
        Int threadCount=4
        Int diskSize=32
        String dockerImage="mobinasri/flagger:v1.1.0"
        Int preemptible=2
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace
        

        FA_PREFIX=$(echo $(basename ~{fasta}) | sed -e 's/\.fa$//' -e 's/\.fa.gz$//' -e 's/\.fasta$//' -e 's/\.fasta.gz$//')

        if [[ "~{fasta}" == *.gz ]]; then
            gunzip -c ~{fasta} > ${FA_PREFIX}.fa
        else
            cp ~{fasta} ${FA_PREFIX}.fa
        fi

        samtools faidx ${FA_PREFIX}.fa

        mkdir -p intermediate
        python3 /home/scripts/get_N_coords.py --inputFasta ${FA_PREFIX}.fa > intermediate/N_coords.bed

        if [ -n "~{includeContigListText}" ]
        then
            cat ${FA_PREFIX}.fa.fai | \
                awk '{print $1"\t0\t"$2}' | \
                bedtools sort -i - | \
                grep -F -f ~{includeContigListText} | \
                bedtools subtract -a - -b intermediate/N_coords.bed > whole_genome.bed
        else
            cat ${FA_PREFIX}.fa.fai | \
                awk '{print $1"\t0\t"$2}' | \
                bedtools sort -i - | \
                bedtools subtract -a - -b intermediate/N_coords.bed > whole_genome.bed
        fi

        cat whole_genome.bed | awk '($3-$2)>~{minContigLen}{print $0}' > whole_genome.gt_~{minContigLenName}.bed

        if [ -n "~{difficultBed_1}" ]
        then
            bedtools intersect -a ~{difficultBed_1} -b whole_genome.bed > ~{difficultString_1}.bed
            bedtools intersect -a ~{difficultString_1}.bed -b whole_genome.gt_~{minContigLenName}.bed > ~{difficultString_1}.gt_~{minContigLenName}.bed
            bedtools subtract -a whole_genome.bed -b ~{difficultString_1}.bed > Non_~{difficultString_1}.bed
        fi

        if [ -n "~{difficultBed_2}" ]
        then
            bedtools intersect -a ~{difficultBed_2} -b whole_genome.bed > ~{difficultString_2}.bed
            bedtools intersect -a ~{difficultString_2}.bed -b whole_genome.gt_~{minContigLenName}.bed > ~{difficultString_2}.gt_~{minContigLenName}.bed
            bedtools subtract -a whole_genome.bed -b ~{difficultString_2}.bed > Non_~{difficultString_2}.bed
        fi

        if [[ -n "~{difficultBed_2}" && -n "~{difficultBed_1}" ]]
        then
            cat ~{difficultString_1}.bed ~{difficultString_2}.bed | bedtools sort -i - | bedtools merge -i - > ~{difficultString_1}_or_~{difficultString_2}.bed   
            bedtools intersect -a Non_~{difficultString_1}.bed -b Non_~{difficultString_2}.bed > Non_~{difficultString_1}_and_Non_~{difficultString_2}.bed
            bedtools intersect -a Non_~{difficultString_1}_and_Non_~{difficultString_2}.bed -b whole_genome.gt_~{minContigLenName}.bed > Non_~{difficultString_1}_and_Non_~{difficultString_2}.gt_~{minContigLenName}.bed
        fi

        if [ -n "~{sexBed}" ]
        then
            bedtools intersect -a ~{sexBed} -b whole_genome.bed > sex.bed
            bedtools subtract -a whole_genome.bed -b sex.bed > autosome.bed
        fi

        if [[ -n "~{sexBed}" && -n "~{difficultBed_1}" ]]
        then
            bedtools intersect -a autosome.bed -b Non_~{difficultString_1}.bed > autosome_Non_~{difficultString_1}.bed
        fi

        if [[ -n "~{sexBed}" && -n "~{difficultBed_2}" ]]
        then
            bedtools intersect -a autosome.bed -b Non_~{difficultString_2}.bed > autosome_Non_~{difficultString_2}.bed
        fi

        if [[ -n "~{sexBed}" && -n "~{difficultBed_1}" && -n "~{difficultBed_2}" ]]
        then
            bedtools intersect -a autosome.bed -b Non_~{difficultString_1}_and_Non_~{difficultString_2}.bed > autosome_Non_~{difficultString_1}_and_Non_~{difficultString_2}.bed
            bedtools intersect -a autosome_Non_~{difficultString_1}_and_Non_~{difficultString_2}.bed \
                           -b whole_genome.gt_~{minContigLenName}.bed > autosome_Non_~{difficultString_1}_and_Non_~{difficultString_2}.gt_~{minContigLenName}.bed
        fi

        touch biased_annotation_names.txt

        if [ ~{length(biasAnnotationsBedArray)} -gt "0" ]
        then
            for BED_FILE in ~{sep=" " biasAnnotationsBedArray}
            do
                bedtools intersect -a ${BED_FILE} -b whole_genome.bed > $(basename ${BED_FILE})
                echo $(basename ${BED_FILE%%.bed}) >> biased_annotation_names.txt
            done
        fi

         
        if [ ~{length(otherAnnotationsBedArray)} -gt "0" ]
        then
            for BED_FILE in ~{sep=" " otherAnnotationsBedArray}
            do
                bedtools intersect -a ${BED_FILE} -b whole_genome.bed > $(basename ${BED_FILE})
            done
        fi

    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output{
        Array[File] annotationBedArray = glob("*.bed")
        Array[String] biasAnnotationNameArray = read_lines("biased_annotation_names.txt")
    }
}

