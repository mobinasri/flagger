version 1.0

import "../tasks/coverage/cov2counts.wdl" as cov2counts_t
import "../tasks/coverage/cov2counts_by_window.wdl" as cov2counts_by_window_t
import "../tasks/mixture_model/fit_model.wdl" as fit_model_t
import "../tasks/mixture_model/fit_model_by_window.wdl" as fit_model_by_window_t
import "../tasks/mixture_model/find_blocks.wdl" as find_blocks_t
import "../tasks/mixture_model/find_blocks_by_window.wdl" as find_blocks_by_window_t
import "../tasks/mixture_model/pdf_generator.wdl" as pdf_generator_t
import "../tasks/other/bedtools.wdl" as bedtools_t
import "../tasks/mixture_model/fit_model_bed.wdl" as fit_model_bed_t
import "../tasks/coverage/bam_coverage.wdl" as bam_coverage_t

workflow runFlagger{
    input {
        # biasedRegionBedArray:
        #     An array of bed files pointing to regions prone to have
        #     systematic coverage bias
        # biasedRegionNameArray:
        #     An array of names associated with the given array of bed 
        #     files (e.g. ["hifi_low_bias", "hifi_high_bias"])
        # biasedRegionFactorArray:
        #     An array of coverage factors to adjust the expected coverage
        #     for the given bed files (e.g. [0.75, 1.25])
        #     Each factor will be used to be multiplied by covFloat and
        #     obtain the expected coverage
        Array[File] biasedRegionBedArray = [] 
        Array[String] biasedRegionNameArray = []
        Array[Float] biasedRegionFactorArray = []
        File coverageGz
        File highMapqCoverageGz
        File fai
        Float covFloat # the coverage with the highest frequency (most of the time same as mean coverage)
        Boolean isDiploid=false # This is only used for pdf generation and separating the pages for each haplotype
        Int windowSize = 5000000 # Size of windows for spliting assembly and calculating coverage dist
        String sampleName
        String suffix = "flagger"
        Int cov2countsDiskSizeGB = 512
    }

    if (length(biasedRegionBedArray) != 0){
        scatter (biasedRegionData in zip(biasedRegionBedArray, zip(biasedRegionNameArray, biasedRegionFactorArray))){
            File biasedRegionBed = biasedRegionData.left
            String biasedRegionName = biasedRegionData.right.left
            Float biasedRegionFactor = biasedRegionData.right.right
            call bedtools_t.merge {
                input:
                    bed = biasedRegionBed,
                    margin = 50000,
                    outputPrefix = basename(biasedRegionBed, ".bed")
            }
            call fit_model_bed_t.runFitModelBed as biasedRegionModels {
                input:
                    bed = merge.mergedBed,
                    suffix = biasedRegionName,
                    coverageGz = coverageGz,
                    covFloat = covFloat * biasedRegionFactor
            }
        }
    }
    call mergeHsatBeds {
        input:
            bedsTarGzArray = biasedRegionModels.bedsTarGz
    }
    call cov2counts_t.cov2counts {
        input:
            coverageGz = coverageGz,
            diskSize = cov2countsDiskSizeGB
    }
    call fit_model_t.fitModel {
        input:
            counts = cov2counts.counts,
            cov = covFloat
    }
    call find_blocks_t.findBlocks {
        input:
            coverageGz = coverageGz,
            table = fitModel.probabilityTable
    }
    call cov2counts_by_window_t.cov2countsByWindow {
        input:
            coverageGz = coverageGz,
            excludeBedArray = biasedRegionBedArray,
            fai = fai,
            windowSize = windowSize,
            diskSize = cov2countsDiskSizeGB
    }
    call fit_model_by_window_t.fitModelByWindow {
        input:
            windowsText = cov2countsByWindow.windowsText,
            countsTarGz = cov2countsByWindow.windowCountsTarGz,
            cov = covFloat 
    }
    call find_blocks_by_window_t.findBlocksByWindow {
        input:
            windowCovsTarGz = cov2countsByWindow.windowCovsTarGz,
            windowProbTablesTarGz = fitModelByWindow.windowProbTablesTarGz,
            windowsText = cov2countsByWindow.windowsText
    }
    call pdf_generator_t.pdfGenerator {
        input:
            windowProbTablesTarGz = fitModelByWindow.windowProbTablesTarGz,
            genomeProbTable = fitModel.probabilityTable,
            isDiploid = isDiploid
    }
    call combineBeds as combineWindowBased{
        input:
            outputPrefix = "window_corrected",
            firstPrefix = "whole_genome",
            secondPrefix = "window_based",
            firstBedsTarGz = findBlocks.bedsTarGz,
            secondBedsTarGz = findBlocksByWindow.windowBedsTarGz
    }
    call combineBeds as combineHsatBased{
       input:
            outputPrefix = "bias_corrected",
            firstPrefix = "window_corrected",
            secondPrefix = "hsat_based",
            firstBedsTarGz = combineWindowBased.combinedBedsTarGz,
            secondBedsTarGz = mergeHsatBeds.bedsTarGz
    }    
    call dupCorrectBeds {
        input:
            covGz = coverageGz,
            highMapqCovGz = highMapqCoverageGz,
            bedsTarGz = combineHsatBased.combinedBedsTarGz,
            prefix="bias_corrected"
    }
    call filterBeds {
        input:
            fai = fai,
            dupCorrectedBedsTarGz = dupCorrectBeds.dupCorrectedBedsTarGz
    }

    call getFinalBed {
        input:
            bedsTarGz = filterBeds.filteredBedsTarGz,
            sampleName = sampleName,
            suffix = suffix
    }

    call gatherFiles {
        input:
            files = [cov2counts.counts, fitModel.probabilityTable, findBlocks.bedsTarGz, cov2countsByWindow.windowCountsTarGz, cov2countsByWindow.windowCovsTarGz, fitModelByWindow.windowProbTablesTarGz, findBlocksByWindow.windowBedsTarGz, combineWindowBased.combinedBedsTarGz, dupCorrectBeds.dupCorrectedBedsTarGz, filterBeds.filteredBedsTarGz, combineHsatBased.combinedBedsTarGz],
            outputName = "${sampleName}.${suffix}.miscellaneous"
    }

    output {
        File miscFilesTarGz = gatherFiles.outputTarGz
        File pdf = pdfGenerator.pdf
        File finalBed = getFinalBed.finalBed 
    }
}

task gatherFiles {
    input {
        Array[File] files
        String outputName
        # runtime configurations
        Int memSize=8
        Int threadCount=4
        Int diskSize=128
        String dockerImage="mobinasri/flagger:v0.3.2"
        Int preemptible=2
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        mkdir ~{outputName}
        cp ~{sep=" " files} ~{outputName}
        tar -cf ~{outputName}.tar ~{outputName}
        pigz -p~{threadCount} ~{outputName}.tar
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File outputTarGz = glob("*.tar.gz")[0]
    }
}

task String2Float {
    input {
        String str
    }
    command <<<
        echo ~{str} > file.txt
    >>>
    runtime {
        docker: "mobinasri/flagger:v0.3.2"
        memory: "1 GB"
        cpu: 1
        disks: "local-disk 1 SSD"
    }
    output {
        Float number = read_float("file.txt")
    }
}
task combineBeds {
    input {
        File firstBedsTarGz
        File secondBedsTarGz
        String firstPrefix
        String secondPrefix
        String outputPrefix = "combined"
        # runtime configurations
        Int memSize=8
        Int threadCount=4
        Int diskSize=128
        String dockerImage="mobinasri/flagger:v0.3.2"
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
        
        mkdir first second
        tar --strip-components 1 -xvzf ~{firstBedsTarGz} --directory first
        tar --strip-components 1 -xvzf ~{secondBedsTarGz} --directory second
                
        FILENAME=~{firstBedsTarGz}
        PREFIX=$(basename ${FILENAME%.*.*.tar.gz})
        
        cat second/*.bed | sort -k1,1 -k2,2n | bedtools merge -i - > second_all.bed 
        mkdir first_minus_second ~{outputPrefix}
        for c in error duplicated haploid collapsed
        do
            bedtools subtract -sorted -a first/${PREFIX}.~{firstPrefix}.${c}.bed -b second_all.bed > first_minus_second/${PREFIX}.${c}.bed
            cat first_minus_second/*.${c}.bed second/*.${c}.bed | sort -k1,1 -k2,2n | bedtools merge -i - > ~{outputPrefix}/${PREFIX}.~{outputPrefix}.${c}.bed
        done

        tar -cf ${PREFIX}.beds.~{outputPrefix}.tar ~{outputPrefix}
        gzip ${PREFIX}.beds.~{outputPrefix}.tar

    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File combinedBedsTarGz = glob("*.beds.${outputPrefix}.tar.gz")[0]
    }
}


task dupCorrectBeds {
    input {
        File covGz
        File highMapqCovGz
        File bedsTarGz
        String prefix
        Int minCov=4
        # runtime configurations
        Int memSize=16
        Int threadCount=8
        Int diskSize=128
        String dockerImage="mobinasri/flagger:v0.3.2"
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
 
        FILENAME=$(basename ~{highMapqCovGz})
        PREFIX=${FILENAME%.cov.gz}

        mkdir ~{prefix}
        tar --strip-components 1 -xvzf ~{bedsTarGz} --directory ~{prefix}

        zcat ~{highMapqCovGz} | \
            awk '{if(substr($1,1,1) == ">") {contig=substr($1,2,40)} else if($3 >= ~{minCov}) {print contig"\t"$1-1"\t"$2}}' | \
            sort -k1,1 -k2,2n | \
            bedtools merge -i - > high_mapq.bed


        zcat ~{covGz} | \
            awk '{if(substr($1,1,1) == ">") {contig=substr($1,2,40)} else if($3 < ~{minCov}) {print contig"\t"$1-1"\t"$2}}' | \
            sort -k1,1 -k2,2n | \
            bedtools merge -i - > extremely_low.bed

        mkdir dup_corrected

        cat high_mapq.bed extremely_low.bed | sort -k1,1 -k2,2n | bedtools merge -i - > exclude_dup.bed

        # do the correction
        bedtools subtract -sorted -a ~{prefix}/${PREFIX}.~{prefix}.duplicated.bed -b exclude_dup.bed > dup_corrected/${PREFIX}.dup_corrected.duplicated.bed
        bedtools intersect -sorted -a ~{prefix}/${PREFIX}.~{prefix}.duplicated.bed -b high_mapq.bed > dup_to_hap.bed
        ##bedtools intersect -sorted -a ~{prefix}/${PREFIX}.~{prefix}.duplicated.bed -b extremely_low.bed > dup_to_err.bed
        cat dup_to_hap.bed ~{prefix}/${PREFIX}.~{prefix}.haploid.bed | sort -k1,1 -k2,2n | bedtools merge -i - > dup_corrected/${PREFIX}.dup_corrected.haploid.bed
        cat extremely_low.bed ~{prefix}/${PREFIX}.~{prefix}.error.bed | sort -k1,1 -k2,2n | bedtools merge -i - > dup_corrected/${PREFIX}.dup_corrected.error.bed
        
        # just copy collapsed comp
        cp ~{prefix}/${PREFIX}.~{prefix}.collapsed.bed dup_corrected/${PREFIX}.dup_corrected.collapsed.bed

        tar -cf ${PREFIX}.beds.dup_corrected.tar dup_corrected
        gzip ${PREFIX}.beds.dup_corrected.tar

    >>>

    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File dupCorrectedBedsTarGz = glob("*.beds.dup_corrected.tar.gz")[0]
    }
}

task filterBeds {
    input {
        File fai
        File dupCorrectedBedsTarGz
        Int mergeLength=100
        Int minBlockLength=1000
        # runtime configurations
        Int memSize=8
        Int threadCount=4
        Int diskSize=32
        String dockerImage="mobinasri/flagger:v0.3.2"
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
        
        mkdir dup_corrected
        tar --strip-components 1 -xvzf ~{dupCorrectedBedsTarGz} --directory dup_corrected

        FILENAME=~{dupCorrectedBedsTarGz}
        PREFIX=$(basename ${FILENAME%.*.*.tar.gz})

        mkdir initial filtered
        for c in error duplicated haploid collapsed
        do
            bedtools merge -d ~{mergeLength} -i dup_corrected/*.${c}.bed | awk '($3-$2) >= ~{minBlockLength}' > initial/${PREFIX}.filtered.${c}.bed
        done

        
        # Gather ambiguous overlaps
        bedtools intersect -sorted -a initial/${PREFIX}.filtered.error.bed -b initial/${PREFIX}.filtered.haploid.bed > err_hap.overlap.bed
        bedtools intersect -sorted -a initial/${PREFIX}.filtered.duplicated.bed -b initial/${PREFIX}.filtered.haploid.bed > dup_hap.overlap.bed
        bedtools intersect -sorted -a initial/${PREFIX}.filtered.collapsed.bed -b initial/${PREFIX}.filtered.haploid.bed > col_hap.overlap.bed
        bedtools intersect -sorted -a initial/${PREFIX}.filtered.error.bed -b initial/${PREFIX}.filtered.duplicated.bed > err_dup.overlap.bed

        # Gather small blocks and assign as unknown
        cat initial/* | sort -k1,1 -k2,2n | bedtools merge -i - > initial/all.bed
        cat ~{fai} | awk '{print $1"\t"0"\t"$2}' | sort -k1,1 -k2,2n > asm.bed
        bedtools subtract -sorted -a asm.bed -b initial/all.bed > initial/${PREFIX}.filtered.unknown.bed
        # Add ambiguous overlaps as unknown
        cat *.overlap.bed initial/${PREFIX}.filtered.unknown.bed | sort -k1,1 -k2,2n | bedtools merge -i - > filtered/${PREFIX}.filtered.unknown.bed         

        # Subtract unknown regions
        for c in error duplicated haploid collapsed
        do
            bedtools subtract -sorted -a initial/${PREFIX}.filtered.${c}.bed -b filtered/${PREFIX}.filtered.unknown.bed > filtered/${PREFIX}.filtered.${c}.bed
        done

        tar -cf ${PREFIX}.beds.filtered.tar filtered
        gzip ${PREFIX}.beds.filtered.tar
        
    >>>

    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File filteredBedsTarGz = glob("*.beds.filtered.tar.gz")[0]
    }
}


task mergeHsatBeds {
    input {
        Array[File]? bedsTarGzArray
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=32
        String dockerImage="mobinasri/flagger:v0.3.2"
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

        if [[ -n "~{sep="" bedsTarGzArray}" ]]; then 
            FILENAMES=(~{sep=" " bedsTarGzArray})
            FILENAME=${FILENAMES[0]}
            PREFIX=$(basename ${FILENAME%.*.*.tar.gz})

            mkdir hsat_unmerged hsat_based
            for s in ~{sep=" " bedsTarGzArray}; do
                tar --strip-components 1 -xvzf $s --directory hsat_unmerged
            done
 
            for comp in error haploid duplicated collapsed; do
                cat hsat_unmerged/*.${comp}.bed | sort -k1,1 -k2,2n | bedtools merge -i - > hsat_based/$PREFIX.hsat_based.${comp}.bed
            done
        else
            mkdir hsat_based
            PREFIX="empty"
            for comp in error haploid duplicated collapsed; do
                touch hsat_based/${PREFIX}.hsat_based.${comp}.bed
            done
        fi
        tar -cf ${PREFIX}.beds.hsat_based.tar hsat_based
        gzip ${PREFIX}.beds.hsat_based.tar
    >>>

    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File bedsTarGz = glob("*.beds.hsat_based.tar.gz")[0]
    }
}

task getFinalBed {
    input {
        File bedsTarGz
        String sampleName
        String suffix
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=32
        String dockerImage="mobinasri/flagger:v0.3.2"
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
        bash /home/scripts/combine_comp_beds.sh \
            -b ~{bedsTarGz} \
            -m /home/scripts/colors.txt \
            -t ~{sampleName}.~{suffix} \
            -o output/~{sampleName}.~{suffix}.flagger_final.bed
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File finalBed = glob("output/*.flagger_final.bed")[0]
    }
}
