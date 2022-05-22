version 1.0

import "../tasks/cov2counts.wdl" as cov2counts_t
import "../tasks/bam_coverage.wdl" as bamCoverage_t
import "../tasks/bam2paf.wdl" as bam2paf_t
import "../tasks/subset_coverage.wdl" as subsetCoverage_t

workflow runChimericAnalysis{
    input {
        File bamFile
        String sampleSuffix
        String platform
        String sampleName
        File assemblyFastaGz 
    }
    call bamCoverage_t.bamCoverage {
        input:
            sampleName = sampleName,
            sampleSuffix = sampleSuffix, 
            platform = platform,
            bamFiles = [bamFile],
            minMAPQ = 0,
            extraOptions = "-G 2048",
            assemblyFastaGz = assemblyFastaGz
    }
    call bam2paf_t.bam2paf {
        input:
            bamFile = bamFile,
            minMAPQ = 0,
            primaryOnly = "yes",
            extraOptions = "-f2048"
    }
    call getFlankingBed {
        input:
            paf = bam2paf.pafFile,
            margin = 1000
    }
    call subsetCoverage_t.subsetCoverage {
        input:
            coverageGz = bamCoverage.coverageGz,
            blocksBed = getFlankingBed.bed,
            suffix = "supp_ends"
    }
    call cov2counts_t.cov2counts{
         input:
            coverageGz = subsetCoverage.outputCoverageGz
    }
    output {
        File suppFlankingBeds = getFlankingBed.bed
        File suppFlankingCovGz = subsetCoverage.outputCoverageGz
        File suppFlankingCounts = cov2counts.counts
    }
}

task getFlankingBed {
    input {
        File paf
        Int margin=1000
        # runtime configurations
        Int memSize=16
        Int threadCount=8
        Int diskSize=128
        String dockerImage="quay.io/masri2019/hpp_coverage:latest"
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
        
        FILENAME=$(basename ~{paf})
        PREFIX=${FILENAME%.paf}
        
        mkdir output
        cat ~{paf} | awk -v m=~{margin} '{s1 = $8 - m; e1 = $8 + m; \
                                          s2 = $9 - m; e2 = $9 + m; \
                                          if(s1 < 0){s1 = 0}; if(s2 < 0){s2 = 0}; \
                                          if(e1 > $7){e1 = $7}; if(e2 > $7){e2 = $7}; \
                                          print $6"\t"s1"\t"e1; print $6"\t"s2"\t"e2}' | \
                      bedtools sort -i - | bedtools merge -i - > output/${PREFIX}.bed

    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File bed = glob("output/*.bed")[0]
    }
}


