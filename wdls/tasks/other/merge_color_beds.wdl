version 1.0 

workflow runMergeColorBeds{
    call mergeColorBeds
    output {
        File mergedBed = mergeColorBeds.mergedBed
    }
}

task mergeColorBeds {
    input {
        File bedsTarGz
        String trackName
        String suffix
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=32
        String dockerImage="mobinasri/bio_base:v0.1"
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
       
        mkdir beds output
 
        tar --strip-components 1 -xvzf ~{bedsTarGz} --directory beds
         
        ERROR_BED=$(ls beds | grep "error")
        DUPLICATED_BED=$(ls beds | grep "duplicated")
        HAPLOID_BED=$(ls beds | grep "haploid")
        COLLAPSED_BED=$(ls beds | grep "collapsed")

        PREFIX=${ERROR_BED%%.error.*}
        SUFFIX=${ERROR_BED##*.error.}
        OUTPUT_NAME=${PREFIX}.${SUFFIX}

        cat beds/${ERROR_BED} | awk '{print $0"\terror"}' >> temp.bed
        cat beds/${DUPLICATED_BED} | awk '{print $0"\tduplicated"}' >> temp.bed
        cat beds/${HAPLOID_BED} | awk '{print $0"\thaploid"}' >> temp.bed
        cat beds/${COLLAPSED_BED} | awk '{print $0"\tcollapsed"}' >> temp.bed
        
        if [ -z ~{suffix} ]; then
            TRACK_NAME=~{trackName}
        else
            TRACK_NAME=~{trackName}.~{suffix}
        fi
        echo "track name=\"${TRACK_NAME}\" visibility=1 itemRgb="On"" > output/${OUTPUT_NAME}
        bedtools sort -i temp.bed | \
            awk 'BEGIN{c["error"]="255,0,0"; c["duplicated"]="210,153,0"; c["haploid"]="0,102,0"; c["collapsed"]="102,0,102"}{print $0"\t0\t+\t"$2"\t"$3"\t"c[$4]}' >> output/${OUTPUT_NAME}
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File mergedBed = glob("output/*.bed")[0]
    }
}

