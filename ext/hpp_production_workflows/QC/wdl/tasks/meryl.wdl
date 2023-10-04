version 1.0

import "extract_reads.wdl" as extractReads_t
import "shard_reads.wdl" as shardReads_t
import "arithmetic.wdl" as arithmetic_t

workflow runMeryl {

    input {
        Array[File] sampleReadsILM
        Array[File] maternalReadsILM
        Array[File] paternalReadsILM
        File? referenceFasta
        Int kmerSize = 21
        Boolean compress = false
        Int merylCountMemSizeGB = 42
        Int merylCountThreadCount = 64
        Int merylUnionSumMemSizeGB = 32
        Int merylUnionSumThreadCount = 32
        Int merylHapmerMemSizeGB = 24
        Int merylHapmerThreadCount = 32
        Int fileExtractionDiskSizeGB = 256
        String dockerImage = "juklucas/hpp_merqury@sha256:ce62531539733eb3f81c8434b70473c18096ed7e2eae605cfcb333ede2437e1c"
    }

    # extract reads
    scatter (readFile in sampleReadsILM) {
        call extractReads_t.extractReads as sampleReadsExtracted {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSize=4,
                threadCount=4,
                diskSizeGB=fileExtractionDiskSizeGB,
                dockerImage=dockerImage
        }
    }
    scatter (readFile in maternalReadsILM) {
        call extractReads_t.extractReads as maternalReadsExtracted {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSize=4,
                threadCount=4,
                diskSizeGB=fileExtractionDiskSizeGB,
                dockerImage=dockerImage
        }
    }
    scatter (readFile in paternalReadsILM) {
        call extractReads_t.extractReads as paternalReadsExtracted {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSize=4,
                threadCount=4,
                diskSizeGB=fileExtractionDiskSizeGB,
                dockerImage=dockerImage
        }
    }

    # get file size of results (for union sum)
    call arithmetic_t.sum as sampleReadSize {
        input:
            integers=sampleReadsExtracted.fileSizeGB
    }
    call arithmetic_t.sum as maternalReadSize {
        input:
            integers=maternalReadsExtracted.fileSizeGB
    }
    call arithmetic_t.sum as paternalReadSize {
        input:
            integers=paternalReadsExtracted.fileSizeGB
    }
    call arithmetic_t.sum as allReadSize {
        input:
            integers=[sampleReadSize.value, maternalReadSize.value, paternalReadSize.value]
    }


    # get max file size (for meryl count sum)
    call arithmetic_t.max as sampleReadSizeMax {
        input:
            integers=sampleReadsExtracted.fileSizeGB
    }
    call arithmetic_t.max as maternalReadSizeMax {
        input:
            integers=maternalReadsExtracted.fileSizeGB
    }
    call arithmetic_t.max as paternalReadSizeMax {
        input:
            integers=paternalReadsExtracted.fileSizeGB
    }

    # do the meryl counting
    scatter (readFile in sampleReadsExtracted.extractedRead) {
        call merylCount as sampleMerylCount {
            input:
                readFile=readFile,
                kmerSize=kmerSize,
                compress=compress,
                threadCount=merylCountThreadCount,
                memSize=merylCountMemSizeGB,
                diskSizeGB=sampleReadSizeMax.value * 4,
                dockerImage=dockerImage
        }
    }
    scatter (readFile in maternalReadsExtracted.extractedRead) {
        call merylCount as maternalMerylCount {
            input:
                readFile=readFile,
                kmerSize=kmerSize,
                compress=compress,
                threadCount=merylCountThreadCount,
                memSize=merylCountMemSizeGB,
                diskSizeGB=maternalReadSizeMax.value * 4,
                dockerImage=dockerImage
        }
    }
    scatter (readFile in paternalReadsExtracted.extractedRead) {
        call merylCount as paternalMerylCount {
            input:
                readFile=readFile,
                kmerSize=kmerSize,
                compress=compress,
                threadCount=merylCountThreadCount,
                memSize=merylCountMemSizeGB,
                diskSizeGB=paternalReadSizeMax.value * 4,
                dockerImage=dockerImage
        }
    }

    # do the meryl merging
    call merylUnionSum as sampleMerylUnionSum {
        input:
            merylCountFiles=sampleMerylCount.merylDb,
            identifier="sample",
            threadCount=merylUnionSumThreadCount,
            memSize=merylUnionSumMemSizeGB,
            diskSizeGB=sampleReadSize.value * 4,
            dockerImage=dockerImage
    }
    call merylUnionSum as maternalMerylUnionSum {
        input:
            merylCountFiles=maternalMerylCount.merylDb,
            identifier="maternal",
            threadCount=merylUnionSumThreadCount,
            memSize=merylUnionSumMemSizeGB,
            diskSizeGB=maternalReadSize.value * 4,
            dockerImage=dockerImage
    }
    call merylUnionSum as paternalMerylUnionSum {
        input:
            merylCountFiles=paternalMerylCount.merylDb,
            identifier="paternal",
            threadCount=merylUnionSumThreadCount,
            memSize=merylUnionSumMemSizeGB,
            diskSizeGB=paternalReadSize.value * 4,
            dockerImage=dockerImage
    }

    call merylHapmer as merylHapmer {
        input:
            sampleMerylDB=sampleMerylUnionSum.merylDb,
            maternalMerylDB=maternalMerylUnionSum.merylDb,
            paternalMerylDB=paternalMerylUnionSum.merylDb,
            threadCount=merylHapmerThreadCount,
            memSize=merylHapmerMemSizeGB,
            diskSizeGB=allReadSize.value * 2,
            dockerImage=dockerImage
    }

	output {
		File sampleMerylDB = sampleMerylUnionSum.merylDb
		File maternalHapmer = merylHapmer.maternalHapmers
		File paternalHapmer = merylHapmer.paternalHapmers
		File hapmerImages = merylHapmer.hapmerImages
	}
}


task merylCount {
    input {
        File readFile
        Int kmerSize=21
        Boolean compress = false
        Int memSize = 42
        Int threadCount = 64
        Int diskSizeGB = 64
        String dockerImage = "juklucas/hpp_merqury:latest"
    }

    String compress_arg = if compress then "compress" else ""

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
        OMP_NUM_THREADS=~{threadCount}

        # generate meryl db for each read
        ID=`basename ~{readFile} | sed 's/.gz$//' | sed 's/.f[aq]\(st[aq]\)*$//'`
        meryl ~{compress_arg} k=~{kmerSize} threads=~{threadCount} memory=$((~{memSize}-10)) count output $ID.meryl ~{readFile}

        # package
        tar cvf $ID.meryl.tar $ID.meryl

        # cleanup
        rm -rf $ID.meryl
	>>>
	output {
		File merylDb = glob("*.meryl.tar")[0]
	}
    runtime {
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}


task merylUnionSum {
    input {
        Array[File] merylCountFiles
        String identifier
        Int memSize = 32
        Int threadCount = 32
        Int diskSizeGB = 64
        String dockerImage = "juklucas/hpp_merqury:latest"
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
        OMP_NUM_THREADS=~{threadCount}

        # extract meryl dbs
        mkdir extracted/
        cd extracted/
        for m in ~{sep=" " merylCountFiles} ; do
            tar xf $m &
        done
        wait
        cd ../

        # merge meryl dbs
        meryl union-sum output ~{identifier}.meryl extracted/*

        # package
        tar cvf ~{identifier}.meryl.tar ~{identifier}.meryl

        # cleanup
        rm -rf extracted/
	>>>
	output {
		File merylDb = identifier + ".meryl.tar"
	}
    runtime {
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}


task merylHapmer {
    input {
        File maternalMerylDB
        File paternalMerylDB
        File sampleMerylDB
        Int memSize = 24
        Int threadCount = 32
        Int diskSizeGB = 64
        String dockerImage = "juklucas/hpp_merqury:latest"
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
        OMP_NUM_THREADS=~{threadCount}

        # extract meryl dbs
        tar xvf ~{maternalMerylDB} &
        tar xvf ~{paternalMerylDB} &
        tar xvf ~{sampleMerylDB} &
        wait

        # generate hapmers
        bash $MERQURY/trio/hapmers.sh maternal.meryl paternal.meryl sample.meryl

        # package images
        tar czvf hapmers_img.tar.gz *.png *.hist

        # our desired destination is a softlink, need to move files to a real directory to tar
        mv maternal.hapmer.meryl maternal.tmp
        mkdir maternal.hapmers.meryl
        mv maternal.tmp/* maternal.hapmers.meryl/
        mv paternal.hapmer.meryl paternal.tmp
        mkdir paternal.hapmers.meryl
        mv paternal.tmp/* paternal.hapmers.meryl/
        tar cvf maternal.hapmers.meryl.tar maternal.hapmers.meryl &
        tar cvf paternal.hapmers.meryl.tar paternal.hapmers.meryl &
        wait
	>>>
	output {
		File maternalHapmers = "maternal.hapmers.meryl.tar"
		File paternalHapmers = "paternal.hapmers.meryl.tar"
		File hapmerImages = "hapmers_img.tar.gz"
	}
    runtime {
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}

