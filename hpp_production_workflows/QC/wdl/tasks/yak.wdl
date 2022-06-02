version 1.0

import "extract_reads.wdl" as extractReads_t
import "shard_reads.wdl" as shardReads_t
import "arithmetic.wdl" as arithmetic_t

workflow runYakAssemblyStats {

    input {
        Array[File] maternalReadsILM
        Array[File] paternalReadsILM
        Array[File] sampleReadsILM
        File assemblyFastaPat
        File assemblyFastaMat
        File? referenceFasta
        Int shardLinesPerFile = 256000000
        Int fileExtractionDiskSizeGB = 256
        String dockerImage = "juklucas/hpp_yak:latest"
    }

    # extract reads
    scatter (readFile in maternalReadsILM) {
        call extractReads_t.extractReads as maternalReadsExtracted {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSizeGB=4,
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
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=fileExtractionDiskSizeGB,
                dockerImage=dockerImage
        }
    }
    scatter (readFile in sampleReadsILM) {
        call extractReads_t.extractReads as sampleReadsExtracted {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=fileExtractionDiskSizeGB,
                dockerImage=dockerImage
        }
    }

    # get file size of results (for yak counting)
    call arithmetic_t.sum as maternalReadSize {
        input:
            integers=maternalReadsExtracted.fileSizeGB
    }
    call arithmetic_t.sum as paternalReadSize {
        input:
            integers=paternalReadsExtracted.fileSizeGB
    }
    call arithmetic_t.sum as sampleReadSize {
        input:
            integers=sampleReadsExtracted.fileSizeGB
    }

    # do counting
    call yakCount as yakCountMat {
        input:
            readFiles=maternalReadsExtracted.extractedRead,
            sampleName="mat",
            diskSizeGB=maternalReadSize.value * 2,
            dockerImage=dockerImage
    }
    call yakCount as yakCountPat {
        input:
            readFiles=paternalReadsExtracted.extractedRead,
            sampleName="pat",
            diskSizeGB=paternalReadSize.value * 2,
            dockerImage=dockerImage
    }
    call yakCount as yakCountSample {
        input:
            readFiles=sampleReadsExtracted.extractedRead,
            sampleName="sample",
            diskSizeGB=sampleReadSize.value * 2,
            dockerImage=dockerImage
    }

    # get stats
    call yakAssemblyStats {
        input:
            assemblyFastaPat=assemblyFastaPat,
            assemblyFastaMat=assemblyFastaMat,
            patYak=yakCountPat.outputYak,
            matYak=yakCountMat.outputYak,
            sampleYak=yakCountSample.outputYak,
            dockerImage=dockerImage
    }

	output {
		File outputTarball = yakAssemblyStats.outputTarball
		File outputSummary = yakAssemblyStats.outputSummary
		File maternalYak = yakCountMat.outputYak
		File paternalYak = yakCountPat.outputYak
	}

}


task yakCount {
    input{
        Array[File] readFiles
        String sampleName
        Int bloomSize=37
        # runtime configurations
        Int memSizeGB=128
        Int threadCount=16
        Int diskSizeGB=256
        String dockerImage="juklucas/hpp_yak:latest"
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

        # Kmer counting with https://github.com/lh3/yak.
        yak count -t~{threadCount} -b~{bloomSize} -o ~{sampleName}.yak <(cat ~{sep=" " readFiles}) <(cat ~{sep=" " readFiles})
    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible: 1
    }

    output {
        File outputYak = "~{sampleName}.yak"
    }
}


task yakAssemblyStats {
    input {
        File assemblyFastaPat
        File assemblyFastaMat
        File patYak
        File matYak
        File sampleYak
        String genomeSize = "3.2g"
        String minSequenceLength = "100k"
        # runtime configurations
        Int memSizeGB = 128
        Int threadCount = 32
        Int diskSizeGB = 256
        String dockerImage = "juklucas/hpp_yak:latest"
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

        # name
        PREFIX=$(basename ~{assemblyFastaPat} | sed 's/.gz$//' | sed 's/.fa\(sta\)*$//' | sed 's/[._][pm]at\(ernal\)*//')

        # Computing error rates
        yak trioeval -t ~{threadCount} ~{patYak} ~{matYak} ~{assemblyFastaPat} > $PREFIX.pat.yak.switch-error.txt
        yak trioeval -t ~{threadCount} ~{patYak} ~{matYak} ~{assemblyFastaMat} > $PREFIX.mat.yak.switch-error.txt

        # QV
        yak qv -t ~{threadCount} -p -K ~{genomeSize} -l ~{minSequenceLength} ~{sampleYak} ~{assemblyFastaPat} > $PREFIX.pat.yak.qv.txt
        yak qv -t ~{threadCount} -p -K ~{genomeSize} -l ~{minSequenceLength} ~{sampleYak} ~{assemblyFastaMat} > $PREFIX.mat.yak.qv.txt

        # condense
        SUMMARY=$PREFIX.summary.txt
        echo "# mat qv" >>$SUMMARY
        tail -n4 $PREFIX.mat.yak.qv.txt >>$SUMMARY
        echo "# pat qv" >>$SUMMARY
        tail -n4 $PREFIX.pat.yak.qv.txt >>$SUMMARY
        echo "# mat switch" >>$SUMMARY
        tail -n3 $PREFIX.mat.yak.switch-error.txt >>$SUMMARY
        echo "# pat switch" >>$SUMMARY
        tail -n3 $PREFIX.pat.yak.switch-error.txt >>$SUMMARY

        # tar
        tar czvf $PREFIX.yak-qc.tar.gz *txt
    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible: 1
    }

    output {
        File outputTarball = glob("*.yak-qc.tar.gz")[0]
        File outputSummary = glob("*.summary.txt")[0]
    }
}
