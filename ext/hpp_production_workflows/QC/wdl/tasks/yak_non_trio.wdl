version 1.0

import "extract_reads.wdl" as extractReads_t
import "shard_reads.wdl" as shardReads_t
import "arithmetic.wdl" as arithmetic_t

workflow runNonTrioYakAssemblyStats {

    input {
        Array[File] sampleReadsILM
        File assemblyFastaHap1
        File assemblyFastaHap2        
        File? referenceFasta
        Int shardLinesPerFile = 256000000
        Int fileExtractionDiskSizeGB = 256
        String dockerImage = "juklucas/hpp_yak:latest"
    }

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

    call arithmetic_t.sum as sampleReadSize {
        input:
            integers=sampleReadsExtracted.fileSizeGB
    }

    call yakCount as yakCountSample {
        input:
            readFiles=sampleReadsExtracted.extractedRead,
            sampleName="sample",
            diskSizeGB=sampleReadSize.value * 2,
            dockerImage=dockerImage
    }

    # get stats
    call yakNonTrioAssemblyStats {
        input:
            assemblyFastaHap1=assemblyFastaHap1,
            assemblyFastaHap2=assemblyFastaHap2,            
            sampleYak=yakCountSample.outputYak,
            dockerImage=dockerImage
    }

	output {
		File outputTarball = yakNonTrioAssemblyStats.outputTarball
		File outputSummary = yakNonTrioAssemblyStats.outputSummary
	}

}


task yakCount {
    input{
        Array[File] readFiles
        String sampleName
        Int bloomSize=37
        # runtime configurations
        Int memSize=128
        Int threadCount=16
        Int diskSizeGB=256
        String dockerImage="juklucas/hpp_yak:latest"
    }
    command <<<
        set -eux -o xtrace

        # Kmer counting with https://github.com/lh3/yak.
        yak count -t~{threadCount} -b~{bloomSize} -o ~{sampleName}.yak <(cat ~{sep=" " readFiles}) <(cat ~{sep=" " readFiles})
    >>>

    runtime {
        docker: dockerImage
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible: 1
    }

    output {
        File outputYak = "~{sampleName}.yak"
    }
}


task yakNonTrioAssemblyStats {
    input {
        File assemblyFastaHap2
        File assemblyFastaHap1
        File sampleYak
        String genomeSize = "3.2g"
        String minSequenceLength = "100k"
        # runtime configurations
        Int memSize = 128
        Int threadCount = 32
        Int diskSizeGB = 256
        String dockerImage = "juklucas/hpp_yak:latest"
    }
    command <<<
        set -eux -o xtrace

        # name
        PREFIX=$(basename ~{assemblyFastaHap2} | sed 's/.gz$//' | sed 's/.fa\(sta\)*$//' | sed 's/[._][pm]at\(ernal\)*//')

        # QV
        yak qv -t ~{threadCount} -p -K ~{genomeSize} -l ~{minSequenceLength} ~{sampleYak} ~{assemblyFastaHap1} > $PREFIX.hap1.yak.qv.txt
        yak qv -t ~{threadCount} -p -K ~{genomeSize} -l ~{minSequenceLength} ~{sampleYak} ~{assemblyFastaHap2} > $PREFIX.hap2.yak.qv.txt

        # condense
        SUMMARY=$PREFIX.summary.txt
        echo "# hap1 qv" >>$SUMMARY
        tail -n4 $PREFIX.hap1.yak.qv.txt >>$SUMMARY
        echo "# hap2 qv" >>$SUMMARY
        tail -n4 $PREFIX.hap2.yak.qv.txt >>$SUMMARY

        # tar
        tar czvf $PREFIX.yak-qc.tar.gz *txt
    >>>

    runtime {
        docker: dockerImage
        memory: memSize
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible: 1
    }

    output {
        File outputTarball = glob("*.yak-qc.tar.gz")[0]
        File outputSummary = glob("*.summary.txt")[0]
    }
}
