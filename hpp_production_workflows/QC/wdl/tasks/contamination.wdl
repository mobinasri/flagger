version 1.0


workflow runContamination {

    input {
        String? assemblyIdentifier
        File assemblyFasta
        File? eukContaminationDatabase
        File? mitoContaminationDatabase
        File? plastidsContaminationDatabase
        File? vecscreenContaminationDatabase
        File? rrnaContaminationDatabase
        Array[File] refseqContaminationDatabases
        String dockerImage="tpesout/hpp_contamination:latest"
    }

    if (defined(eukContaminationDatabase)) {
        call contaminationEuk {
            input:
                assemblyFasta=assemblyFasta,
                eukContaminationDatabase=eukContaminationDatabase,
                dockerImage=dockerImage
        }
    }
    if (defined(mitoContaminationDatabase)) {
        call contaminationMito {
            input:
                assemblyFasta=assemblyFasta,
                mitoContaminationDatabase=mitoContaminationDatabase,
                dockerImage=dockerImage
        }
    }
    if (defined(plastidsContaminationDatabase)) {
        call contaminationPlastids {
            input:
                assemblyFasta=assemblyFasta,
                plastidsContaminationDatabase=plastidsContaminationDatabase,
                dockerImage=dockerImage
        }
    }
    if (defined(vecscreenContaminationDatabase)) {
        call contaminationVecscreen {
            input:
                assemblyFasta=assemblyFasta,
                vecscreenContaminationDatabase=vecscreenContaminationDatabase,
                dockerImage=dockerImage
        }
    }
    if (defined(rrnaContaminationDatabase)) {
        call contaminationRRNA {
            input:
                assemblyFasta=assemblyFasta,
                rrnaContaminationDatabase=rrnaContaminationDatabase,
                dockerImage=dockerImage
        }
    }
    if (length(refseqContaminationDatabases) > 0) {
        call contaminationWindowmasker {
            input:
                assemblyFasta=assemblyFasta,
                dockerImage=dockerImage
        }
    }
    scatter (refseqContaminationDatabase in refseqContaminationDatabases) {
        call contaminationRefseq {
            input:
                assemblyWindowmaskedFasta=contaminationWindowmasker.outputWindowmasker,
                refseqContaminationDatabase=refseqContaminationDatabase,
                dockerImage=dockerImage
        }
    }
    call mergeContaminationResults {
        input:
            assemblyIdentifier=assemblyIdentifier,
            assemblyFasta=assemblyFasta,
            eukOut=contaminationEuk.outputEuk,
            mitoOut=contaminationMito.outputMito,
            plastidsOut=contaminationPlastids.outputPlastids,
            rrnaOut=contaminationRRNA.outputRRNA,
            refseqOuts=contaminationRefseq.outputRefseq,
            vecscreenOut=contaminationVecscreen.outputVecscreen,
            dockerImage=dockerImage
    }
    call createContaminationBed {
        input:
            assemblyIdentifier=assemblyIdentifier,
            assemblyFasta=assemblyFasta,
            eukOut=contaminationEuk.outputEuk,
            mitoOut=contaminationMito.outputMito,
            plastidsOut=contaminationPlastids.outputPlastids,
            rrnaOut=contaminationRRNA.outputRRNA,
            refseqOuts=contaminationRefseq.outputRefseq,
            vecscreenOut=contaminationVecscreen.outputVecscreen,
            dockerImage=dockerImage
    }
	output {
	    File? eukContamination = contaminationEuk.outputEuk
	    File? mitoContamination = contaminationMito.outputMito
	    File? plastidsContamination = contaminationPlastids.outputPlastids
	    File? rrnaContamination = contaminationRRNA.outputRRNA
	    File? windowmaskedFasta=contaminationWindowmasker.outputWindowmasker
	    Array[File] refseqContamination = contaminationRefseq.outputRefseq
	    File? vecscreenContamination = contaminationVecscreen.outputVecscreen
	    File mergedResult = mergeContaminationResults.outputSummary
	    File fullBed = createContaminationBed.fullContaminationBed
	    File mergedBed = createContaminationBed.mergedContaminationBed
	    File contaminationCoverage = createContaminationBed.contaminationCoverage
	    File contaminationContigs = createContaminationBed.contaminationContigs
	}
}


task contaminationEuk {
    input {
        File assemblyFasta
        File? eukContaminationDatabase
        String eukExtraArguments=""
        Int memSizeGB = 8
        Int threadCount = 8
        Int diskSizeGB = 64
        String dockerImage = "tpesout/hpp_contamination:latest"
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

        # initilization
        ASM_FILENAME=$(basename -- "~{assemblyFasta}")
        if [[ $ASM_FILENAME =~ \.gz$ ]]; then
            cp ~{assemblyFasta} .
            gunzip $ASM_FILENAME
            ASM_FILENAME="${ASM_FILENAME%.gz}"
        else
            ln -s ~{assemblyFasta}
        fi
        PREFIX="${ASM_FILENAME%.*}"

        # make db index
        DB_FILENAME=$(basename -- "~{eukContaminationDatabase}")
        if [[ $DB_FILENAME =~ \.gz$ ]]; then
            cp ~{eukContaminationDatabase} .
            gunzip $DB_FILENAME
            DB_FILENAME="${DB_FILENAME%.gz}"
        else
            ln -s ~{eukContaminationDatabase}
        fi
        makeblastdb -in $DB_FILENAME -input_type fasta -dbtype nucl

        # blastoff
        time blastn -query $ASM_FILENAME \
            -db $DB_FILENAME \
            -task megablast \
            -num_threads ~{threadCount} \
            -word_size 28 \
            -best_hit_overhang 0.1 \
            -best_hit_score_edge 0.1 \
            -dust yes \
            -evalue 0.0001 \
            -outfmt '7' \
            -perc_identity 90.0  ~{eukExtraArguments} \
            | awk '($3>=98.0 && $4>=50)||($3>=94.0 && $4>=100)||($3>=90.0 && $4>=200)' \
            > $PREFIX.contam_in_euks.awk

        echo -e "#query acc.ver\tsubject acc.ver\t% id\tlength\tmsmatch\tgaps\tq start\tq end\ts start\ts end\tevalue\tscore" >$PREFIX.contam_in_euks.tmp
        cat $PREFIX.contam_in_euks.awk | { grep -v "^#" || true; } >>$PREFIX.contam_in_euks.tmp

        update_megablast_output_with_subject_description.py $DB_FILENAME $PREFIX.contam_in_euks.tmp $PREFIX.contam_in_euks.tsv

	>>>
	output {
		File outputEuk = glob("*.contam_in_euks.tsv")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}


task contaminationMito {
    input {
        File assemblyFasta
        File? mitoContaminationDatabase
        String mitoExtraArguments=""
        Int memSizeGB = 8
        Int threadCount = 8
        Int diskSizeGB = 64
        String dockerImage = "tpesout/hpp_contamination:latest"
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

        # initilization
        ASM_FILENAME=$(basename -- "~{assemblyFasta}")
        if [[ $ASM_FILENAME =~ \.gz$ ]]; then
            cp ~{assemblyFasta} .
            gunzip $ASM_FILENAME
            ASM_FILENAME="${ASM_FILENAME%.gz}"
        else
            ln -s ~{assemblyFasta}
        fi
        PREFIX="${ASM_FILENAME%.*}"

        # make db index
        DB_FILENAME=$(basename -- "~{mitoContaminationDatabase}")
        if [[ $DB_FILENAME =~ \.gz$ ]]; then
            cp ~{mitoContaminationDatabase} .
            gunzip $DB_FILENAME
            DB_FILENAME="${DB_FILENAME%.gz}"
        else
            ln -s ~{mitoContaminationDatabase}
        fi
        makeblastdb -in $DB_FILENAME -input_type fasta -dbtype nucl

        # blastoff
        time blastn -query $ASM_FILENAME \
            -db $DB_FILENAME \
            -task megablast \
            -num_threads ~{threadCount} \
            -word_size 28 \
            -best_hit_overhang 0.1 \
            -best_hit_score_edge 0.1 \
            -dust yes \
            -evalue 0.0001 \
            -perc_identity 90 \
            -soft_masking true \
            -outfmt 7 ~{mitoExtraArguments} \
            | awk '$4>=500' \
            > $PREFIX.mito.awk

        echo -e "#query acc.ver\tsubject acc.ver\t% id\tlength\tmsmatch\tgaps\tq start\tq end\ts start\ts end\tevalue\tscore" >$PREFIX.mito.tmp
        cat $PREFIX.mito.awk | { grep -v "^#" || true; } >>$PREFIX.mito.tmp

        update_megablast_output_with_subject_description.py $DB_FILENAME $PREFIX.mito.tmp $PREFIX.mito.tsv
	>>>
	output {
		File outputMito = glob("*.mito.tsv")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}


task contaminationPlastids {
    input {
        File assemblyFasta
        File? plastidsContaminationDatabase
        String plastidsExtraArguments=""
        Int memSizeGB = 8
        Int threadCount = 8
        Int diskSizeGB = 64
        String dockerImage = "tpesout/hpp_contamination:latest"
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

        # initilization
        ASM_FILENAME=$(basename -- "~{assemblyFasta}")
        if [[ $ASM_FILENAME =~ \.gz$ ]]; then
            cp ~{assemblyFasta} .
            gunzip $ASM_FILENAME
            ASM_FILENAME="${ASM_FILENAME%.gz}"
        else
            ln -s ~{assemblyFasta}
        fi
        PREFIX="${ASM_FILENAME%.*}"

        # make db index
        DB_FILENAME=$(basename -- "~{plastidsContaminationDatabase}")
        if [[ $DB_FILENAME =~ \.gz$ ]]; then
            cp ~{plastidsContaminationDatabase} .
            gunzip $DB_FILENAME
            DB_FILENAME="${DB_FILENAME%.gz}"
        else
            ln -s ~{plastidsContaminationDatabase}
        fi
        makeblastdb -in $DB_FILENAME -input_type fasta -dbtype nucl

        # blastoff
        time blastn -query $ASM_FILENAME \
            -db $DB_FILENAME \
            -task megablast \
            -num_threads ~{threadCount} \
            -word_size 28 \
            -best_hit_overhang 0.1 \
            -best_hit_score_edge 0.1 \
            -dust yes \
            -evalue 0.0001 \
            -perc_identity 90 \
            -soft_masking true \
            -outfmt 7 ~{plastidsExtraArguments} \
            | awk '$4>=500' \
            > $PREFIX.plastids.awk

        echo -e "#query acc.ver\tsubject acc.ver\t% id\tlength\tmsmatch\tgaps\tq start\tq end\ts start\ts end\tevalue\tscore" >$PREFIX.plastids.tmp
        cat $PREFIX.plastids.awk | { grep -v "^#" || true; } >>$PREFIX.plastids.tmp

        update_megablast_output_with_subject_description.py $DB_FILENAME $PREFIX.plastids.tmp $PREFIX.plastids.tsv
	>>>
	output {
		File outputPlastids = glob("*.plastids.tsv")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}


task contaminationRRNA {
    input {
        File assemblyFasta
        File? rrnaContaminationDatabase
        Int chunkSize = 500000
        Int chunkOverlap = 5000
        String rrnaExtraArguments=""
        Int memSizeGB = 8
        Int threadCount = 8
        Int diskSizeGB = 64
        String dockerImage = "tpesout/hpp_contamination:latest"
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

        # initilization
        ASM_FILENAME=$(basename -- "~{assemblyFasta}")
        if [[ $ASM_FILENAME =~ \.gz$ ]]; then
            cp ~{assemblyFasta} .
            gunzip $ASM_FILENAME
            ASM_FILENAME="${ASM_FILENAME%.gz}"
        else
            ln -s ~{assemblyFasta}
        fi
        PREFIX="${ASM_FILENAME%.*}"

        # chunk it
        chunk_assembly.py $ASM_FILENAME $ASM_FILENAME.chunked ~{chunkSize} ~{chunkOverlap}

        # make db index
        DB_FILENAME=$(basename -- "~{rrnaContaminationDatabase}")
        if [[ $DB_FILENAME =~ \.gz$ ]]; then
            cp ~{rrnaContaminationDatabase} .
            gunzip $DB_FILENAME
            DB_FILENAME="${DB_FILENAME%.gz}"
        else
            ln -s ~{rrnaContaminationDatabase}
        fi
        makeblastdb -in $DB_FILENAME -input_type fasta -dbtype nucl

        # blastoff
        time blastn -query $ASM_FILENAME.chunked \
            -db $DB_FILENAME \
            -task megablast \
            -num_threads ~{threadCount} \
            -template_length 18 \
            -template_type coding \
            -window_size 120 \
            -word_size 12 \
            -xdrop_gap 20 \
            -no_greedy \
            -best_hit_overhang 0.1 \
            -best_hit_score_edge 0.1 \
            -dust yes \
            -evalue 1E-9 \
            -gapextend 2 \
            -gapopen 4 \
            -penalty -5 \
            -perc_identity 93 \
            -reward 2 \
            -soft_masking true \
            -outfmt 7 ~{rrnaExtraArguments} \
            | awk '$4>=100' \
            > $PREFIX.rrna.awk.chunked

        # unchunk it
        unchunk_blast.py $PREFIX.rrna.awk.chunked $PREFIX.rrna.awk 250000

        echo -e "#query acc.ver\tsubject acc.ver\t% id\tlength\tmsmatch\tgaps\tq start\tq end\ts start\ts end\tevalue\tscore" >$PREFIX.rrna.tmp
        cat $PREFIX.rrna.awk | { grep -v "^#" || true; } >>$PREFIX.rrna.tmp

        update_megablast_output_with_subject_description.py $DB_FILENAME $PREFIX.rrna.tmp $PREFIX.rrna.tsv

	>>>
	output {
		File outputRRNA = glob("*.rrna.tsv")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}


task contaminationWindowmasker {
    input {
        File assemblyFasta
        Int memSizeGB = 8
        Int threadCount = 8
        Int diskSizeGB = 64
        String dockerImage = "tpesout/hpp_contamination:latest"
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

        # initilization
        ASM_FILENAME=$(basename -- "~{assemblyFasta}")
        if [[ $ASM_FILENAME =~ \.gz$ ]]; then
            cp ~{assemblyFasta} .
            gunzip $ASM_FILENAME
            ASM_FILENAME="${ASM_FILENAME%.gz}"
        else
            ln -s ~{assemblyFasta}
        fi
        PREFIX="${ASM_FILENAME%.*}"

        windowmasker -mk_counts -in $ASM_FILENAME -out $PREFIX.stage1
        windowmasker -ustat $PREFIX.stage1 -in $ASM_FILENAME -out $PREFIX.masked.fa -dust true -outfmt fasta

	>>>
	output {
		File outputWindowmasker = glob("*.masked.fa")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}

task contaminationRefseq {
    input {
        File? assemblyWindowmaskedFasta
        File? refseqContaminationDatabase
        String refseqExtraArguments=""
        Int memSizeGB = 8
        Int threadCount = 8
        Int diskSizeGB = 128
        String dockerImage = "tpesout/hpp_contamination:latest"
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

        # initilization
        ASM_FILENAME=$(basename -- "~{assemblyWindowmaskedFasta}")
        ln -s ~{assemblyWindowmaskedFasta}
        PREFIX="${ASM_FILENAME%.masked.fa}"
        RSDB_FILE=$(basename -- "~{refseqContaminationDatabase}")
        RSDB_FILE_ID="${RSDB_FILE%.*}"
        RSDB_ID="${RSDB_FILE%%.*}"

        # make db index
        DB_FILENAME=$(basename -- "~{refseqContaminationDatabase}")
        if [[ $DB_FILENAME =~ \.gz$ ]]; then
            cp ~{refseqContaminationDatabase} .
            gunzip $DB_FILENAME
            DB_FILENAME="${DB_FILENAME%.gz}"
        else
            ln -s ~{refseqContaminationDatabase}
        fi
        makeblastdb -in $DB_FILENAME -input_type fasta -dbtype nucl

        # blastoff
        blastn -query $ASM_FILENAME \
            -db $DB_FILENAME \
            -task megablast \
            -num_threads ~{threadCount} \
            -word_size 28 \
            -best_hit_overhang 0.1 \
            -best_hit_score_edge 0.1 \
            -dust yes \
            -evalue 0.0001 \
            -min_raw_gapped_score 100 \
            -penalty -5 \
            -perc_identity 98 \
            -soft_masking true \
            -outfmt 7 ~{refseqExtraArguments} \
            > $PREFIX.refseq.$RSDB_FILE_ID.out
#            -negative_taxidlist ${negative_tax_id_list} \

        echo -e "========== RefSeq: $RSDB_FILE_ID ==========" >$PREFIX.refseq.$RSDB_FILE_ID.tmp
        echo -e "#query acc.ver\tsubject acc.ver\t% id\tlength\tmsmatch\tgaps\tq start\tq end\ts start\ts end\tevalue\tscore" >>$PREFIX.refseq.$RSDB_FILE_ID.tmp
        cat $PREFIX.refseq.$RSDB_FILE_ID.out | { grep -v "^#" || true; } >>$PREFIX.refseq.$RSDB_FILE_ID.tmp

        update_megablast_output_with_subject_description.py $DB_FILENAME $PREFIX.refseq.$RSDB_FILE_ID.tmp $PREFIX.refseq.$RSDB_FILE_ID.tsv $RSDB_ID

	>>>
	output {
		File outputRefseq = glob("*.refseq.*.tsv")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}



task contaminationVecscreen {
    input {
        File assemblyFasta
        File? vecscreenContaminationDatabase
        String vecscreenExtraArguments=""
        Int chunkSize = 1000000
        Int chunkOverlap = 10000
        Int memSizeGB = 8
        Int threadCount = 8
        Int diskSizeGB = 64
        String dockerImage = "tpesout/hpp_contamination:latest"
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

        # initilization
        ASM_FILENAME=$(basename -- "~{assemblyFasta}")
        if [[ $ASM_FILENAME =~ \.gz$ ]]; then
            cp ~{assemblyFasta} .
            gunzip $ASM_FILENAME
            ASM_FILENAME="${ASM_FILENAME%.gz}"
        else
            ln -s ~{assemblyFasta}
        fi
        PREFIX="${ASM_FILENAME%.*}"

        # chunk it (vecscreen cannot handle long sequences)
        chunk_assembly.py $ASM_FILENAME $ASM_FILENAME.chunked ~{chunkSize} ~{chunkOverlap}

        # index it (old version of blast.. don't know if we NEED 2.7.1 but the latest doesn't work)
        DB_FILENAME=$(basename -- "~{vecscreenContaminationDatabase}")
        if [[ $DB_FILENAME =~ \.gz$ ]]; then
            cp ~{vecscreenContaminationDatabase} .
            gunzip $DB_FILENAME
            DB_FILENAME="${DB_FILENAME%.gz}"
        else
            ln -s ~{vecscreenContaminationDatabase}
        fi
        /opt/blast/ncbi-blast-2.7.1+/bin/makeblastdb -in $DB_FILENAME -input_type fasta -dbtype nucl

        # screen it
        vecscreen -d $DB_FILENAME \
            -f3 \
            -i $ASM_FILENAME.chunked \
            -o $PREFIX.vecscreen.out ~{vecscreenExtraArguments}
        VSlistTo1HitPerLine.awk \
            suspect=0 \
            weak=0 \
            $PREFIX.vecscreen.out \
            > $PREFIX.vecscreen.filtered.out

        # unchunk it
        unchunk_vecscreen.py $PREFIX.vecscreen.filtered.out $PREFIX.vecscreen.tsv ~{chunkSize}

	>>>
	output {
		File outputVecscreen = glob("*.vecscreen.tsv")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}


task mergeContaminationResults {
    input {
        String? assemblyIdentifier
        File assemblyFasta
        File? eukOut
        File? mitoOut
        File? rrnaOut
        File? plastidsOut
        Array[File] refseqOuts
        File? rrnaOut
        File? vecscreenOut
        Int memSizeGB = 2
        Int threadCount = 1
        Int diskSizeGB = 64
        String dockerImage = "tpesout/hpp_contamination:latest"
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

        if [[ -z "~{assemblyIdentifier}" ]] ; then
            OUT="$(basename ~{assemblyFasta} | sed 's/.fa\(sta\)\?\(.gz\)\?//').contamination.txt"
        else
            OUT="~{assemblyIdentifier}.contamination.txt"
        fi
        touch $OUT

        # contam in euks
        if [[ -f "~{eukOut}" ]] ; then
            echo "========== Eukaryote ==========" >>$OUT
            cat ~{eukOut} >>$OUT
        fi

        # mitochondria
        if [[ -f "~{mitoOut}" ]] ; then
            echo "" >>$OUT
            echo "========== Mitocondria ==========" >>$OUT
            cat ~{mitoOut} >>$OUT
        fi

        # plastids
        if [[ -f "~{plastidsOut}" ]] ; then
            echo "" >>$OUT
            echo "========== Plastids ==========" >>$OUT
            cat ~{plastidsOut} >>$OUT
        fi

        # ribosomal rna
        if [[ -f "~{rrnaOut}" ]] ; then
            echo "" >>$OUT
            echo "========== RRNA ==========" >>$OUT
            cat ~{rrnaOut} >>$OUT
        fi

        # vecscreen
        if [[ -f "~{vecscreenOut}" ]] ; then
            echo "" >>$OUT
            echo "========== Vecscreen ==========" >>$OUT
            cat ~{vecscreenOut} >>$OUT
        fi

        # refseq
        if [[ -z "~{true="" false="NO_REFSEQ" length(refseqOuts) > 0}" ]] ; then
            echo "" >>$OUT
            echo "========== RefSeq ==========" >>$OUT
            cat ~{sep=" " refseqOuts} >>$OUT
        fi

	>>>
	output {
		File outputSummary = glob("*.contamination.txt")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}


task createContaminationBed {
    input {
        String? assemblyIdentifier
        File assemblyFasta
        File? eukOut
        File? mitoOut
        File? rrnaOut
        File? plastidsOut
        Array[File] refseqOuts
        File? rrnaOut
        File? vecscreenOut
        Int bedtoolsMergeDistance = 0
        Float discardContigCoverageThreshold = 0.95
        Int memSizeGB = 2
        Int threadCount = 1
        Int diskSizeGB = 64
        String dockerImage = "tpesout/hpp_contamination:latest"
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

        # get filenames
        if [[ -z "~{assemblyIdentifier}" ]] ; then
            ASM_ID=`basename ~{assemblyFasta} | sed 's/.fa\(sta\)\?\(.gz\)\?//'`
        else
            ASM_ID="~{assemblyIdentifier}"
        fi
        OUT="${ASM_ID}.full_contamination.bed"

        # output header
        echo -e "#contig\tstart\tstop\tcontamination_source\tcontamination_contig\tcontamination_description" >$OUT

        # save all data (if there is a file for it)
        # blast output is 1-based, bed is 0-based, so we subtract 1 from start pos

        # contam in euks
        if [[ -f "~{eukOut}" ]] ; then
            cat ~{eukOut}       | { grep -v "^#" || true; } | awk -F "\t" '{print $1 "\t" $7 - 1 "\t" $8 "\tcontam_in_euk\t" $2 "\t" $13}' >>$OUT
        fi

        # mito
        if [[ -f "~{mitoOut}" ]] ; then
            cat ~{mitoOut}      | { grep -v "^#" || true; } | awk -F "\t" '{print $1 "\t" $7 - 1 "\t" $8 "\tmitochondria\t" $2 "\t" $13}' >>$OUT
        fi

        # plastids
        if [[ -f "~{plastidsOut}" ]] ; then
            cat ~{plastidsOut}  | { grep -v "^#" || true; } | awk -F "\t" '{print $1 "\t" $7 - 1 "\t" $8 "\tplastids\t" $2 "\t" $13}' >>$OUT
        fi

        # ribosomal RNA
        if [[ -f "~{rrnaOut}" ]] ; then
            cat ~{rrnaOut}      | { grep -v "^#" || true; } | awk -F "\t" '{print $1 "\t" $7 - 1 "\t" $8 "\trrna\t" $2 "\t" $13}' >>$OUT
        fi

        # vecscreen
        if [[ -f "~{vecscreenOut}" ]] ; then
            cat ~{vecscreenOut} | { grep -v "^#" || true; } | awk -F "\t" '{print $1 "\t" $2 - 1 "\t" $3 "\tvecscreen\t" $4 "\t" $4}' >>$OUT
        fi

        # refseq
        if [[ -z "~{true="" false="NO_REFSEQ" length(refseqOuts) > 0}" ]] ; then
            cat ~{sep=" " refseqOuts} | { grep -v "^#" || true; } | { grep -v "==========" || true; } | awk -F "\t" '{print $1 "\t" $7 - 1 "\t" $8 "\trefseq\t" $2 "\t" $13}' >>$OUT
        fi

        # merge results (unless we don't have results)
        cat $OUT | { grep -v "^$" || true; } | bedtools sort >tmp
        ls -lah
        if [[ -s "tmp" ]] ; then
            bedtools merge -d ~{bedtoolsMergeDistance} -delim ";" -c 4,5,6 -o distinct,distinct,distinct -i tmp > ${ASM_ID}.merged_contamination.bed

            # genome coverage
            ASM=$(basename ~{assemblyFasta})
            if [[ ~{assemblyFasta} =~ \.gz$ ]]; then
                cp ~{assemblyFasta} .
                gunzip $ASM
                ASM="${ASM%.gz}"
            else
                ln -s ~{assemblyFasta} .
            fi
            samtools faidx $ASM
            bedtools genomecov -i ${ASM_ID}.merged_contamination.bed -g $ASM.fai | awk '$2 == 1 {print $1 "\t" $4 "\t" $5}' | sort -nr -k 3 >tmp

            # annotate coverage with types of matches
            while read LINE; do
                CTG=$(echo "$LINE" | cut -f1)
                MATCHES=$(cat ${ASM_ID}.merged_contamination.bed | awk -v CTG=$CTG '$1 == CTG' | cut -f4 | sed 's/;/\n/' | sort | uniq | sed -z 's/\n/,/g;s/,$//')
                echo -e "$LINE\t$MATCHES" >>${ASM_ID}.contamination_coverage.tsv
            done <tmp

            # annotate contigs to delete
            cat ${ASM_ID}.contamination_coverage.tsv | awk '$3 >= ~{discardContigCoverageThreshold}' | cut -f1 >${ASM_ID}.contamination_contigs.txt
        else
            touch ${ASM_ID}.merged_contamination.bed
            touch ${ASM_ID}.contamination_coverage.tsv
            touch ${ASM_ID}.contamination_contigs.txt
        fi

	>>>
	output {
		File fullContaminationBed = glob("*.full_contamination.bed")[0]
		File mergedContaminationBed = glob("*.merged_contamination.bed")[0]
		File contaminationCoverage = glob("*.contamination_coverage.tsv")[0]
		File contaminationContigs = glob("*.contamination_contigs.txt")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
