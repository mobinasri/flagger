version 1.0

workflow runAssemblyCleanUp{
    call assemblyCleanUp
    output {
        File mitoExemplarFastaGz = assemblyCleanUp.mitoExemplarFastaGz
        File mitoFastaGz = assemblyCleanUp.mitoFastaGz
        File cleanedAssemblyFastaGz = assemblyCleanUp.cleanedAssemblyFastaGz
    }
}


task assemblyCleanUp {
    input{
        File haplotype1GfaTarGz
        File haplotype2GfaTarGz
        String childID
        String suffix
        # runtime configurations
        Int memSizeGB=32
        Int threadCount=8
        Int diskSizeGB=128
        Int preemptible=1
        String dockerImage="quay.io/biocontainers/verkko:1.4--h48217b1_0"
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

        mkdir outputs
        tar -xvzf ~{haplotype1GfaTarGz} --strip-components 1 --directory outputs
        tar -xvzf ~{haplotype2GfaTarGz} --strip-components 1 --directory outputs
        cd outputs 

        cat *.hap[12].p_ctg.gfa | awk '{if (match($1, "^S")) { print ">"$2; print $3}}' > combined.fasta
        cat *.hap[12].p_ctg.gfa > combined.gfa
        cat *.hap[12].p_ctg.noseq.gfa | awk '{if (match($1, "^S")) { print "path "$2"\t"$2; print "piece000002"; print "end"; }}' > combined.scfmap
        cat *.hap[12].p_ctg.noseq.gfa | awk -v NAME="" '{if (match($1, "^S")) { if (NAME != "") print NAME"\t"SUM/LEN; SUM=0; NAME=$2; LEN=substr($4, 6, length($4)); } if (match($1, "^A")) { SUM+=$7-$6}} END {print NAME"\t"SUM/LEN}' > combined.csv

        # Fix error:
        # mashmap: error while loading shared libraries: libmkl_rt.so.2: cannot open shared object file: No such file or directory
        ## Seen in 1.4--h48217b1_0 Biocontainer for Verkko
        pip install mkl

        ## Download newer version of script that doesn't miss as many contigs
        ## lowers segment length in mashmap to 5kb (from 10kb) and enforces more overlap for circularization
        wget https://raw.githubusercontent.com/marbl/verkko/ff64f65a1cfebbd70207f17ac80723ecd3dc9030/src/scripts/screen-assembly.pl
        chmod +x screen-assembly.pl

        ./screen-assembly.pl \
          --assembly      combined.fasta \
          --graph         combined.gfa \
          --graphmap      combined.scfmap \
          --hifi-coverage combined.csv \
          --minlength     100000 \
          --contaminant   ebv /usr/local/lib/verkko/data/human-ebv-AJ507799.2.fasta.gz mito /usr/local/lib/verkko/data/human-mito-NC_012920.1.fasta.gz rdna /usr/local/lib/verkko/data/human-rdna-KY962518.1.fasta.gz \
          --output        assembly \
          --mashmap       mashmap \
          --threads       ~{threadCount} \
          > ./screen-assembly.out 2> ./screen-assembly.err

          #  Try to circularize, the script will output to the same file as input since it caches the sequences
          #  if it fails to circularize it will output sequence unchanged
          for xx in assembly.*.exemplar.fasta ; do
              if [ -e $xx ]; then
                  python3 /usr/local/lib/verkko/scripts/circularize_ctgs.py -f 0.90 -p 100 -o $xx --min-ovl 2500 $xx
              fi
          done
          cd ../
         
          mv outputs/assembly.fasta ~{childID}.~{suffix}.cleaned.fa
          gzip ~{childID}.~{suffix}.cleaned.fa
 
          mv outputs/assembly.mito.fasta ~{childID}.~{suffix}.mito.fa
          gzip ~{childID}.~{suffix}.mito.fa

          mv outputs/assembly.mito.exemplar.fasta ~{childID}.~{suffix}.mito.exemplar.fa
          gzip ~{childID}.~{suffix}.mito.exemplar.fa
          
    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible: preemptible
    }

    output {
        File mitoFastaGz = "${childID}.${suffix}.mito.fa.gz"
        File mitoExemplarFastaGz = "${childID}.${suffix}.mito.exemplar.fa.gz"
        File cleanedAssemblyFastaGz = "${childID}.${suffix}.cleaned.fa.gz"
    }
}

