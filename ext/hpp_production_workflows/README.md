# HPP Production Workflows

This repository holds WDL workflows and Docker build scripts for 
production workflows used by the [Human Pangenome Project](https://humanpangenome.org/).

All WDLs and containers created in this repository are licensed under the MIT license. The underlying tools (that the WDLs and containers run) are likely covered under one or more Free and Open Source Software licenses, but we cannot make any guarantees to that fact.

### Organization
The repo is split into assembly and QC workflows with the following organization:

```
 ── docker/
    └── toolName/
        └── Dockerfile
        └── Makefile
        └── scripts/
            └── toolName/
                └── scriptName.py
 ── wdl/
    └── tasks/
    │   └── taskName.wdl
    └── workflows/ 
        └── workFlowName.wdl
```

------------------


## Assembly

Assemblies are produced with Hifiasm from reads filtered with HiFiAdapterFilt:
* Filter HiFi Data w/ [HiFiAdapterFilt](https://github.com/sheinasim/HiFiAdapterFilt/tree/master/DB)
* Run [Hifiasm](https://github.com/chhylp123/hifiasm)

------------------

## QC

### Automated Assembly QC

Assembly QC is run with one of two workflows:
* standard_qc.wdl
* standard_qc_no_qv.wdl (for samples with no child Illumina data)

The following tools are included:
* [asmgene](https://github.com/lh3/minimap2)
* [dipcall v0.1](https://github.com/lh3/dipcall/tree/v0.1)
* [dipcall v0.2](https://github.com/lh3/dipcall/tree/v0.2)
* [merqury](https://github.com/marbl/merqury)
* [QUAST](https://sourceforge.net/projects/quast/files/)
* [YAK](https://github.com/lh3/yak)

### Miscellaneous QC

In addition to the automated QC of assemblies, the following tasks are included:
* read_stats.wdl: outputs read length statistics for bam/fastq files
* mask_assembly.wdl: masks assemblies in regions found in bed file
* dropFastaContigs.wdl: drops contigs from assemblies based on text file list of contigs to drop
* contamination.wdl: detects contamination in assemblies

#### Contamination

The contamination workflow is modeled after an NCBI [document](https://https.ncbi.nlm.nih.gov/tools/vecscreen/contam/) 
describing methods for detection and flagging of common
contaminants in assemblies.  

The workflow takes as input an assembly and
contamination reference sequences to be used in the following subtasks.
Each subtask is optional and is enabled by providing the reference sequence.
* Eukaryotic Contamination: cloning artifacts that are likely to show up 
as contaminants across all eukaryotic species: vector sequences, 
E.coli genome, phage genomes, bacterial Insertion Sequences and 
transposons [here](https://ftp.ncbi.nlm.nih.gov/pub/kitts/contam_in_euks.fa.gz)
* Mitochondria: mitochondrial sequences end in genomic.fna.gz [here](https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/)
* Plastids: plastid sequences end in genomic.fna.gz [here](https://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/)
* Vecscreen: sequencing adaptor screens [here](https://ftp.ncbi.nlm.nih.gov/pub/kitts/adaptors_for_screening_euks.fa)
* RRNA: ribosomal RNA [here](https://ftp.ncbi.nlm.nih.gov/pub/kitts/rrna.gz)
* RefSeq: reference sequences from other species organized 
[here](https://ftp.ncbi.nlm.nih.gov/refseq/release/) ending in genomic.fna.gz.
Large files may need to be split into smaller shards.

The output includes the direct output from each subtask, BED files 
summarizing all the genomic locations where contamination was found, and 
a file detailing the ratio of contamination sequence per assembly contig.

------------------ 


