## SecPhase: Phasing the long reads aligned to dual/diploid assemblies

### Motivation
When we align long reads (HiFi or ONT at this time!) to the dual/diploid assembly of the same sample, the reads coming from homozygous regions may
not be aligned to the correct haplotype. Other possible locations of such reads are usually reported in the secondary alignments by the aligner. Finding
the correct haplotype becomes even harder if our assemblies are erroneous. Breaks and indel errors in the assembly may mislead the aligner.

In the Figure below you can see an example of two haplotypes that are highly similar and different only in one base. One of these haplotypes is assembled
correctly and the other one has a long indel or a break point (could be split into two separate contigs). If we align a read from the first haplotype to
the diploid assembly, the aligner may align the read to the both haplotypes and report the one to the false haplotype as the primary alignment. This happens
because of the error in the assembly. Even if the assembly is perfect it may also happen for highly similar haplotypes since
the aligner chooses one haplotype randomly.
<img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/docs/phasing/images/phase_reads_1.png" width="700" height="400">

Given a set of secondary alignments (most of the time only one) and one primary alignment we want to make sure if
the primary one is to the correct haplotype and if it is not find the correct one among the secondary alignments.

### Approach

#### 1. Find the initial set of markers
One way to achieve this aim is to find the single-base markers that can navigate us toward the correct haplotype. An example of such marker is shown
in the top figure by green and red bars. All mismatched bases in the alignments of a single read form the initial set of candidate markers. After
projecting the intial markers into the coordinate of the read (instead of reference) it is possible to determine the match/mismatch status of a marker in
the other alignments of the same read. A marker that is a mismatch in all of the alignments is removed immediately.

#### 2. Filter markers
Among the remaining markers we perform two main filtering:

#### 2.1 Filter Markers within insertions

If a marker, which is a mismatch in at least one alignment, appears within an insertion in another alignment we remove it in this step. We remove it
because it can be either a misassembly on the haplotype that induces the insertion or an error on the read. So that marker can be misleading on either cases. In the figure below you can see an example of such a marker (shown with a blue bar).

<img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/docs/phasing/images/phase_reads_2.png" width="700" height="400">

#### 2.2 Filter markers with low BAQ

After filtering the markers within insertions the remaining markers are all either match or mismatch in any of the alignments.
The other issue is that sometimes the alignment around a marker is not reliable due to the errors in the reads or the assembly especially homopolymer run errors. To measure the reliability of the alignment around each marker we calculate
Base Alignment Quality ( BAQ ). This is an adjustment to the raw base quality of the marker and works by realiging the reads to where they were
already aligned to. This realignment is performed through a banded HMM whose parameters have to be tuned before hand.
There are three parameters that have to be tuned for each sequencing platform; gap opening probability, gap extenstion probability and the bandwidth of the HMM.

After tuning the parameters based on the platform (which is either HiFi or ONT here) we filtered the markers with BAQs lower than a specific threshold (20 for HiFi and 10 for ONT).

There are two points worth noting:
- BAQ cannot be larger than the raw base quality reported by the sequencer (or base caller)
- For better performance the bases far from the markers (farther than 500) are not included in the BAQ calculation.

More information about BAQ can be found in [this paper](https://academic.oup.com/bioinformatics/article/27/8/1157/227268). BAQ calculation is already implemented in htslib and [that implmentation](https://github.com/samtools/htslib/blob/9672589346459d675d62851d5b7b5f2e5c919076/probaln.c) has been imported to this pipeline.

<img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/docs/phasing/images/phase_reads_3.png" width="700" height="400">

#### 3. Select the best alignment using Marker Consistency Score

After finding the confident markers we calculate marker consistency score for each alignment. For calcultaing this score for each alignment we take the markers that appeared as mismatches on that alignment and take a summation of their BAQ values with a negative sign. After calcualating the marker consistency score for all primary and secondary alignments of the same read we select the one with the largest score (usually zero is the largest) and report it as the alignment to the correct haplotype. If the selected alignment is primary we do nothing. You can see an example of calculating marker consistency score in the figure above.

Two heuristics are applied for increasing specificity:
- If the selected alignment is having a score lower than a specific threshold (in this pipeline `-50`) it is not reported as the correct one.
- If the selected alignment is secondary its score should be at least `20` (`10` for ONT) units larger than the score of the primary alignment otherwise it is not reported as the correct one.
(These parameters will be tuned more systematically in the later versions)

### Running Secphase in marker mode

To run Secphase it is recommended to use the docker image `mobinasri/secphase:v0.4.2`.

Here are the parameters `secphase` can accept:
```
./secphase -h

Usage: secphase  -i <INPUT_BAM> -f <FASTA> 
Options:
         --inputBam, -i         Input BAM file
         --inputFasta, -f         Input FASTA file
         --inputVcf, -v         Input phased VCF file
         --variantBed, -B         Input BED file for subsetting phased variants
         --outDir, -o         Output dir for saving outputs [Default = "secphase_out_dir"]
         --prefix, -P         Prefix of the output files [Default = "secphase"]
         --disableMarkerMode, -M         If alignments do not overlap with variants Secphase will not switch to marker mode
         --hifi, -x         hifi preset params (only for marker mode) [-q -c -t10 -d 1e-4 -e 0.1 -b20 -m10 -s40 -p40 -r0 -n -10] (Only one of --hifi or --ont should be enabled)
         --ont, -y        ont preset params (only for marker mode) [-q -c -t20 -d 1e-3 -e 0.1 -b20 -m10 -s20 -p20 -r0 -n -10] (Only one of --hifi or --ont should be enabled) 
         --baq, -q         Calculate BAQ [Disabled by default]
         --gapOpen, -d         Gap prob [Default: 1e-4, (for ONT use 1e-2)]
         --gapExt, -e         Gap extension [Default: 0.1]
         --bandwidth, -b         DP bandwidth [Default: 20]
         --consensus, -c         Use consensus confident blocks [Disabled by default]
         --indelThreshold, -t         Indel size threshold for confident blocks [Default: 10 (for ONT use 20)]
         --initQ, -s         Before calculating BAQ set all base qualities to this number [Default: 40 (for ONT use 20)]
         --minQ, -m         Minimum base quality (or BAQ if -q is set) to be considered as a marker  [Default: 20 (for ONT use 10)]
         --primMarginScore, -p         Minimum margin between the consistency score of primary and secondary alignment to select the secondary alignment [Default: 40]
         --primMarginRandom, -r         Maximum margin between the consistency score of primary and secondary alignment to select one randomly [Default: 0]
         --minScore, -n         Minimum marker score of the selected secondary alignment [Default: -10]
         --minVariantMargin, -g         Minimum margin for creating blocks around phased variants [Default: 50]
         --minGQ, -G         Minimum genotype quality of the phased variants [Default: 10]
```
Some notes:

- By default Secphase will not calculate BAQ for markers. For enabling BAQ the parameteres `--baq` and `--consensus` have to be set. For HiFi reads preset `--hifi` and for ONT `--ont` is recommended. These presets will enable `--baq --consensus` automatically. 
- The input bam file must be sorted **by read name and contain the** `cs` **tag**. 
- The sorted bam file should be indexed with `secphase_index` prior to running `secphase`; it will create an index file with the suffix `.secphase.index` which is necessary for multi-threading. 

In summary given a bam file (usually sorted by reference position) you can run these lines:

```
## Sort by read name
samtools sort -n -@8 ${INPUT_DIR}/${BAM_PREFIX}.bam > ${INPUT_DIR}/${BAM_PREFIX}.sorted_qname.bam

## Create secphase index
docker run \
	-v ${INPUT_DIR}:${INPUT_DIR} \
	mobinasri/secphase:v0.4.2 \
	secphase_index \
	-i ${INPUT_DIR}/${BAM_PREFIX}.sorted.bam

## Run Secphase for HiFi alignments
docker run \
	-v ${INPUT_DIR}:${INPUT_DIR} \
	-v ${OUTPUT_DIR}:${OUTPUT_DIR}
	mobinasri/secphase:v0.4.2 \
	secphase --hifi \
	-i ${INPUT_DIR}/${BAM_PREFIX}.sorted.bam \
	-f ${INPUT_DIR}/${FASTA_PREFIX}.fa \
	--outDir ${OUTPUT_DIR} \
	--prefix ${SAMPLE_NAME} \
	--threads ${THREAD_COUNT}

## Run Secphase for ONT alignments
docker run \
	-v ${INPUT_DIR}:${INPUT_DIR} \
	-v ${OUTPUT_DIR}:${OUTPUT_DIR}
	mobinasri/secphase:v0.4.2 \
	secphase --ont \
	-i ${INPUT_DIR}/${BAM_PREFIX}.sorted.bam \
	-f ${INPUT_DIR}/${FASTA_PREFIX}.fa \
	--outDir ${OUTPUT_DIR} \
	--prefix ${SAMPLE_NAME} \
	--threads ${THREAD_COUNT}
```
`${OUTPUT_DIR}/${SAMPLE_NAME}.out.log` conatins the names of the reads whose secondary and primary alignments have to be swapped.
Here is an example of a record in this file:

```
$	m64043_200710_174426/2353/ccs
*	-35.00	HG00438#1#JAHBCB010000044.1	4904283	4924283
@	-0.00	HG00438#2#JAHBCA010000036.1	3898995	3918995
```

Based on the initial letter of each line we can identify the correct phasing of each read:
- `$` proceeds the read name which is `m64043_200710_174426/1311005/ccs` in this example
- `*` proceeds the score and the start/end position of the primary alignment which is not to the correct haplotype
- `@` proceeds the score and the start/end position of the secondary alignment which is to the correct haplotype
- `!` proceeds the score and the start/end position of any other alignments (This example has only two alignments)

### Running Secphase in dual mode

It is also possible to give a phased vcf file to Secphase along with the desired regions in a bed file. These variants will help phasing the read alignments to a diploid assembly. Using the command below Secphase will use variants (variant mode) for phasing whenever a read spans over at least one variant in at least one of its primary/secondary alignments. Whenever a read has no overlap with variants it will use the original marker mode. This mode is called dual mode since it is switching between variant and marker mode. It is possible to disable the dual mode and use only the variant mode by using `--disableMarkerMode`.

For producing a phased vcf file one approach can be aligning all HiFi reads to each haplotype, call variants (by DeepVariant for example) and phasing them with ONT Ultra Long reads (by Whatshap for example). It is recommended to narrow down this process to homozygous regions. To obtain the homozygous blocks the script `find_homozygous_regions.py` can be employed (described below).

```
## Run Secphase for HiFi alignments in dual mode
docker run \
	-v ${INPUT_DIR}:${INPUT_DIR} \
	-v ${OUTPUT_DIR}:${OUTPUT_DIR}
	mobinasri/secphase:v0.4.2 \
	secphase --ont \
	-i ${INPUT_DIR}/${BAM_PREFIX}.sorted.bam \
	-f ${INPUT_DIR}/${FASTA_PREFIX}.fa \
	--inputVcf ${INPUT_DIR}/phased_variants.vcf.gz \
	--variantBed ${INPUT_DIR}/variant_blocks.bed \
	--outDir ${OUTPUT_DIR} \
	--prefix ${SAMPLE_NAME} \
	--threads ${THREAD_COUNT}
```
### Correcting bam file with secphase output

To swap the pri/sec tags of the reads reported in `*.out.log` and produce a modified bam file you can run the program  `correct_bam`.
Again it is recommended to run it using the docker image `mobinasri/secphase:v0.4.2`.

Here are the parameters `correct_bam` can accept:
```
Usage: correct_bam  -i <INPUT_BAM> -o <OUTPUT_BAM> -p <PHASING_LOG> -m <MAPQ_TABLE>
	Modify the input bam file:
	* Apply the phasing log by swapping the primary and secondary alignments whenever necessary(stdout log of ./phase_reads)
	* Set the MAPQs to the values given in the mapq table
		mapq table is a tab delimited text containing 4 columns:
		1. read name
		2. contig name
		3. left-most coordinate on contig (1-based)
		4. adjusted mapq
	* Filter secondary alignments (After applying the phasing log)
	* Skip outputing the optional fields (like cs and MD tags)
	* Filter the reads shorter than the given threshold
	* Filter the alignments shorter than the given threshold

Options:
         --inputBam,	-i         input bam file
         --outputBam,	-o         output bam file
         --maxMapq,	-x         maximum mapq [default:100]
         --phasingLog,	-P         the phasing log path (output of secphase) [optional]
         --mapqTable,	-M         the adjusted mapq table (tab-delimited) path (4 columns with no header: read_name, contig_name, 1_based_contig_start , new_mapq) [optional]
         --exclude,	-e         Path to a file containing the read names that have to be excluded [optional]
         --noTag,	-t         output no optional fields
         --primaryOnly,	-p         output only primary alignments
         --minReadLen,	-m         min read length [default: 5k]
         --minAlignmentLen,	-a         min alignment length [default: 5k]
         --maxDiv,	-d         min gap-compressed divergence ("de" tag) [default: 0.12]
         --threads,	-n         number of threads (for bam I/O)[default: 2]
```

To produce the modified bam file: (Here the input bam file can be sorted by reference position but must contain the `cs` tag)

```
docker run \
	-v ${INPUT_DIR}:${INPUT_DIR} \
	-v ${OUTPUT_DIR}:${OUTPUT_DIR} \
	mobinasri/secphase:v0.4.2 \
	correct_bam \
	-i ${INPUT_DIR}/${BAM_PREFIX}.bam \
	-P ${OUTPUT_DIR}/${SAMPLE_NAME}.out.log \
	-o ${OUTPUT_DIR}/${BAM_PREFIX}.corrected.bam \
	--primaryOnly
```

Note the default values for `--minReadLen` and `--minAlignmentLen` are both `5k` and should be changed if not desired.

### Detecting homozygous regions  

The script `programs/src/find_homozygous_regions.py` can be used to detect homozygous regions in a diploid assembly. It works by parsing the paf file of an alignment between two haplotypes and returning windows of 100% identity.

You can use minimap2 or winnowmap to get a paf file from the alignment of the two haplotype assemblies, but you must include the `--eqx` flag so that matches/mismatches are included in the cigar string (`X | =`).

Example:
```
minimap2 -t12 --eqx -c --cs -x asm5 HG002.paternal.f1_assembly_v2_genbank.fa HG002.maternal.f1_assembly_v2_genbank.fa > HG002_pat2mat_mm2.paf
```
How to run `find_homozygous_regions.py` :

```
Usage:
python3 find_homozygous_regions.py -p <paf_file> -m <min_length> -e <extend_windows> -o <prefix>

options:
  -h, --help            show this help message and exit
  -p PAF_FILE, --paf_file PAF_FILE
                        paf file of aligments between two haplotype assemblies. --eqx flag must be used in alignment
  -m MIN_LENGTH, --min_length MIN_LENGTH
                        Minimum window size in bp for homozygous regions
  -e EXTEND_WINDOWS, --extend_windows EXTEND_WINDOWS
                        (Optional) Number of bp to extend windows by, to capture surrounding heterozygosity
  -o OUT_BED, --out_bed OUT_BED
                        Output prefix
```

### Workflows

Each of the phasing and correction programs is wdlized separately. The phasing wdl file is named [`secphase.wdl`](https://dockstore.org/workflows/github.com/mobinasri/secphase/Secphase) and the correction wdl file is named [`correct_bam.wdl`](https://dockstore.org/my-workflows/github.com/mobinasri/secphase/CorrectBam). These links can be used for importing the workflows to Terra or running locally with cromwell.
Some notes about the inputs to these workflows:
- For `secphase.wdl` set `secphase.options` to `--ont` for ONT data and change it to `--hifi` for HiFi data.

### Acknowledgements

Thanks to Heng Li for sharing the core idea of [`mmphase`](https://github.com/lh3/minimap2/tree/master/misc) which is incorporated in Secphase with the modifications mentioned above.
