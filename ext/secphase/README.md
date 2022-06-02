## SecPhase: Phasing the long reads aligned to the diploid assemblies

### Motivation
When we align long reads (HiFi or ONT at this time!) to the diploid assembly of the same sample, the reads coming from homozygous regions may
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
already aligned to. This realignment is performed through a banded HMM which parameters have to be tuned before hand. 
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

### How To Run The Phasing Program

To run the phasing program it is recommended to use the docker image `quay.io/masri2019/hpp_coverage:latest`.

Here are the parameters `phase_reads` can accept:
```
./phase_reads -h

Usage: phase_reads  -i <INPUT_BAM> -f <FASTA> 
Options:
         -i         Input bam file (sorted by qname and must contain cs tag) [required]
         -f         Input fasta file [required]
         -q         Calculate BAQ [Default: false]
         -d         Gap prob [Default: 1e-4, (for ONT use 1e-2)]
         -e         Gap extension [Default: 0.1]
         -b         DP bandwidth [Default: 20]
         -c         Use consensus confident blocks [Default: false]
         -t         Indel size threshold for confident blocks [Default: 10 (for ONT use 20)]
         -s         Before calculating BAQ set all base qualities to this number [Default: 40 (for ONT use 20)]
         -m         Minimum base quality (or BAQ if -q is set) to be considered as a marker  [Default: 20 (for ONT use 10)]
```

The default values are set for the HiFi reads and the input bam file must be sorted by read name and contain the `cs` tag. So given a bam file (usually sorted by reference position) you can run these lines:

```
## Sort by read name
samtools sort -n -@8 ${INPUT_DIR}/${BAM_PREFIX}.bam > ${INPUT_DIR}/${BAM_PREFIX}.sorted_qname.bam

## Phase reads for HiFi
docker run \
	-v ${INPUT_DIR}:${INPUT_DIR} \
	quay.io/masri2019/hpp_coverage:latest \
	phase_reads -q -c -t10 -d 1e-4 -e 0.1 -b20 -m20 -s40 \
	-i ${INPUT_DIR}/${BAM_PREFIX}.sorted.bam \
	-f ${INPUT_DIR}/${FASTA_PREFIX}.fa > ${PHASING_OUT}.log

## Phase reads for ONT
docker run \
	-v ${INPUT_DIR}:${INPUT_DIR} \
	quay.io/masri2019/hpp_coverage:latest \
	phase_reads -q -c -t20 -d 1e-2 -e 0.1 -b20 -m10 -s20 \
	-i ${INPUT_DIR}/${BAM_PREFIX}.sorted.bam \
	-f ${INPUT_DIR}/${FASTA_PREFIX}.fa > ${OUTPUT_DIR}/${PHASING_OUT}.log
```

`${PHASING_OUT}.log` conatins the names of the reads which secondary and primary alignments have to be swapped. 
Here is an example of a record in the `${PHASING_OUT}.log`:

```
$	m64043_200710_174426/2353/ccs
*	-35.00	HG00438#1#JAHBCB010000044.1	4904283
@	-0.00	HG00438#2#JAHBCA010000036.1	3898995
```

Based on the initial letter of each line we can indentify the correct phasing of each read:
- `$` proceeds the read name which is `m64043_200710_174426/1311005/ccs` in this example
- `*` proceeds the score and the start position of the primary alignment which is not to the correct haplotype
- `@` proceeds the score and the start position of the secondary alignment which is to the correct haplotype
- `!` proceeds the score and the start position of any other alignments (This example has only two alignments)

The source code of the phasing program is available in [`phase_reads.c`](https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/docker/coverage/programs/src/phase_reads.c).

### How To Run The Correction Program

To swap the pri/sec tags of the reads reported in `${PHASING_OUT}.log` and produce a modified bam file you can run the program  `correct_bam`.
Again it is recommended to run it using the docker image `quay.io/masri2019/hpp_coverage:latest`.

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
         --inputBam,	-i         input bam file (must contain cs tag)
         --outputBam,	-o         output bam file
         --phasingLog,	-P         the phasing log path [optional]
         --mapqTable,	-M         the adjusted mapq table path [optional]
         --noTag,	-t         output no optional fields
         --primaryOnly,	-p         output only primary alignments
         --minReadLen,	-m         min read length [default: 5k]
         --minAlignmentLen,	-a        min alignment length [default: 5k]
```

To produce the modified bam file: (Here the input bam file can be sorted by reference position but must contain the `cs` tag)

```
docker run \
	-v ${INPUT_DIR}:${INPUT_DIR} \
	-v ${OUTPUT_DIR}:${OUTPUT_DIR} \
	quay.io/masri2019/hpp_coverage:latest \
	correct_bam \
	-i ${INPUT_DIR}/${BAM_PREFIX}.bam \
	-P ${INPUT_DIR}/${PHASING_OUT}.log \
	-o ${OUTPUT_DIR}/${BAM_PREFIX}.corrected.bam \
	--primaryOnly
```

Note the default values for `--minReadLen` and `--minAlignmentLen` are both `5k` and should be changed if not desired.

The source code of the correction program is written in C and is available in [`correct_bam.c`](https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/docker/coverage/programs/src/correct_bam.c).

### Workflows

Each of the phasing and correction programs is wdlized separately and the wdl files can be found  [here](https://github.com/human-pangenomics/hpp_production_workflows/tree/asset/coverage/wdl/tasks). The phasing wdl file is named [`phase_reads.wdl`](https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/wdl/tasks/phase_reads.wdl) and the correction wdl file is named [`correct_bam.wdl`](https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/wdl/tasks/correct_bam.wdl).

### Acknowledgements

Thanks to Heng Li for sharing the core idea of [`mmphase`](https://github.com/lh3/minimap2/tree/master/misc) which is incorporated in `phase_reads` with some modifications.
