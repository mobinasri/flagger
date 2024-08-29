## Flagger: A mixture model-based coverage analysis for assessing diploid assemblies

### Overview
The main purpose of the **Flagger** pipeline is to assess a diploid assembly. To use this pipeline a BAM file containing the read alignments to the diploid assembly should be prepared in advance. Using the BAM file this pipeline is able to flag mis-assemblies by detecting anomalies in the coverage distribution along the assembly. It can also categorize the mis-assemblies into 3 main groups: erroneous, (falsely) duplicated, and collapsed. 

The examples shown here are from the [Human Pan-Genome Project](https://humanpangenome.org/) since it was the main motivation for developing this pipeline.

The pipeline has 3 core steps:
- Calculate the read coverage of each assembly base 
- Fit a mixture model to the coverage distributions
- Extract the blocks assigned to the model's 4 main components: erroneous, duplicated, haploid, and collapsed.

### Docker
All programs used in this analysis are available in the docker image `mobinasri/flagger:v0.3.2`. It is recommended to use this image for running the programs.

## The Pipeline

### 1. Calculating Depth of Coverage
Given the read alignments in the BAM format it is possible to calculate the the depth of coverage for each base by `samtools depth`. The output of `samtools depth -aa` is like below. (`-aa` option allows outputing the bases with zero coverage)
````
contig_1  1 0
contig_1  2 1
contig_1  3 1
contig_1  4 1
contig_1  5 2
contig_1  6 2
````
In the order in which they appear, the columns are showing the contig name, the base coordinate and the coverage. 
Each base has a separate line even if consecutive bases are having the same coverage. 
In order to make this output more compact it can be converted to the format below
````
>contig_1 6
1 1 0
2 4 1
5 6 2
````
In which each contig's name appears only once before the first block of that contig and the consecutive bases with the same coverage take only one line.
The number that comes after the name of the contig is the contig size. The coverage files with such a format have the suffix `.cov`.
Here are the main commands for producing the coverage files:
````
## Find base-level coverages
samtools depth -aa -Q 0 ${INPUT_DIR}/read_alignment.bam > ${INPUT_DIR}/read_alignment.depth

## Convert depth to cov
docker run \
 -v ${INPUT_DIR}:${INPUT_DIR} \
 -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
 mobinasri/flagger:v0.4.0 \
 depth2cov \
 -d ${INPUT_DIR}/read_alignment.depth \
 -f ${INPUT_DIR}/asm.fa.fai \
 -o ${OUTPUT_DIR}/read_alignment.cov
````
`depth2cov` is a program that converts the output of `samtools depth` to `.cov` format. Its source code is available [here](https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/docker/coverage/scripts/depth2cov.c).

Since the reads aligned to the homozygous regions are expected to have low mapping qualities we don't filter reads based on their mapping qualities.
In the figure below you can see the histograms of mapping qualities and the distributions of alignment indentities for HG00438 as an example. The diploid assembly and HiFi data set used here are from HPRC-Year1 repositories ([assembly repo](https://github.com/human-pangenomics/HPP_Year1_Assemblies) and [read set repo](https://github.com/human-pangenomics/HPP_Year1_Data_Freeze_v1.0)).  Three sets of alignments are shown here; the alignments to the diploid assembly and to each haploid assembly (maternal and paternal) separately. The alignments to the haploid assemblies are shown here just for comparison and are not used for the current analysis.

<img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/docs/coverage/images/HG00438_mapq_hist.png" width="700" height="275">

In this example about 20% of the diploid alignments are having MAPQs lower than 20. To correctly phase the reads with low MAPQ alignments it is recommended to run [Secphase](https://github.com/mobinasri/secphase) before calculating the coverage.

### 2. Coverage Distribution and Fitting The Mixture Model


The frequencies of coverages can be calculated w/ `cov2counts`. Its source code is available in [cov2counts.c](https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/docker/coverage/scripts/cov2counts.c)
````
docker run \
 -v ${INPUT_DIR}:${INPUT_DIR} \
 -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
 mobinasri/flagger:v0.4.0 \
 cov2counts \
 -i ${INPUT_DIR}/read_alignment.cov \
 -o ${OUTPUT_DIR}/read_alignment.counts
````
The output file with `.counts` suffix is a 2-column tab-delimited file; the first column shows coverages and the second column shows the frequencies of
those coverages. Using `.counts` files we can produce distribution plots easily. For example below we are showing the coverage distribution for HG00438 diploid assembly.
<img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/docs/coverage/images/HG00438_dip_hifi_cov_dist.png" width="700" height="275">

The python script `fit_gmm.py` is able to take a file `.counts` suffix and fit a mixture model and find the best parameters through 
[Expectation Maximization (EM)](https://en.wikipedia.org/wiki/Expectation%E2%80%93maximization_algorithm). This mixture model
consists of 4 main components and each component represents a specific type of regions.

The 4 components:

1. **Erroneous component**, which is modeled by a poisson distribution. To avoid overfitting, this mode only uses the coverages below 10 so its mean is limited to be between 0 and 10. It represents the regions with very low read support.
2. **(Falsely) Duplicated component**, which is modeled by a gaussian distribution which mean is constrained to be half of the haploid component's mean. It should mainly represents the falsely duplicated regions. It is worth noting that according to the recent 
[T2T paper, The complete sequence of a human genome,](https://www.biorxiv.org/content/10.1101/2021.05.26.445798v1.abstract) there exist
 some satellite arrays (especially HSAT1) where the ONT and HiFi coverage drops systematically due to bias in sample preparation and sequencing.
 As a result this mode should contain a mix of duplicated and coverage-biased blocks. To avoid overestimating this component a correction step is added to the pipeline.
3. **Haploid component**, which is modeled by a gaussian distribution. It represents blocks with the coverages that we expect for the blocks of an error-free assembly.
4. **Collpased component**, which is actually a set of components each of which follows a gaussian distribution and their means are constrained to be
multiples of the haploid component's mean. Similar to the duplicated component there are some regions where HiFi coverage increases systematically. This issue is addressed in one of the correction steps.

Here is the command that fits the model:
````
## Run fit_model_extra.py to fit the model
docker run \
 -v ${INPUT_DIR}:${INPUT_DIR} \
 -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
 mobinasri/flagger:v0.4.0 \
 python3 /home/programs/src/fit_gmm.py \
 --counts ${INPUT_DIR}/read_alignment.counts \
 --cov ${EXPECTED_COVERAGE} \
 --output ${OUTPUT_DIR}/read_alignment.table \
````

Its output is a file with `.table` suffix. It contains a TAB-delimited table with the 7 fields described below:
|Col|Type  |Description                               |
|:--|:----:|:-----------------------------------------|
|coverage  |int|A depth of coverage                       |
|freq  |float   |The frequency of the coverage value in the first column                     |
|fit  |float   |The frequency value fit to the model    |
|error  |float   |The weight of the error component       |
|duplicated  |float  |The weight of the duplicated component               |
|haploid  |float|The weight of the haploid component                      |
|collapsed  |float   |The weight of the collapsed component                    |

Here is an example of such a table:
````
#coverage	freq	fit	error	duplicated	haploid	collapsed
0	1633030	2479818.7542	0.9973	0.0009	0.0019	0.0000
1	552306	1014768.3370	0.9890	0.0044	0.0067	0.0000
2	311274	425804.6727	0.9564	0.0207	0.0229	0.0000
3	241294	195958.8474	0.8434	0.0856	0.0710	0.0000
4	201587	117491.9394	0.5708	0.2616	0.1676	0.0000
5	207799	108998.0736	0.2497	0.4973	0.2530	0.0000
...
...
...
64	1211029	680304.0405	0.0000	0.0000	0.4493	0.5507
65	1111166	635906.0944	0.0000	0.0000	0.3724	0.6276
66	1034004	604733.6596	0.0000	0.0000	0.3004	0.6996
67	948793	584283.3204	0.0000	0.0000	0.2362	0.7638
68	849748	572280.0759	0.0000	0.0000	0.1814	0.8186
69	800708	566695.7729	0.0000	0.0000	0.1364	0.8636
70	766435	565754.1385	0.0000	0.0000	0.1008	0.8992
````

The figure below shows the model components. It is again for the HG00438 diploid assembly:

<img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/docs/coverage/images/HG00438_fit_model.png" width="700" height="275">

### 3. Extracting Blocks 

Now we have to assign each coverage value to one of the four components (erroneous, duplicated, haploid, and collapsed). To do so, for each coverage value,
we pick the component with the highest probability. For example for the coverage value, 0, is assigned to the erroneous component most of the times (the red line).
In the figure below the coverage intervals are colored based on their assigned component.

<img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/docs/coverage/images/HG00438_fit_model_colored.png" width="700" height="275">

After assigning the coverage values we then assign the bases of the corresponding assembly to the most probable component. Finally we will have 4 bed 
files each of which points to the regions assinged to a single component.
`find_blocks` program is able to receive a coverage file (suffix `.cov`) and its corresponding table file (suffix `.table`) and produce the
4 mentioned bed files.
````
## Run find_blocks_from_table
docker run \
 -v ${INPUT_DIR}:${INPUT_DIR} \
 -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
 mobinasri/flagger:v0.4.0 \
 find_blocks_from_table \
 -c ${INPUT_DIR}/read_alignment.cov \
 -t ${INPUT_DIR}/read_alignment.table \
 -p ${OUTPUT_DIR}/${OUTPUT_PREFIX}
````

It outputs four bed files:
````
${prefix}.error.bed
${prefix}.duplicated.bed
${prefix}.haploid.bed
${prefix}.collapsed.bed
````

## Additional Correction Steps
It has been noticed that the results of the previous steps need corrections. In the neccessary steps below 3 main issues are explained and 
their solutions are also provided. The output of each correction step is used as an input for the next one.

### 1. Incorporating HSATs Coverage Bias

As it was mentioned in the 2nd step of the pipeline, there are some HSats in the genome where the HiFi or ONT coverage is systematically increased or decreased.
Such platform-specific biases mislead the pipeline. For example the HiFi coverage of Hsat2 in chr1 is about 1.5 times larger than the average sequencing coverage and it is very probable to wrongly flag this region as collapsed.

In the figure below the read bars are showing the coverage of the HiFi alignments to the chm13v1.0.

   <img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/docs/coverage/images/chm13v1.1_hifi_hsat2_chr1.png" width="700" height="150">
To incorporate such coverage biases and correct the results in the corresponding regions, the steps below are neccessary:

 1. Align the diploid assembly to an (HSat-)annotated reference like chm13v1.1
 2. Find the assembly blocks aligned to the HSat of interest. To do so we need a bed file pointing to the HSat in the reference then we can run the script `project_blocks_multi_thread.py` to project it back to the assembly.
 3. Make a separate coverage file that covers the HSat blocks
 4. Run the pipeline again for the HSat-specific coverage file and obtain the 4 bed files
 5. Use this partitioning instead of the one obtained by whole-genome model
 
````
## ${ASM2REF}.bam is a bam file containing the assembly-to-ref alignments
## First convert bam to paf since project_blocks_multi_thread.py works with paf
k8 paftools.js sam2paf <(samtools view -h -F4 -F256 ) ${ASM2REF}.bam > ${ASM2REF}.paf

## Given the paf file extract the HSat regions in the assembly 
## and save them in ${HSAT_PROJECTION}.bed
docker run \
 -v ${INPUT_DIR}:${INPUT_DIR} \
 -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
 mobinasri/flagger:v0.4.0 \
 python3 /home/programs/src/project_blocks_multi_thread.py 
 --threads 8
 --mode 'ref2asm' \
 --paf ${INPUT_DIR}/${ASM2REF}.paf \
 --blocks ${HSAT}.bed \
 --outputProjectable ${OUTPUT_DIR}/${PROJECTABLE}.bed \
 --outputProjection ${OUTPUT_DIR}/${HSAT_PROJECTION}.bed

## Sort and merge ${PROJECTION}.bed  
bedtools sort -i  ${OUTPUT_DIR}/${HSAT_PROJECTION}.bed | bedtools merge -i - > ${OUTPUT_DIR}/${HSAT_PROJECTION}.merged.bed

## Given the HSat blocks in the assembly, ${HSAT_PROJECTION}.merged.bed
## and the whole genome coverage file, ${COVERAGE}.cov.gz
## make a smaller coverage file covering only the HSat blocks
zcat ${COVERAGE}.cov.gz | \
            awk '{if(substr($1,1,1) == ">") {contig=substr($1,2,40); len_contig=$2} else {print contig"\t"$1-1"\t"$2"\t"$3"\t"len_contig}}' | \
            bedtools intersect -a - -b ${HSAT_PROJECTION}.merged.bed | \
            awk '{if(contig != $1){contig=$1; print ">"contig" "$5}; print $2+1"\t"$3"\t"$4}' | pigz -p4 > ${HSAT_COVERAGE}.cov.gz
````
`${HSAT_COVERAGE}.cov.gz` should be used as for the 2nd and 3rd steps of the pipeline. Please note that while calling `fit_gmm.py` the parameter
`--coverage` (The starting point of the fitting process) is recommended to be adjusted based on the expected coverage in the corresponding HSat.

The process should be repeated for each HSat of interest that can be prone to covarage bias. For example for HiFi data it is recommended to do it for each of HSat1, HSat2 and Hsat3 and adjust `--coverage` as a function of the average sequencing coverage `${AVG_COVERAGE}`.

 - HSat1 -> `--coverage  0.75 * ${AVG_COVERAGE}`
 - HSat2 -> `--coverage  1.25 * ${AVG_COVERAGE}`
 - HSat3 -> `--coverage  1.25 * ${AVG_COVERAGE}`

Please note that there may exist no bias for a specific sample and fitting a separate model means adding more parameters in the big picture not assuming the HSat coverage to be always greater (or less) than the whole-genome coverage. The factors provided above are just the initial values to navigate the EM algorithm if there exist a bias. The EM algorithm should be able to find the correct parameter values regardless of having a bias or not.
The HSat corrected bed files contain `hsat_corrected` in the their file names. (Among the outputs of flagger.wdl)

### 2. Window-Specific Models

In step 3 a single model is fit for the whole diploid assembly. In step 4 that model is used to partition the assembly into 4 main components. It is noticed that the model components may change for different regions and it may affect the accuracy of the partitioning process. In order to make the coverage thresholds more sensitive to the local patterns the diploid assembly is split into windows of length (5-10Mb). For each window a separate model  should be fit. To do so first we split the whole-genome coverage file produced in step 3 into multiple coverage files one for each window.
```

## First we make a cov file excluding the regions that have may have coverage biases (Read the following correction step)
cat asm.fa.fai | awk '{print $1"\t0\t"$2}' | sort -k1,1V -k2,2n > asm.bed
## exclude.bed contains all regions that have may have coverage biases
bedtools subtract -a asm.bed -b exclude.bed > asm.no_bias.bed

## Make a tab-delimited file including the contig names and their effective length (not excluded)
cat asm.excluded.bed | awk '{ctg_len[$1] += $3-$2}END{for (c in ctg_len){print c"\t"ctg_len[c]}}' > ctg_lens.txt

## Remove excluded regions from the cov file
cat read_alignment.cov | \
    awk '{if(substr($1,1,1) == ">") {contig=substr($1,2); len_contig=$2} else {print contig"\t"$1-1"\t"$2"\t"$3"\t"len_contig}}' | \
    bedtools intersect -a - -b asm.no_bias.bed | \
    awk '{if(contig != $1){contig=$1; print ">"contig" "$5}; print $2+1"\t"$3"\t"$4}' > read_alignment.no_bias.cov

## Split cov into multiple covs (each covering 5-10 Mb)
docker run \
 -v ${INPUT_DIR}:${INPUT_DIR} \
 -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
 mobinasri/flagger:v0.4.0 \
 split_cov_by_window \
 -c ${INPUT_DIR}/read_alignment.no_bias.cov \
 -f ${INPUT_DIR}/ctg_lens.txt \
 -s 5000000 \
 -p ${OUTPUT_DIR}/${OUTPUT_PREFIX}
```
`split_cov_by_window` program will produce a list of coverage files like below.
```
HG00438.diploid.hifi.HG00438#1#JAHBCB010000001.1_10000001_15000000.cov
HG00438.diploid.hifi.HG00438#1#JAHBCB010000001.1_15000001_20000000.cov
HG00438.diploid.hifi.HG00438#1#JAHBCB010000001.1_1_5000000.cov
HG00438.diploid.hifi.HG00438#1#JAHBCB010000001.1_20000001_25000000.cov
...
HG00438.diploid.hifi.HG00438#2#JAHBCA010000255.1_1_31620.cov
HG00438.diploid.hifi.HG00438#2#JAHBCA010000256.1_1_29665.cov
HG00438.diploid.hifi.HG00438#2#JAHBCA010000257.1_1_30877.cov
HG00438.diploid.hifi.HG00438#2#JAHBCA010000258.1#MT_1_16569.cov
```
To show how the window-specific models may be different from the whole-genome model, here is an example of a false duplication in the maternal assembly of HG01175 that couldn't be detected using the whole-genome model.
The coverage distributions of three contigs along with the whole assembly are shown below. Two of those contigs are maternal and the other one is paternal. They are containing equivalent regions in the genome but one of the maternal contigs (`2#JAHBCE010000091.1` the red one ) is having a few mega bases of a false duplication of the paternal contig (`1#JAHBCF010000027.1` the yellow one).

<img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/docs/coverage/images/HG00438_dip_hifi_contig_cov_dist.png" width="700" height="275">

After fitting the model the duplicated components can reveal such false duplications as you can see in the figure below. `2#JAHBCE010000105.1` does not have any duplicated component but the two others do. Since we know the maternal haplotype has a correct copy of this region in `2#JAHBCE010000105.1` we can deduce that `1#JAHBCF010000027.1` is assembled correctly and `2#JAHBCE010000091.1` is actually containing the false duplication. Without having more information it is not always easy to conclude which copy is the real one and which copy is the false one especially for segmental duplications within a single haplotype.

<img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/docs/coverage/images/HG00438_contig_fit_model.png" width="700" height="400">

One important observation is that for short contigs we don't have a smooth coverage distribution and it is not possible to fit the mixture model. To address this issue we have done the window-specific coverage analysis only for the contigs longer than 5Mb and for the shorter contigs we use the results of the whole-genome analysis described previously. For HSat regions we use the result of the previous correction step.

The combined bed files contain `combined` in their file names. (Among the outputs of flagger.wdl)

### 3. Correcting The Bed Files Pointing To The False Duplications
In some cases the duplicated component is mixed up with the haploid one. It usually happens when the coverage in the haploid component drops systematically and the contig has long stretches of false duplication. 

One other indicator of a false duplication is the accumulation of alignments with very low MAPQ. We have produced another coverage file only for MAPQ > 20 alignments and intersect it with the `combined` results to correct them. Whenever we see a region flagged as duplicated that has more than 5 high quality alignments we change the flag to haploid. One extreme case happened in `HG00438#2#JAHBCA010000050.1`. More than half of this contig is falsely duplicated but the `combined` results (see the previous correction step) reported almost the whole contig as `duplicated`. In the figure below it is shown how this correction step could
fix this issue.

   <img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/docs/coverage/images/dup_correction.png" width="700" height="325">

````
## Get blocks with more than 5 high quality alignments
cat high_mapq.cov | \
            awk '{if(substr($1,1,1) == ">") {contig=substr($1,2,40)} else if($3 >= 5) {print contig"\t"$1-1"\t"$2}}' | \
            bedtools merge -i - > high_mapq.bed

## Do the correction
bedtools subtract -a ${PREFIX}.combined.duplicated.bed -b high_mapq.bed > ${PREFIX}.dup_corrected.duplicated.bed
bedtools intersect -a ${PREFIX}.combined.duplicated.bed -b high_mapq.bed > dup_to_hap.bed
cat dup_to_hap.bed ${PREFIX}.combined.haploid.bed | bedtools sort -i - | bedtools merge -i - > ${PREFIX}.dup_corrected.haploid.bed
````
The duplication corrected bed files contain `dup_corrected` in the file name. (Among the outputs of flagger.wdl)

In the final bed files there is a noticeable number of very short blocks (a few bases or a few tens of bases long). They are very prone to be falsely categorized into one of the components. To increase the specificity we have merged the blocks closer than 100 and then removed the ones shorter than 1Kb. The filtered bed files contain `filtered` in the file name. (Among the outputs of flagger.wdl)

### Notes
1. Some regions are falsely flagged as collapsed. The reason is that the equivalent region in the other haplotype is not assembled correctly so the reads from two haplotypes are aligned to only one of them. This flagging can be useful since it points to a region whose counterpart in the other haplotype is not assembled correctly or not assembled at all. 

    One approach is to correct the coverage in the correctly assembled region. The coverage can be corrected by detecting the marker snps and removing the reads from the wrong halpotype or segment. `filter_alt_reads` can be used for this aim (Read [this](https://github.com/mobinasri/flagger#4-remove-the-alignments-with-alternative-alleles-optional)). Here is an example of a region with ~40X coverage but after detecting the marker snps (by variant calling) and removing the wrong alignments the coverage has decreased to ~17X which is much closer to the expected coverage (~20X).

   <img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/docs/coverage/images/coverage_correction.png" width="700" height="275">


### Data, Source Code and Workflows Availability

The haplotype-resolved assemblies of the HPRC-Y1 samples and their corresponding data sets are available in

https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=working/HPRC/

For more details read this github page:

https://github.com/human-pangenomics/HPP_Year1_Assemblies

We have used the Genbank version of the HPRC-Y1 assemblies.

The Python scripts, C source codes are available in

https://github.com/mobinasri/flagger/tree/main/programs

The wdl files that have been used for this analysis are available in

https://github.com/mobinasri/flagger/tree/main/wdls/workflows

### Results Availability

The results are available in 

https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/e9ad8022-1b30-11ec-ab04-0a13c5208311--COVERAGE_ANALYSIS_Y1_GENBANK/FLAGGER/

### Acknowledgements

Thanks to Jordan M Eizenga for sharing the source code of EM that fits the mixture model and his useful support for doing this analysis.


