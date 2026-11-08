## Running HMM-Flagger on a small test (Linux and macOS)

The test data covers the first 10 Mb of HG002 chr15 in both haplotypes, with 40x PacBio HiFi
alignments (winnowmap). It runs in a few seconds with the prebuilt binaries.

### 1. Download the binaries

```
VERSION=v1.2.1    # pick a release

# Linux (amd64 or arm64)
ARCH=$(uname -m | sed 's/x86_64/amd64/; s/aarch64/arm64/')
PLATFORM=linux-${ARCH}

# macOS (Apple Silicon, macOS 11 or newer)
PLATFORM=macos-arm64

curl -fLO https://github.com/mobinasri/flagger/releases/download/${VERSION}/hmm-flagger-${VERSION}-${PLATFORM}.zip
unzip hmm-flagger-${VERSION}-${PLATFORM}.zip
export PATH=${PWD}/hmm-flagger-${VERSION}-${PLATFORM}:${PATH}
```

Set only the `PLATFORM` line for your system. On macOS, download with `curl` as above. A zip
downloaded with a browser is blocked until you run
`xattr -dr com.apple.quarantine hmm-flagger-${VERSION}-macos-arm64`.

### 2. Download the test data

```
curl -fLO https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/e093fd72-e31a-11ee-b020-27964ee37032--flagger_test_files/flagger_v0.4.0/test_files/test_flagger_end_to_end/small_test.tar.gz
tar -xzf small_test.tar.gz
cd small_test
```

### 3. Run bam2cov and hmm_flagger

```
# whole-genome BED and the annotation JSON for bam2cov
awk '{print $1"\t0\t"$2}' fasta_files/*.dip.chr15_0_to_10Mb.fa.fai > whole_genome.bed
echo "{\"whole_genome\" : \"${PWD}/whole_genome.bed\"}" > annotations.json

bam2cov --bam bam_files/*.bam \
        --output coverage.cov.gz \
        --annotationJson annotations.json \
        --baselineAnnotation whole_genome \
        --threads 8

# alpha matrix for the HiFi preset
curl -fLO https://raw.githubusercontent.com/mobinasri/flagger/refs/heads/main/misc/alpha_tsv/HiFi_DC_1.2/alpha_optimum_trunc_exp_gaussian_w_16000_n_50.HiFi_DC_1.2_DEC_2024.v1.1.0.tsv

mkdir -p hmm_flagger_outputs
hmm_flagger --input coverage.cov.gz \
            --outputDir hmm_flagger_outputs \
            --alphaTsv alpha_optimum_trunc_exp_gaussian_w_16000_n_50.HiFi_DC_1.2_DEC_2024.v1.1.0.tsv \
            --windowLen 16000 \
            --labelNames Err,Dup,Hap,Col \
            --threads 8
```

### 4. Check the result

The predictions are in `hmm_flagger_outputs/final_flagger_prediction.bed`. The summary should show
about 99.8% of the bases as `Hap` (columns: Err, Dup, Hap, Col, Unk):

```
grep base_level hmm_flagger_outputs/prediction_summary_final.tsv | grep percentage | grep ALL_SIZES | grep whole_genome
PREDICTION	base_level	percentage	annotation	whole_genome	ALL_SIZES	ALL_LABELS	0.00	0.00	99.76	0.24	0.00
```
