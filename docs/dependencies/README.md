## Dependencies

This page lists the main dependencies of HMM-Flagger v1.2.1 and the versions installed in its Docker images. Each image is published for linux/amd64 and linux/arm64 under a fixed version tag, so a given tag always contains the same versions. The WDL workflows use these tags by default.

| Image | Dockerfile |
|---|---|
| `mobinasri/bio_base:v1.2.1` | [misc/docker/bio_base/Dockerfile](../../misc/docker/bio_base/Dockerfile) |
| `mobinasri/flagger:v1.2.1` (built on `bio_base`) | [Dockerfile](../../Dockerfile) |
| `mobinasri/long_read_aligner:v1.2.1` (built on `bio_base`) | [misc/docker/long_read_aligner/Dockerfile](../../misc/docker/long_read_aligner/Dockerfile) |

Base system: Ubuntu 20.04.6 LTS.

### C libraries (used to build HMM-Flagger and bam2cov)

| Dependency | Version |
|---|---|
| [htslib](https://github.com/samtools/htslib) | 1.13 |
| [sonLib](https://github.com/benedictpaten/sonLib) | commit `5cbc158` (2022-09-14) |

### Command-line tools

| Tool | Version | Image |
|---|---|---|
| [samtools](https://github.com/samtools/samtools) | 1.23.1 | bio_base |
| [bcftools](https://github.com/samtools/bcftools) | 1.23.1 | bio_base |
| [bedtools](https://github.com/arq5x/bedtools2) | 2.31.1 | bio_base |
| [seqkit](https://github.com/shenwei356/seqkit) | 2.13.0 | bio_base |
| [seqtk](https://github.com/lh3/seqtk) | 1.5-r133 (commit `94e7070`) | bio_base |
| [picard](https://github.com/broadinstitute/picard) | 3.4.0 | bio_base |
| bedGraphToBigWig (UCSC) | v369 (amd64), 482 (arm64) | flagger |

### Aligners

| Tool | Version | Image |
|---|---|---|
| [minimap2](https://github.com/lh3/minimap2) | 2.31 (r1302) | long_read_aligner |
| [Winnowmap](https://github.com/marbl/Winnowmap) | 2.03 | long_read_aligner |
| [k8](https://github.com/attractivechaos/k8) | 0.2.4 (amd64), 1.2 (arm64) | long_read_aligner |
| [VerityMap](https://github.com/ablab/VerityMap) | 2.1.2-alpha | long_read_aligner |

### Python packages (Python 3.8.10, `flagger` image)

| Package | Version |
|---|---|
| numpy | 1.24.4 |
| scipy | 1.10.1 |
| pandas | 2.0.3 |
| matplotlib | 3.7.5 |
| biopython | 1.83 |
| pysam | 0.24.0 |
| whatshap | 2.5 |
| smt | 2.6.3 |
| Cython | 3.2.9 |

### Building from source

To build the binaries without Docker, see [Building from source](../build_from_source/README.md).
