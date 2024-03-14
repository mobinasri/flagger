FROM biocontainers/biocontainers:v1.0.0_cv5

MAINTAINER Julian Lucas, juklucas@ucsc.edu

LABEL base_image="biocontainers:v1.0.0_cv5"
LABEL version="0.1"
LABEL software="primrose"
LABEL software.version="1.3.0"
LABEL about.summary="Primrose is a tool to predict 5mC in CpG motifs of hifi_reads.bam files (that have kinetics data). Outputs another bam with 5mC predictions."
LABEL about.home="https://github.com/PacificBiosciences/primrose"
LABEL about.documentation="https://github.com/PacificBiosciences/primrose"
LABEL about.license_file="https://github.com/PacificBiosciences/primrose/blob/main/LICENSE"


USER root

RUN apt-get clean all && \
    apt-get update && \
    conda install -c bioconda primrose && \
    apt-get update && \
    apt-get clean all && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

USER biodocker

WORKDIR /data/