FROM quay.io/masri2019/qc-stats:latest
MAINTAINER Mobin Asri masri@ucsc.edu


COPY ./scripts  /home/scripts

ENV QC_STATS_GENERATOR_PATH="/home/scripts/haploidQCStatsGenerator.py"
ENV DIP_QC_STATS_GENERATOR_PATH="/home/scripts/diploidQCStatsGenerator.py"

WORKDIR /data
