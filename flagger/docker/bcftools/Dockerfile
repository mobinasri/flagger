FROM quay.io/masri2019/hpp_base:latest
MAINTAINER Mobin Asri, masri@ucsc.edu

RUN pip3 install numpy
RUN mkdir -p /home/apps

## bcftools
RUN cd /home/apps && \
        wget https://github.com/samtools/bcftools/releases/download/1.12/bcftools-1.12.tar.bz2 && \
        tar -vxjf bcftools-1.12.tar.bz2 && \
        rm -rf samtools-1.12.tar.bz2 && \
        cd bcftools-1.12 && \
        make
ENV PATH="/home/apps/bcftools-1.12:$PATH"

COPY ./scripts /home/scripts
ENV SPLIT_BED_PY=/home/scripts/split_bed.py
ENV SPLIT_BED_CONTIG_WISE_PY=/home/scripts/split_bed_contig_wise.py
ENV FILTER_ALT_READS_BIN=/home/scripts/bin/filter_alt_reads

## UCSC convention is to work in /data
RUN mkdir -p /data
WORKDIR /data
