FROM quay.io/masri2019/hpp_base:latest
MAINTAINER Mobin Asri, masri@ucsc.edu

RUN mkdir -p /home/apps

## cutadapt

RUN python3 -m pip install --upgrade pip && pip install git+https://github.com/marcelm/cutadapt.git

## blast

RUN cd /home/apps && \
    wget https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.11.0+-x64-linux.tar.gz && \
    tar -xzvf ncbi-blast-2.11.0+-x64-linux.tar.gz && \
    rm -rf ncbi-blast-2.11.0+-x64-linux.tar.gz
ENV PATH="/home/apps/ncbi-blast-2.11.0+/bin:${PATH}"

## seqtk

RUN cd /home/apps && \
    git clone https://github.com/lh3/seqtk.git && \
    cd seqtk && \
    make
ENV PATH="/home/apps/seqtk:${PATH}"

## HiFiAdapterFilt

COPY ./HiFiAdapterFilt /home/apps/HiFiAdapterFilt
# pbadapterfilt_fastq.sh is the modified version of pbadapterfilt.sh that recieves fastq as input
ENV HIFI_ADAPTER_FILTER_BASH=/home/apps/HiFiAdapterFilt/pbadapterfilt_fastq.sh
ENV PATH="/home/apps/HiFiAdapterFilt/DB:/home/apps/HiFiAdapterFilt:${PATH}"

## UCSC convention is to work in /data
RUN mkdir -p /data
WORKDIR /data
