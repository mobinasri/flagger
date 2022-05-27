FROM ubuntu:18.04
MAINTAINER Mobin Asri, masri@ucsc.edu

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -y build-essential git zlib1g-dev wget libbz2-dev pkg-config libncurses5-dev libncursesw5-dev liblzma-dev libcurl4-openssl-dev pigz python python3 python3-dev curl

RUN mkdir /home/apps

RUN cd /home/apps && wget https://bootstrap.pypa.io/pip/3.6/get-pip.py && python3 get-pip.py


RUN cd /home/apps && \
	wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2 && \
	tar -vxjf samtools-1.12.tar.bz2 && \
	rm -rf samtools-1.12.tar.bz2 && \
	cd samtools-1.12 && \
	make
ENV PATH="/home/apps/samtools-1.12:$PATH"

RUN cd home/apps && \
	wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz && \
	tar -zxvf bedtools-2.29.1.tar.gz && \ 
        rm -rf bedtools-2.29.1.tar.gz && \
	cd bedtools2 && \
	make
ENV PATH="/home/apps/bedtools2/bin:$PATH"

## UCSC convention is to work in /data
RUN mkdir -p /data
WORKDIR /data
