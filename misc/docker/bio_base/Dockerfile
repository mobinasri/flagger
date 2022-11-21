FROM ubuntu:18.04
MAINTAINER Mobin Asri, masri@ucsc.edu

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -y build-essential git zlib1g-dev wget libbz2-dev pkg-config libncurses5-dev libncursesw5-dev liblzma-dev libcurl4-openssl-dev pigz python python3 python3-dev curl libssl-dev

RUN mkdir /home/apps

RUN cd /home/apps && wget https://bootstrap.pypa.io/pip/3.6/get-pip.py && python3 get-pip.py

#install cmake
RUN cd /home/apps && \
    wget https://github.com/Kitware/CMake/releases/download/v3.23.2/cmake-3.23.2.tar.gz && \
    tar -xzvf cmake-3.23.2.tar.gz && \
    cd cmake-3.23.2 && \
    ./bootstrap && \
    make && \
    make install

# install samtools
RUN cd /home/apps && \
	wget https://github.com/samtools/samtools/releases/download/1.15.1/samtools-1.15.1.tar.bz2 && \
	tar -vxjf samtools-1.15.1.tar.bz2 && \
	rm -rf samtools-1.15.1.tar.bz2 && \
	cd samtools-1.15.1 && \
	make
ENV PATH="/home/apps/samtools-1.15.1:$PATH"

#install seqkit
RUN cd /home/apps && \
        mkdir seqkit && \ 
        cd seqkit && \
        wget https://github.com/shenwei356/seqkit/releases/download/v2.2.0/seqkit_linux_amd64.tar.gz && \
        tar -xvzf seqkit_linux_amd64.tar.gz
ENV PATH="/home/apps/seqkit:$PATH"
        

#install seqtk
RUN cd /home/apps && \
        git clone https://github.com/lh3/seqtk.git && \
        cd seqtk && \
        make
ENV PATH="/home/apps/seqtk:$PATH"

# install bedtools
RUN cd home/apps && \
	wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz && \
	tar -zxvf bedtools-2.29.1.tar.gz && \ 
        rm -rf bedtools-2.29.1.tar.gz && \
	cd bedtools2 && \
	make
ENV PATH="/home/apps/bedtools2/bin:$PATH"

#download picard.jar
RUN cd /home/apps && \
    wget https://github.com/broadinstitute/picard/releases/download/2.26.10/picard.jar

#install bcftools
RUN cd /home/apps && \
        wget https://github.com/samtools/bcftools/releases/download/1.12/bcftools-1.12.tar.bz2 && \
        tar -vxjf bcftools-1.12.tar.bz2 && \
        rm -rf bcftools-1.12.tar.bz2 && \
        cd bcftools-1.12 && \
        make
ENV PATH=/home/apps/bcftools-1.12:$PATH

## UCSC convention is to work in /data
RUN mkdir -p /data
WORKDIR /data
