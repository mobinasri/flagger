FROM mobinasri/bio_base:v0.3.0
MAINTAINER Mobin Asri, masri@ucsc.edu

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -y autoconf bzip2 lzma-dev zlib1g-dev libcurl4-openssl-dev libcrypto++-dev libpthread-stubs0-dev libbz2-dev liblzma-dev

RUN cd /home/apps && \
    git clone https://github.com/benedictpaten/marginPhase && \
    cd marginPhase && \
    git submodule update --init && \
    mkdir build && \ 
    cd build && \
    cmake .. && \
    make 
ENV PATH=/home/apps/marginPhase/build:$PATH


## UCSC convention is to work in /data
RUN mkdir -p /data
WORKDIR /data
