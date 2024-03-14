FROM ubuntu:18.04
MAINTAINER Julian Lucas, juklucas@ucsc.edu

RUN apt-get update \
    && DEBIAN_FRONTEND="noninteractive" \
    apt-get install -y \
    build-essential git zlib1g-dev wget libbz2-dev pkg-config libncurses5-dev \
    libncursesw5-dev liblzma-dev libcurl4-openssl-dev python3-dev pigz zip

RUN mkdir -p /home/apps

RUN cd /home/apps && \
    wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 && \
    tar -vxjf samtools-1.10.tar.bz2 && \
    rm -rf samtools-1.10.tar.bz2 && \
    cd samtools-1.10 && \
    make
ENV PATH="/home/apps/samtools-1.10:$PATH"

RUN cd /home/apps && \
    git clone https://github.com/lh3/yak && \
    cd yak && \
    git reset --hard f37704a97b68ad08e9129dd7e1ea60072107bda4 && \
    make
ENV PATH="/home/apps/yak:$PATH"

# seqtk
RUN cd /home/apps && \
    git clone https://github.com/lh3/seqtk.git && \
    cd seqtk && \
    git reset --hard 367728dff0f653df87aab51b0911dc1adeb17f6d && \
    make
ENV PATH="/home/apps/seqtk:${PATH}"

# Hifiasm-0.19.8
RUN cd /home/apps && \
    wget https://github.com/chhylp123/hifiasm/archive/refs/tags/0.19.8.zip && \
    unzip 0.19.8.zip  && \
        cd hifiasm-0.19.8 && \ 
    make
ENV PATH="/home/apps/hifiasm-0.19.8:$PATH"

RUN cd /home/apps && \
    git clone --branch v0.4 https://github.com/lh3/gfatools && \
    cd gfatools && \
    git reset --hard 9d4f360538717bad2ea53dd2eef389072430e913 && \
    make
ENV PATH="/home/apps/gfatools:$PATH"

COPY ./scripts /home/scripts

