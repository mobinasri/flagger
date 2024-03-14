FROM ubuntu:20.04
MAINTAINER Julian Lucas, juklucas@ucsc.edu

LABEL    software="verkko_hic" \ 
    container="verkko-hic" \ 
    about.summary="HiC Phasing For Verkko Assembly Runs" \ 
    about.home="https://github.com/Dmitry-Antipov/verkkohic"


RUN apt-get update && apt-get -y upgrade \
  && apt-get install -y --no-install-recommends \
    git \
    wget \
    g++ \
    pigz \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*


## Install Verkko, Mashmap, & Networkx w/ Miniconda
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
ADD environment.yml /tmp/environment.yml

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh \
    && conda init bash \
    && . /root/.bashrc \
    && conda update conda \
    && conda env create -f /tmp/environment.yml


## verkko hic
WORKDIR /opt
RUN git clone https://github.com/Dmitry-Antipov/verkkohic.git && \
    cd verkkohic && \
    git reset --hard b8fdda6aac3245f2ee39e3c0d6c01b58debcc092

ENV PATH="/opt/verkkohic:${PATH}"

## Verkko HiC requires the following parameters to be set:
ENV VERKKO="/root/miniconda3/envs/verkko_hic"


## Make sure to enter the environment when starting
RUN echo "source activate $(head -1 /tmp/environment.yml | cut -d' ' -f2)" > ~/.bashrc
ENV PATH /opt/conda/envs/$(head -1 /tmp/environment.yml | cut -d' ' -f2)/bin:$PATH

WORKDIR /data
ENTRYPOINT [ "/bin/bash", "-l", "-c" ]

