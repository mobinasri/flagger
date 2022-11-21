FROM mobinasri/bio_base:dev-v0.1
MAINTAINER Mobin Asri, masri@ucsc.edu

RUN mkdir -p /home/apps
RUN pip3 install scipy pandas matplotlib
RUN apt-get update
RUN apt-get install -y build-essential python3-dev autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev wget unzip


#install sonLib
RUN cd /home/apps && \
    git clone https://github.com/benedictpaten/sonLib && \
    cd sonLib && \
    make

#install hstlib  C API
RUN cd /home/apps && \
    wget https://github.com/samtools/htslib/releases/download/1.13/htslib-1.13.tar.bz2 && \
    tar -xvjf htslib-1.13.tar.bz2 && \
    cd htslib-1.13 && \
    autoreconf -i && \
    ./configure && \
    make && \
    make install
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib


COPY ./programs /home/programs
#COPY ./scripts  /home/scripts
RUN cd /home/programs && make
ENV PATH="$PATH:/home/programs/bin"

ENV PARTITION_SECPHASE_READS_PY=/home/programs/src/partition_secphase_reads.py

## UCSC convention is to work in /data
RUN mkdir -p /data
WORKDIR /data
