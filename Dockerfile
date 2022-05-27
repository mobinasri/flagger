FROM mobinasri/bio_base:latest
MAINTAINER Mobin Asri, masri@ucsc.edu

RUN mkdir -p /home/apps
RUN pip3 install scipy pandas matplotlib
RUN apt-get update
RUN apt-get install -y build-essential python3-dev autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev wget unzip


#install Java
RUN cd /home/apps && \
    wget https://download.oracle.com/java/17/latest/jdk-17_linux-x64_bin.tar.gz && \
    tar -xvzf jdk-17_linux-x64_bin.tar.gz
ENV PATH=$PATH:/home/apps/jdk-17.0.1/bin

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


#intstall IGV
RUN cd /home/apps && \
    wget https://data.broadinstitute.org/igv/projects/downloads/2.11/IGV_2.11.3.zip &&\
    unzip IGV_2.11.3.zip
ENV PATH=$PATH:/home/apps/IGV_2.11.3

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

#install tabix
RUN apt-get install -y tabix 

#install whatshap
ENV LANG="C.UTF-8" 
RUN pip3 install cython && pip3 install whatshap

COPY ./programs /home/programs
# Add cigar_it to the submodules dir
COPY ./ext/secphase/programs/submodules/cigar_it /home/programs/submodules
COPY ./scripts  /home/scripts
RUN cd /home/programs && make
ENV PATH="$PATH:/home/programs/bin"

ENV CALC_MEAN_SD_PY=/home/programs/src/calc_mean_sd.py
ENV FIT_MODEL_PY=/home/programs/src/fit_model.py
ENV PROJECT_BLOCKS_PY=/home/programs/src/project_blocks.py
ENV ASSIGN_GAPS_PY=/home/programs/src/assign_gaps.py
ENV FIT_MODEL_EXTRA_PY=/home/programs/src/fit_model_extra.py
ENV PDF_GENERATOR_PY=/home/programs/src/pdf_generator.py
ENV PARTITION_SECPHASE_READS_PY=/home/programs/src/partition_secphase_reads.py

## UCSC convention is to work in /data
RUN mkdir -p /data
WORKDIR /data
