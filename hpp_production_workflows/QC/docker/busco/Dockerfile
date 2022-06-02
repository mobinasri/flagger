FROM tpesout/hpp_base:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

ENV PATH="/root/bin/python_3.6.0/bin:${PATH}"

RUN apt-get update && \
    apt-get install -y fort77 gfortran build-essential xorg-dev  liblzma-dev libblas-dev gcc-multilib gobjc++ aptitude \
    libreadline-dev bamtools libbamtools-dev libboost-dev

### R
# 3.2.1
WORKDIR /opt/R
RUN wget https://cran.r-project.org/src/base/R-3/R-3.2.1.tar.gz && \
    tar xvf R-3.2.1.tar.gz && \
    rm -r /opt/R/R-3.2.1.tar.gz && \
    cd R-3.2.1 && \
    ./configure --enable-utf8  --prefix /root/bin/R_3.2.1 && \
    make && \
    make install

### blast
# 2.2.30+
WORKDIR /opt/blast
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.30/ncbi-blast-2.2.30+-x64-linux.tar.gz && \
    tar -xvf ncbi-blast-2.2.30+-x64-linux.tar.gz && \
    rm -r /opt/blast/ncbi-blast-2.2.30+-x64-linux.tar.gz && \
    mkdir /root/bin/blast_2.2.30+ && \
    mv /opt/blast/ncbi-blast-2.2.30+/bin/* /root/bin/blast_2.2.30+/

### hmmer
# 3.2.1
WORKDIR /opt/hmmer
RUN wget https://github.com/EddyRivasLab/hmmer/archive/hmmer-3.2.1.tar.gz && \
    tar -xvf hmmer-3.2.1.tar.gz && \
    rm -r /opt/hmmer/hmmer-3.2.1.tar.gz && \
    cd hmmer-hmmer-3.2.1 && \
    git clone https://github.com/EddyRivasLab/easel.git && \
    cd easel && \
    git checkout 924d7efff9a765d1d5807904f9686e2513a1e4f2 && \
    cd .. && \
    autoconf && \
    ./configure --prefix /root/bin/hmmer_3.2.1/ && \
    make && \
    make install

### BUSCO
# 3.0.2
WORKDIR /root/tools/BUSCO
RUN wget https://gitlab.com/ezlab/busco/-/archive/3.0.2/busco-3.0.2.tar.gz && \
    tar xvf busco-3.0.2.tar.gz && \
    mv busco-3.0.2 busco/ && \
    rm busco-3.0.2.tar.gz && \
    cd busco && \
    python3 setup.py install
COPY config.ini /root/tools/BUSCO/busco/config/config.ini
# vertebrata dataset
WORKDIR /root/tools/BUSCO/dataset/
RUN wget http://busco.ezlab.org/v2/datasets/vertebrata_odb9.tar.gz && \
    tar xvf vertebrata_odb9.tar.gz && \
    rm vertebrata_odb9.tar.gz

### Augustus
# 3.3.1
WORKDIR /root/tools/Augustus/
RUN wget https://github.com/Gaius-Augustus/Augustus/archive/v3.3.1-tag1.tar.gz && \
    tar xvf v3.3.1-tag1.tar.gz && \
    mv Augustus-3.3.1-tag1/ Augustus-3.3.1/ && \
    rm v3.3.1-tag1.tar.gz && \
    cd Augustus-3.3.1/ && \
    mkdir tooldir && \
    ln -s /opt/samtools/samtools-1.9/ tooldir/samtools  && \
    ln -s /opt/bcftools/bcftools-1.9/ tooldir/bcftools && \
    ln -s /opt/htslib/htslib-1.9/ tooldir/htslib && \
    TOOLDIR=/root/tools/Augustus/Augustus-3.3.1/tooldir make


ENV PATH="/root/tools/Augustus/Augustus-3.3.1/bin:/root/tools/Augustus/Augustus-3.3.1/scripts:$PATH"

WORKDIR /data
