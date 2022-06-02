FROM ubuntu:16.04
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

# update and install dependencies
RUN apt-get update && \
    apt-get -y install time git make wget autoconf gcc g++ vim sudo build-essential bzip2 zlib1g-dev libbz2-dev \
        libcurl4-gnutls-dev liblzma-dev libncurses5-dev libncursesw5-dev libssl-dev

#    autopoint tcl8.6-dev python-dev \
#                       perl \
#                       gfortranfort77 xorg-dev libblas-dev gcc-multilib gobjc++ aptitude \
#                       libreadline-dev openjdk-8-jdk libpcre3 libpcre3-dev libbamtools-dev libsuitesparse-dev \
#                       liblpsolve55-dev libboost-iostreams-dev libhdf5-dev libboost-all-dev libpng-dev \
#                       graphviz sudo
#&& \
#    apt-get clean && \
#    apt-get purge && \
#    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

### cmake build tool
WORKDIR /opt/cmake_install
RUN mkdir /opt/cmake && \
    wget https://cmake.org/files/v3.13/cmake-3.13.5-Linux-x86_64.sh && \
    sh cmake-3.13.5-Linux-x86_64.sh --prefix=/opt/cmake --skip-license && \
    ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake && \
    rm -r /opt/cmake_install

### htslib
# 1.9
WORKDIR /opt/htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
    tar -xvf htslib-1.9.tar.bz2 && \
    rm -r /opt/htslib/htslib-1.9.tar.bz2 && \
    cd htslib-1.9 && \
    ./configure && \
    make

### samtools
# 1.9
WORKDIR /opt/samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar xvf samtools-1.9.tar.bz2 && \
    rm -r /opt/samtools/samtools-1.9.tar.bz2 && \
    cd samtools-1.9/ && \
    autoheader && \
    autoconf -Wno-header && \
    ./configure && \
    make && \
    mkdir -p /root/bin/samtools_1.9 && \
    cp samtools /root/bin/samtools_1.9

### bcftools
## 1.9
WORKDIR /opt/bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && \
    tar xvf bcftools-1.9.tar.bz2 && \
    rm -r /opt/bcftools/bcftools-1.9.tar.bz2 && \
    cd bcftools-1.9/ && \
    autoheader && \
    autoconf -Wno-header && \
    ./configure --without-curses && \
    make && \
    mkdir /root/bin/bcftools_1.9 && \
    cp bcftools /root/bin/bcftools_1.9

### python
# 2.7
WORKDIR /opt/python
RUN wget https://www.python.org/ftp/python/2.7.13/Python-2.7.13.tgz && \
    tar xvf Python-2.7.13.tgz && \
    rm /opt/python/Python-2.7.13.tgz && \
    cd Python-2.7.13 && \
    ./configure --with-ensurepip=install --prefix=/opt/python/Python-2.7.13 && \
    make install
# 3.6.0
WORKDIR /opt/python
RUN wget https://www.python.org/ftp/python/3.6.0/Python-3.6.0.tgz && \
    tar xvf Python-3.6.0.tgz && \
    rm /opt/python/Python-3.6.0.tgz && \
    cd Python-3.6.0 && \
    ./configure --with-ensurepip=install --prefix=/opt/python/Python-3.6.0  && \
    make install && \
    ln -s /opt/python/Python-3.6.0/bin/python3 /opt/python/Python-3.6.0/bin/python
RUN chmod -R +x /opt/python/Python-3.6.0/


WORKDIR /data