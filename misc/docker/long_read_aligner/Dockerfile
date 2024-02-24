FROM mobinasri/bio_base:v0.4.0
MAINTAINER Mobin Asri, masri@ucsc.edu

RUN pip3 install scipy pandas matplotlib
RUN apt-get update
RUN apt-get install -y build-essential curl python3-dev autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev wget unzip zip
RUN mkdir -p /home/apps

## Winnowmap
RUN cd /home/apps && \
    wget https://github.com/marbl/Winnowmap/archive/refs/tags/v2.03.zip && \
    unzip v2.03.zip && \
    rm v2.03.zip && \
    cd Winnowmap-2.03 && \
    make -j8
ENV PATH="/home/apps/Winnowmap-2.03/bin:${PATH}"

## minimap2
RUN cd /home/apps && \
    #git clone -b 2.23 https://github.com/lh3/minimap2 && \
    wget https://github.com/lh3/minimap2/archive/refs/tags/v2.26.zip && \
    unzip v2.26.zip && \
    rm v2.26.zip && \
    cd minimap2-2.26 && \
    #cd minimap2 && \
    make
ENV PATH="/home/apps/minimap2-2.26:${PATH}"
ENV PAFTOOLS_PATH="/home/apps/minimap2-2.26/misc/paftools.js"

## k8
RUN cd /home/apps && \
    curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf -
RUN cd /home/apps/k8-0.2.4/ && cp k8-`uname -s` k8
ENV PATH="/home/apps/k8-0.2.4:$PATH"


#install veritymap
RUN apt install -y software-properties-common
RUN add-apt-repository -y ppa:ubuntu-toolchain-r/test && \
    apt install -y gcc-9 g++-9
ENV CC=/usr/bin/gcc-9
ENV CXX=/usr/bin/g++-9

RUN cd /home/apps && \
    wget https://github.com/ablab/VerityMap/archive/refs/tags/v2.1.2-alpha.zip && \
    unzip v2.1.2-alpha.zip && \
    cd VerityMap-2.1.2-alpha/ && \
    pip install -r requirements.txt && \
    cd veritymap && \
    make

ENV VERITY_MAP_PY="/home/apps/VerityMap-2.1.2-alpha/veritymap/main.py"

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

## UCSC convention is to work in /data
RUN mkdir -p /data
WORKDIR /data

