FROM mobinasri/bio_base:dev-v0.1
MAINTAINER Mobin Asri, masri@ucsc.edu

RUN apt-get update && apt-get install -y zip unzip curl
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
    wget https://github.com/lh3/minimap2/archive/refs/tags/v2.24.zip && \
    unzip v2.24.zip && \
    rm v2.24.zip && \
    cd minimap2-2.24 && \
    #cd minimap2 && \
    make
ENV PATH="/home/apps/minimap2-2.24:${PATH}"
#ENV PAFTOOLS_PATH="/home/apps/minimap2-2.21/misc/paftools.js"
#ENV PATH="/home/apps/minimap2:${PATH}"
ENV PAFTOOLS_PATH="/home/apps/minimap2-2.24/misc/paftools.js"

## k8
RUN cd /home/apps && \
    curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf -
RUN cd /home/apps/k8-0.2.4/ && cp k8-`uname -s` k8
ENV PATH="/home/apps/k8-0.2.4:$PATH"

## UCSC convention is to work in /data
RUN mkdir -p /data
WORKDIR /data

