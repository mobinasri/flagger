FROM ubuntu:18.04
MAINTAINER Mobin Asri, masri@ucsc.edu

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -y build-essential git zlib1g-dev wget libbz2-dev pkg-config libncurses5-dev libncursesw5-dev liblzma-dev libcurl4-openssl-dev python3-dev pigz zip

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
	make
ENV PATH="/home/apps/yak:$PATH"

# Hifiasm-0.18.2
RUN cd /home/apps && \
	wget https://github.com/chhylp123/hifiasm/archive/refs/tags/0.18.2.zip && \
	unzip 0.18.2.zip  && \
        cd hifiasm-0.18.2 && \ 
	make
ENV PATH="/home/apps/hifiasm-0.18.2:$PATH"

RUN cd /home/apps && \
	git clone --branch v0.4 https://github.com/lh3/gfatools && \
	cd gfatools && \
	make
ENV PATH="/home/apps/gfatools:$PATH"

COPY ./scripts /home/scripts

