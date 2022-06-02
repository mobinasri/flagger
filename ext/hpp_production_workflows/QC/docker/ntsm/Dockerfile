FROM ubuntu:18.04

MAINTAINER Julian Lucas, juklucas@ucsc.edu


RUN apt-get update && apt-get install -y \
	python3.8 \
	python3.8-dev \
	python3.8-distutils \
	python3-pip \
    autotools-dev \
    automake \
    git \
    pigz \
    samtools

RUN ln -s /usr/bin/python3 /usr/bin/python & \
    ln -s /usr/bin/pip3 /usr/bin/pip

RUN pip install pyfaidx	


## Install ntsm
WORKDIR /opt/

RUN git clone https://github.com/JustinChu/ntsm.git \
	&& cd ntsm \
	&& git reset --hard b95e4c6f8d73a95a1dff99c6f01bd7d26821981e \
	&& git submodule update --init \
	&& ./autogen.sh \
	&& ./configure \
	&& make

ENV PATH=/opt/ntsm:$PATH

WORKDIR /data
