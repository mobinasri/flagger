FROM python:3.8

MAINTAINER Julian Lucas, juklucas@ucsc.edu

LABEL base_image="python:3.8"
LABEL software="s5cmd"
LABEL software.version="2.0"
LABEL about.home="https://github.com/peak/s5cmd"
LABEL about.license_file="https://github.com/peak/s5cmd/blob/master/LICENSE"
LABEL about.summary="s5cmd is a very fast S3 and local filesystem execution tool. This container is built to download from open buckets that do not charge egress as part of the HPRC's data management."

## Install s5cmd
WORKDIR /opt/

RUN curl -LO https://github.com/peak/s5cmd/releases/download/v2.0.0/s5cmd_2.0.0_Linux-64bit.tar.gz \
	&& tar -xvf s5cmd_2.0.0_Linux-64bit.tar.gz \
	&& rm s5cmd_2.0.0_Linux-64bit.tar.gz

ENV PATH=/opt:$PATH

WORKDIR /data
