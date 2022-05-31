FROM mobinasri/bio_base:v0.1
MAINTAINER Mobin Asri, masri@ucsc.edu

RUN mkdir -p /home/apps

RUN cd /home/apps && \
	git clone https://github.com/lh3/seqtk.git && \
        cd seqtk && \
        make
ENV PATH="/home/apps/seqtk:$PATH"

## UCSC convention is to work in /data
RUN mkdir -p /data
WORKDIR /data
