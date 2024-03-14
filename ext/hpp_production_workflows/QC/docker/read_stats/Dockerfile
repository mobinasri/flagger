FROM mobinasri/bio_base:v0.3.0
MAINTAINER Julian Lucas, juklucas@ucsc.edu

RUN apt-get update

### FAI READ STATS
WORKDIR /opt
COPY scripts/fai_read_stats/ /opt/fai_read_stats/

RUN chmod +x /opt/fai_read_stats/fai_read_stats.py && \
    pip3 install matplotlib numpy

ENV PATH="/opt/fai_read_stats/:${PATH}"

WORKDIR /data
