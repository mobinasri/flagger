FROM tpesout/hpp_base:latest
MAINTAINER Julian Lucas, juklucas@ucsc.edu

ENV PATH="/root/bin/python_3.6.0/bin/:${PATH}"
ENV PATH="/root/bin/samtools_1.9/:${PATH}"

RUN apt-get update && \
    apt-get install -y python-matplotlib

### FAI READ STATS
WORKDIR /opt
COPY scripts/fai_read_stats/ /opt/fai_read_stats/

RUN chmod +x /opt/fai_read_stats/fai_read_stats.py && \
    pip3 install matplotlib numpy

ENV PATH="/opt/fai_read_stats/:${PATH}"

WORKDIR /data
