FROM python:3.8

MAINTAINER Julian Lucas, juklucas@ucsc.edu

LABEL    software="ont_summary_stats" \ 
    container="ont_summary_stats" \ 
    about.summary="Container with script to calculated coverage metrics from Oxford Nanopore's sequencing summary stats files."

RUN apt-get update

### Script courtesy of Miten Jain
WORKDIR /opt
COPY scripts/calculate_summary_stats.py /opt/calculate_summary_stats.py

RUN chmod +x /opt/calculate_summary_stats.py && \
    pip3 install numpy


WORKDIR /data