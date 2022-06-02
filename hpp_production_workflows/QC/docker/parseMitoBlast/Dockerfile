FROM amancevice/pandas:latest
MAINTAINER Julian Lucas, juklucas@ucsc.edu

### PARSE MITO BLAST RESULTS
WORKDIR /opt
COPY scripts/parse_mito_blast_results/ /opt/parse_mito_blast_results/

RUN chmod +x /opt/parse_mito_blast_results/parse_mito_blast_results.py

ENV PATH="/opt/parse_mito_blast_results/:${PATH}"

WORKDIR /data
