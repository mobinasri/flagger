FROM tpesout/hpp_base:latest
MAINTAINER Julian Lucas, juklucas@ucsc.edu

WORKDIR /opt

RUN wget https://github.com/lh3/minigraph/archive/refs/tags/v0.15.tar.gz && \
	tar -xvf v0.15.tar.gz && \
	rm v0.15.tar.gz && \
	cd minigraph-0.15 && \
	make

ENV PATH="/opt/minigraph-0.15:$PATH"

WORKDIR /data
