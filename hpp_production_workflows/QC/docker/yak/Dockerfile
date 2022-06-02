FROM tpesout/hpp_base:latest
MAINTAINER Julian Lucas, juklucas@ucsc.edu

WORKDIR /opt

RUN git clone https://github.com/lh3/yak && \
	cd yak && \
	git reset --hard b17412df1bcd0121414b6d62e332a00b1a5e11cc && \
	make
ENV PATH="/opt/yak:$PATH"

WORKDIR /data
