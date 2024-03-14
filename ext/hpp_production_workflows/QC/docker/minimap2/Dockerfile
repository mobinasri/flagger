FROM tpesout/hpp_base:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

WORKDIR /opt/minimap2
RUN wget https://github.com/lh3/minimap2/releases/download/v2.23/minimap2-2.23.tar.bz2 && \
	tar xvf minimap2-2.23.tar.bz2 && \
	rm minimap2-2.23.tar.bz2 && \
	cd minimap2-2.23 && \
	make
ENV PATH="/opt/minimap2/minimap2-2.23:/opt/minimap2/minimap2-2.23/misc:$PATH"

WORKDIR /opt/k8
RUN wget https://github.com/attractivechaos/k8/releases/download/0.2.5/k8-0.2.5.tar.bz2 && \
    tar xvf k8-0.2.5.tar.bz2 && \
    rm k8-0.2.5.tar.bz2 && \
    cp k8-0.2.5/k8-`uname -s` k8-0.2.5/k8
ENV PATH="/opt/k8/k8-0.2.5:$PATH"

WORKDIR /data
