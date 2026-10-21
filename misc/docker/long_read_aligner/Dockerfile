ARG BIO_BASE=mobinasri/bio_base:v1.2.1
FROM ${BIO_BASE}
MAINTAINER Mobin Asri, masri@ucsc.edu

# TARGETARCH (amd64/arm64) is set by BuildKit; fall back to dpkg for the legacy builder
ARG TARGETARCH

RUN pip3 install scipy pandas matplotlib
RUN apt-get update
RUN apt-get install -y build-essential curl python3-dev autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev wget unzip zip
RUN mkdir -p /home/apps

## Winnowmap
RUN ARCH=${TARGETARCH:-$(dpkg --print-architecture)} && \
    cd /home/apps && \
    wget https://github.com/marbl/Winnowmap/archive/refs/tags/v2.03.zip && \
    unzip v2.03.zip && \
    rm v2.03.zip && \
    cd Winnowmap-2.03 && \
    if [ "${ARCH}" = "arm64" ]; then make -j8 arm_neon=1 aarch64=1; else make -j8; fi
ENV PATH="/home/apps/Winnowmap-2.03/bin:${PATH}"

## minimap2
RUN ARCH=${TARGETARCH:-$(dpkg --print-architecture)} && \
    cd /home/apps && \
    wget https://github.com/lh3/minimap2/archive/refs/tags/v2.31.zip && \
    unzip v2.31.zip && \
    rm v2.31.zip && \
    cd minimap2-2.31 && \
    if [ "${ARCH}" = "arm64" ]; then make arm_neon=1 aarch64=1; else make; fi
ENV PATH="/home/apps/minimap2-2.31:${PATH}"
ENV PAFTOOLS_PATH="/home/apps/minimap2-2.31/misc/paftools.js"

## k8
# k8 has no prebuilt Linux binary for arm64 so on arm64 it is installed from bioconda with micromamba
# (paftools.js supports both k8 0.x and 1.x)
RUN ARCH=${TARGETARCH:-$(dpkg --print-architecture)} && \
    cd /home/apps && \
    mkdir -p k8 && \
    if [ "${ARCH}" = "arm64" ]; then \
        micromamba create -y -p /home/apps/conda_env -c conda-forge -c bioconda k8=1.2 && \
        micromamba clean -a -y && \
        ln -s /home/apps/conda_env/bin/k8 k8/k8; \
    else \
        curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf - && \
        cp k8-0.2.4/k8-`uname -s` k8/k8; \
    fi
ENV PATH="/home/apps/k8:$PATH"


#install veritymap
RUN apt install -y software-properties-common
RUN add-apt-repository -y ppa:ubuntu-toolchain-r/test && \
    apt install -y gcc-9 g++-9
ENV CC=/usr/bin/gcc-9
ENV CXX=/usr/bin/g++-9

# VerityMap uses x86 SSE intrinsics and does not build on arm64,
# so on arm64 VERITY_MAP_PY points to a script that exits with an error
RUN ARCH=${TARGETARCH:-$(dpkg --print-architecture)} && \
    cd /home/apps && \
    if [ "${ARCH}" = "arm64" ]; then \
        mkdir -p VerityMap-2.1.2-alpha/veritymap && \
        printf 'import sys\nsys.exit("ERROR: VerityMap is not available in the linux/arm64 image; use the linux/amd64 image or another aligner (minimap2/winnowmap).")\n' > VerityMap-2.1.2-alpha/veritymap/main.py; \
    else \
        wget https://github.com/ablab/VerityMap/archive/refs/tags/v2.1.2-alpha.zip && \
        unzip v2.1.2-alpha.zip && \
        cd VerityMap-2.1.2-alpha/ && \
        pip install -r requirements.txt && \
        cd veritymap && \
        make; \
    fi

ENV VERITY_MAP_PY="/home/apps/VerityMap-2.1.2-alpha/veritymap/main.py"

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

## UCSC convention is to work in /data
RUN mkdir -p /data
WORKDIR /data

