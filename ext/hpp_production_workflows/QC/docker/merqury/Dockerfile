FROM ubuntu:20.04
MAINTAINER Julian Lucas, juklucas@ucsc.edu


# 1. Set base ubuntu env
RUN apt update --fix-missing && \
    DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata && \
    apt upgrade -y && \
    apt install -y build-essential=12.8ubuntu1 default-jre=2:1.11-72 bedtools=2.27.1+dfsg-4ubuntu1 samtools=1.10-3 r-base=3.6.3-2 && \
    apt install -y vim htop git wget pigz

# 2. Install R packages
RUN Rscript -e 'install.packages("argparse", version="2.0.1")' && \
    Rscript -e 'install.packages("ggplot2", version="3.3.2")' && \
    Rscript -e 'install.packages("scales", version="1.1.1")'

# 3. Install IGV
WORKDIR /opt
RUN wget https://data.broadinstitute.org/igv/projects/downloads/2.8/IGV_2.8.2.zip && \
    unzip IGV_2.8.2.zip && \
    rm IGV_2.8.2.zip
ENV PATH=/opt/IGV_2.8.2/:$PATH

# 4. Get Meryl
WORKDIR /opt
RUN git clone https://github.com/marbl/meryl.git && \
    cd meryl && \
    git reset --hard 1926400fcb098ab3e4745ea0872d23f91d730340 && \
    cd src && \
    make -j 12
ENV PATH=/opt/meryl/build/bin:$PATH

# 5. Get Merqury
WORKDIR /opt
RUN git clone https://github.com/marbl/merqury.git && \
    cd merqury && \
    git reset --hard 136a9a349af4000831dd15c192e82ad3e0c14c83
ENV MERQURY=/opt/merqury
ENV PATH=/opt/merqury:$PATH

WORKDIR /data

