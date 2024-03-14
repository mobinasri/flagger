FROM tpesout/hpp_base:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

RUN apt-get update && \
    apt-get -y install software-properties-common && \
    add-apt-repository ppa:ubuntugis/ubuntugis-unstable && \
    apt-get update && \
    apt-get -y install libc-dev bzip2 libcurl4-openssl-dev openjdk-8-jdk ant fort77 gfortran build-essential xorg-dev \
        liblzma-dev libblas-dev gcc-multilib gobjc++ aptitude libreadline-dev r-base xvfb gsfonts-x11 xfonts-base \
        xfonts-scalable xfonts-100dpi xfonts-75dpi libssl-dev libxml2-dev libudunits2-dev libcairo2-dev \
        libudunits2-dev libgdal-dev libgeos-dev libproj-dev


### R
# 3.2.1
WORKDIR /opt/R
RUN wget https://cran.r-project.org/src/base/R-3/R-3.6.3.tar.gz && \
    tar xvf R-3.6.3.tar.gz && \
    rm -r /opt/R/R-3.6.3.tar.gz && \
    cd R-3.6.3 && \
    ./configure --enable-utf8  --prefix /root/bin/R_3.6.3 && \
    make && \
    make install && \
  	/root/bin/R_3.6.3/bin/R -e "install.packages('argparse', repos ='http://cran.us.r-project.org', dependencies=TRUE)" && \
  	/root/bin/R_3.6.3/bin/R -e "install.packages('ggplot2', repos ='http://cran.us.r-project.org', dependencies=TRUE)" && \
  	/root/bin/R_3.6.3/bin/R -e "install.packages('scales', repos ='http://cran.us.r-project.org', dependencies=TRUE)" && \
  	/root/bin/R_3.6.3/bin/R -e "install.packages('gplots', repos ='http://cran.us.r-project.org', dependencies=TRUE)"
ENV PATH="/root/bin/R_3.6.3/bin:${PATH}"

### mash
WORKDIR /opt/mash
RUN wget https://github.com/marbl/Mash/releases/download/v2.2/mash-Linux64-v2.2.tar && \
    tar xvf mash-Linux64-v2.2.tar && \
    rm mash-Linux64-v2.2.tar && \
    ./mash-Linux64-v2.2/mash && \
    wget https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh && \
    wget https://raw.githubusercontent.com/VGP/vgp-assembly/c1214afa52c9588e8fa0d4ec212986d5151fef26/pipeline/mash/plot.R
ENV PATH="/opt/mash/mash-Linux64-v2.2:$PATH"
ENV MASH_REFSEQ="/opt/mash/refseq.genomes.k21s1000.msh"
ENV PLOT_R="/opt/mash/plot.R"

WORKDIR /data