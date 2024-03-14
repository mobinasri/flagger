FROM tpesout/hpp_base:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

### blast ###
WORKDIR /opt/blast
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/ncbi-blast-2.11.0+-src.tar.gz && \
    tar -xvf ncbi-blast-2.11.0+-src.tar.gz && \
    rm ncbi-blast-2.11.0+-src.tar.gz && \
    cd ncbi-blast-2.11.0+-src/c++/ && \
    ./configure && \
    make -j 8
ENV PATH="$PATH:/root/bin/samtools_1.9:/opt/blast/ncbi-blast-2.11.0+-src/c++/ReleaseMT/bin"

### vecscreen ###
RUN apt-get install -y ncbi-tools-bin gawk

### blast taxonomy db ###
WORKDIR /opt/blast/taxdb
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz && \
    tar xvf taxdb.tar.gz && \
    rm taxdb.tar.gz
ENV BLASTDB=/opt/blast/taxdb

### vecscreen tool from ncbi
WORKDIR /usr/bin
RUN wget ftp://ftp.ncbi.nlm.nih.gov/pub/kitts/VSlistTo1HitPerLine.awk && \
    chmod +x VSlistTo1HitPerLine.awk

### special blast version needed for makeblastdb and vecscreen
WORKDIR /opt/blast
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/ncbi-blast-2.7.1+-x64-linux.tar.gz && \
    tar xvf ncbi-blast-2.7.1+-x64-linux.tar.gz && \
    rm ncbi-blast-2.7.1+-x64-linux.tar.gz

### use python3, we have some scripts that need this
ENV PATH="/opt/python/Python-3.6.0/bin:${PATH}"
RUN ln -s /opt/python/Python-3.6.0/bin/python3 /usr/bin/python && \
    /opt/python/Python-3.6.0/bin/pip3 install biopython==1.78

### bedtools
WORKDIR /opt/bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz && \
    tar -zxvf bedtools-2.28.0.tar.gz && \
    rm -r /opt/bedtools/bedtools-2.28.0.tar.gz && \
    cd bedtools2 && \
    make && \
    ./bin/bedtools
ENV PATH="$PATH:/opt/bedtools/bedtools2/bin/"

### scripts to chunk/unchunk assembly, vecscreen and rrna require smaller query
COPY chunk_assembly.py /usr/bin
COPY unchunk_vecscreen.py /usr/bin
COPY unchunk_blast.py /usr/bin
### a script which saves species' names to output
COPY update_megablast_output_with_subject_description.py /usr/bin


WORKDIR /data
