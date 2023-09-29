FROM ubuntu:20.04

MAINTAINER Jean Monlong, jean.monlong@gmail.com

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get -y install \
    cmake git make gcc g++ autoconf bzip2 wget tabix libz-dev libncurses5-dev libbz2-dev liblzma-dev && \
    apt-get clean

WORKDIR /opt/winnowmap

RUN wget https://github.com/marbl/Winnowmap/archive/refs/tags/v2.03.tar.gz && \
    tar -xzvf v2.03.tar.gz && \
    cd Winnowmap-2.03 && \
    make && \
    cp bin/winnowmap /usr/bin/winnowmap && \
    cp bin/meryl /usr/bin/meryl && \
    rm ../v2.03.tar.gz

### samtools
# 1.9
WORKDIR /opt/samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar xvf samtools-1.9.tar.bz2 && \
	rm -r /opt/samtools/samtools-1.9.tar.bz2 && \
	cd samtools-1.9/ && \
	autoheader && \
	autoconf -Wno-header && \
	./configure && \
	make && \
	cp samtools /usr/bin/samtools

WORKDIR /home
