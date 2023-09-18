FROM ubuntu:20.04

MAINTAINER Jean Monlong, jean.monlong@gmail.com

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get -y install \
    cmake git make gcc g++ autoconf bzip2 wget tabix libz-dev libncurses5-dev libbz2-dev liblzma-dev && \
    apt-get clean

WORKDIR /opt/ngmlr

RUN wget https://github.com/philres/ngmlr/releases/download/v0.2.7/ngmlr-0.2.7-linux-x86_64.tar.gz && \
    tar -xzvf ngmlr-0.2.7-linux-x86_64.tar.gz && \
    mv ngmlr-0.2.7/ngmlr /bin/ && \
    chmod +x /bin/ngmlr && \
    rm ngmlr-0.2.7-linux-x86_64.tar.gz

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
