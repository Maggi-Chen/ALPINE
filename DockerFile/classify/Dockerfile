FROM python:3.10-slim
LABEL description='Dockerfile for AAV-LR pipeline - Classify reads'

RUN apt-get update  && apt upgrade -y
RUN apt-get install -qy  build-essential 
RUN apt-get install -qy gcc-multilib apt-utils zlib1g-dev r-base wget git 

RUN apt-get install -qy python3 python3-pip zip unzip vim

WORKDIR /opt
RUN wget https://github.com/samtools/samtools/releases/download/1.19.2/samtools-1.19.2.tar.bz2

RUN tar xvjf samtools-1.19.2.tar.bz2
WORKDIR /opt/samtools-1.19.2/
RUN ./configure
RUN make
RUN make install

RUN mkdir /opt/blast && \
    cd /opt/blast && \
    wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz && \
    tar -xvzf ncbi-blast-2.16.0+-x64-linux.tar.gz 
WORKDIR /opt/blast
ENV PATH="/opt/blast/ncbi-blast-2.16.0+/bin:${PATH}"
RUN blastn -help

WORKDIR /opt
COPY . /opt/

RUN python3 -m venv env && \
    env/bin/python3 -m pip install --upgrade pip && \
    env/bin/pip install -r requirements.txt

ENV PATH /opt/env/bin:/opt:$PATH

WORKDIR /opt
RUN git clone https://github.com/lh3/minimap2

RUN apt-get update && \
    apt-get install -y \
    bwa \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/minimap2
RUN make
ENV PATH="$PATH:/opt/minimap2/"

WORKDIR /opt
RUN chmod 777 /opt/classify_read.py

RUN python classify_read.py --help


