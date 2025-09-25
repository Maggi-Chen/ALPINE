FROM python:3.8-slim
LABEL description='Dockerfile for AAV-LR pipeline - Samtools'

RUN apt-get update
RUN apt-get install -qy  build-essential
RUN apt-get install -qy gcc-multilib apt-utils zlib1g-dev r-base wget 

RUN apt-get install -qy python3 python3-pip zip unzip

WORKDIR /opt
RUN wget https://github.com/samtools/samtools/releases/download/1.19.2/samtools-1.19.2.tar.bz2

RUN tar xvjf samtools-1.19.2.tar.bz2
WORKDIR /opt/samtools-1.19.2/
RUN ./configure
RUN make
RUN make install

