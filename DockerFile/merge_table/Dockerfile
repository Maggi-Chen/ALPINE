FROM python:3.8-slim
LABEL description='Dockerfile for AAV-LR pipeline - Merge count table'

RUN apt update -y && apt upgrade -y
RUN apt install -y curl unzip gnupg2 libldap2-dev libsasl2-dev gcc

RUN apt update -y && apt upgrade -y
RUN apt install -y \
  build-essential \
  git \
  libcurl4-openssl-dev \
  libssl-dev

WORKDIR /opt/

# Copy local files and shared module
# Note: Build with context as DockerFile/ directory:
#   docker build -f merge_table/Dockerfile -t aavlr_merge .
COPY merge_table/ /opt/
COPY shared/column_order.py /opt/column_order.py

RUN python3 -m venv env && \
    env/bin/python3 -m pip install --upgrade pip && \
    env/bin/pip install -r requirements.txt

ENV PATH /opt/env/bin:/opt:$PATH

CMD python3 /opt/merge_count_table.py
RUN chmod 777 /opt/merge_count_table.py
