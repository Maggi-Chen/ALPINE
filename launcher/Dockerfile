FROM python:3.11-slim
LABEL description='Dockerfile for Transgene Launcher'

RUN apt update -y && apt upgrade -y
RUN apt install -y \
  build-essential \
  git \
  libcurl4-openssl-dev \
  libssl-dev

WORKDIR /opt/launcher/

COPY Dockerfile /opt/
COPY requirements.txt ./
RUN python3 -m venv env && \
    env/bin/python3 -m pip install --upgrade pip && \
    env/bin/pip install -r requirements.txt

COPY src/*.py ./

ENV PATH /opt/launcher/env/bin:/opt/launcher:$PATH

RUN chmod +x /opt/launcher/launcher.py
CMD python3 /opt/launcher/launcher.py

