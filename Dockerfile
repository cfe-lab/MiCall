# Generate the Docker container to run MiCall on BaseSpace.
FROM python:2.7

MAINTAINER BC CfE in HIV/AIDS https://github.com/cfe-lab/MiCall

## Prerequisites
RUN apt-get update -qq --fix-missing && apt-get install -qq -y \
  unzip \
  wget \
  && rm -rf /var/lib/apt/lists/*

## bowtie2
RUN wget -q -O bowtie2.zip https://github.com/BenLangmead/bowtie2/releases/download/v2.2.8/bowtie2-2.2.8-linux-x86_64.zip; \
  unzip bowtie2.zip -d /opt/; \
  ln -s /opt/bowtie2-2.2.8/ /opt/bowtie2; \
  rm bowtie2.zip

ENV PATH $PATH:/opt/bowtie2

## Gotoh
COPY micall/alignment /opt/alignment
WORKDIR /opt/alignment
RUN python setup.py install
WORKDIR /
RUN rm -r /opt/alignment

## MiCall
COPY micall_basespace.py /opt/micall/
COPY micall/__init__.py micall/project* /opt/micall/micall/
COPY micall/core /opt/micall/micall/core/
COPY micall/utils /opt/micall/micall/utils/
