# Dockerized version of MiCall.

# This Docker container will be used by BaseSpace, as well as for running on a
# Docker-enabled workstation.  Instructions may be found by running with the
# --help option; e.g. if you build this Docker image as micall:version, call
#
# docker run micall:version --help
#
# or for help on the specific subcommands,
#
# docker run micall:version {basespace,folder,sample,hcv_sample} --help
#
# This Dockerfile can be used to build two types of MiCall images:
# - a "production" image, which can be used to deploy and run Micall; and
# - a "dev" image, which contains packages needed for testing
#   and development of MiCall.
# The dev image is slower to build.
#
# To specify which image you want to build, use the `--target` tag to
# `docker build`, e.g.
#
# docker build --target production -t [image name]:[tag] [source directory]
#
# If you omit the `--target` tag altogether, `docker build` will build
# the development image.

FROM python:3.11

MAINTAINER BC CfE in HIV/AIDS https://github.com/cfe-lab/MiCall

## Download package sources
RUN apt-get update -qq -y

## Prerequisites
RUN apt-get update -qq --fix-missing && apt-get install -qq -y \
  unzip \
  wget \
  && rm -rf /var/lib/apt/lists/*

RUN wget -qO rustup.sh https://sh.rustup.rs && \
    chmod +x /rustup.sh && \
    /rustup.sh -y -q && \
    . /root/.cargo/env && \
    rm rustup.sh && \
    cargo install --root / --git https://github.com/jeff-k/merge-mates.git --rev 2fec61363f645e2008a4adff553d098beae21469

## Installing blast
RUN apt-get update -qq --fix-missing && \
    apt-get install -q -y ncbi-blast+

## bowtie2
RUN wget -q -O bowtie2.zip https://github.com/BenLangmead/bowtie2/releases/download/v2.2.8/bowtie2-2.2.8-linux-x86_64.zip && \
  unzip bowtie2.zip -d /opt/ && \
  ln -s /opt/bowtie2-2.2.8/ /opt/bowtie2 && \
  rm bowtie2.zip

ENV PATH $PATH:/opt/bowtie2

## Install Haploflow
RUN apt-get install -y build-essential sudo git ronn cmake && \
    cd /opt/ && \
    git clone https://github.com/hzi-bifo/Haploflow && \
    cd Haploflow && \
    git checkout 9a5a0ff6c3a0435e723e41f98fe82ec2ad19cf50 && \
    yes | sh build.sh && \
    ln -s /opt/Haploflow/build/haploflow /bin/haploflow

## Install dependencies for genetracks/drawsvg
RUN apt-get install -q -y libcairo2-dev
RUN pip install --upgrade pip

## Install just the dependencies of MiCall (for faster build times in development).
COPY pyproject.toml README.md /opt/micall/
RUN pip install /opt/micall[denovo,basespace]

## Trigger matplotlib to build its font cache
RUN python -c 'import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot'

COPY . /opt/micall/

RUN pip install /opt/micall[denovo,basespace]
RUN micall make_blast_db

WORKDIR /data
ENTRYPOINT ["micall", "micall_docker"]
