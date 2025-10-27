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

LABEL maintainer="BC CfE in HIV/AIDS <https://github.com/cfe-lab/MiCall>"

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
    apt-get install -q -y ncbi-blast+ && \
    rm -rf /var/lib/apt/lists/*

## bowtie2
RUN wget -q -O bowtie2.zip https://github.com/BenLangmead/bowtie2/releases/download/v2.2.8/bowtie2-2.2.8-linux-x86_64.zip && \
  unzip bowtie2.zip -d /opt/ && \
  ln -s /opt/bowtie2-2.2.8/ /opt/bowtie2 && \
  rm bowtie2.zip

ENV PATH=$PATH:/opt/bowtie2

## Installing IVA dependencies
RUN apt-get install -q -y zlib1g-dev libncurses5-dev libncursesw5-dev && \
    cd /bin && \
    wget -q http://sun.aei.polsl.pl/kmc/download-2.1.1/linux/kmc && \
    wget -q http://sun.aei.polsl.pl/kmc/download-2.1.1/linux/kmc_dump && \
    chmod +x kmc kmc_dump && \
    cd /opt && \
    wget -q https://sourceforge.net/projects/mummer/files/mummer/3.23/MUMmer3.23.tar.gz && \
    tar -xzf MUMmer3.23.tar.gz --no-same-owner && \
    cd MUMmer3.23 && \
    make --quiet install && \
    rm -r docs src ../MUMmer3.23.tar.gz && \
    ln -s /opt/MUMmer3.23/nucmer \
        /opt/MUMmer3.23/delta-filter \
        /opt/MUMmer3.23/show-coords \
        /bin && \
    cd /opt && \
    wget -q https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 && \
    tar -xf samtools-1.3.1.tar.bz2 --no-same-owner --bzip2 && \
    cd samtools-1.3.1 && \
    ./configure --quiet --prefix=/ && \
    make --quiet && \
    make --quiet install && \
    cd /opt && \
    rm -rf samtools-1.3.1* && \
    wget -q http://downloads.sourceforge.net/project/smalt/smalt-0.7.6-bin.tar.gz && \
    tar -xzf smalt-0.7.6-bin.tar.gz --no-same-owner && \
    ln -s /opt/smalt-0.7.6-bin/smalt_x86_64 /bin/smalt
RUN pip install 'git+https://github.com/cfe-lab/iva.git@v1.1.1'

## Install dependencies for genetracks/drawsvg
RUN apt-get install -q -y libcairo2-dev
RUN pip install --upgrade pip setuptools

## Install just the dependencies of MiCall (for faster build times in development).
COPY pyproject.toml README.md /opt/micall/
RUN pip install /opt/micall[basespace]

## Trigger matplotlib to build its font cache
RUN python -c 'import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot'

COPY . /opt/micall/

RUN pip install /opt/micall[basespace]
RUN micall make_blast_db

WORKDIR /data
ENTRYPOINT ["micall", "analyze"]
