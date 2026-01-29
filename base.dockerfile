#
# Docker image for MiCall's base system.
# Includes all required system dependencies.
#

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

WORKDIR /data
ENTRYPOINT ["micall", "analyze"]
