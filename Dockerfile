# Generate the Docker container to run MiCall on BaseSpace.
FROM python:3.7

MAINTAINER BC CfE in HIV/AIDS https://github.com/cfe-lab/MiCall

## Prerequisites
RUN apt-get update -qq --fix-missing && apt-get install -qq -y \
  unzip \
  wget \
  libpython3.7-dev \
  && rm -rf /var/lib/apt/lists/*

RUN wget -qO rustup.sh https://sh.rustup.rs && \
    chmod +x /rustup.sh && \
    /rustup.sh -y -q && \
    . /root/.cargo/env && \
    rm rustup.sh && \
    cargo install --root / --git https://github.com/jeff-k/merge-mates.git

## Installing blast
RUN apt-get update -qq --fix-missing && \
    apt-get install -q -y ncbi-blast+ && \
    rm -rf /var/lib/apt/lists/*

## bowtie2
RUN wget -q -O bowtie2.zip https://github.com/BenLangmead/bowtie2/releases/download/v2.2.8/bowtie2-2.2.8-linux-x86_64.zip && \
  unzip bowtie2.zip -d /opt/ && \
  ln -s /opt/bowtie2-2.2.8/ /opt/bowtie2 && \
  rm bowtie2.zip

ENV PATH $PATH:/opt/bowtie2

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

## Install dependencies for genetracks/drawSvg
RUN apt-get install -q -y libcairo2-dev

## Gotoh
COPY micall/alignment /opt/micall/alignment
COPY requirements.txt requirements-basespace.txt /opt/micall/

## Python packages, plus trigger matplotlib to build its font cache
WORKDIR /opt
RUN pip install --upgrade pip && \
  pip install -r /opt/micall/requirements-basespace.txt && \
  python -c 'import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot'

## MiCall
COPY micall_basespace.py micall_kive.py micall_kive_resistance.py version.tx[t] /opt/micall/
COPY micall/__init__.py micall/project* /opt/micall/micall/

COPY micall/blast_db/make_blast_db.py    /opt/micall/micall/blast_db/make_blast_db.py
COPY micall/core    /opt/micall/micall/core/
COPY micall/data    /opt/micall/micall/data/
COPY micall/drivers    /opt/micall/micall/drivers/
COPY micall/g2p     /opt/micall/micall/g2p/
COPY micall/resistance   /opt/micall/micall/resistance/
COPY micall/monitor /opt/micall/micall/monitor/
COPY micall/utils   /opt/micall/micall/utils/

RUN python /opt/micall/micall/blast_db/make_blast_db.py

CMD ["python", "/opt/micall/micall_basespace.py", "--link_run", "/input"]
