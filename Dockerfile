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

FROM debian:bookworm-slim

LABEL maintainer="BC CfE in HIV/AIDS <https://github.com/cfe-lab/MiCall>"

## Prerequisites
RUN apt-get update -qq --fix-missing && apt-get install -qq -y \
  unzip \
  wget build-essential

RUN wget -qO rustup.sh https://sh.rustup.rs && \
    chmod +x /rustup.sh && \
    /rustup.sh -y -q && \
    . /root/.cargo/env && \
    rm rustup.sh && \
    cargo install --root / --git https://github.com/jeff-k/merge-mates.git --rev 2fec61363f645e2008a4adff553d098beae21469

## Installing blast
RUN apt-get install -q -y ncbi-blast+

## bowtie2
RUN wget -q -O bowtie2.zip https://github.com/BenLangmead/bowtie2/releases/download/v2.2.8/bowtie2-2.2.8-linux-x86_64.zip && \
  unzip bowtie2.zip -d /opt/ && \
  ln -s /opt/bowtie2-2.2.8/ /opt/bowtie2 && \
  rm bowtie2.zip

ENV PATH=$PATH:/opt/bowtie2

## Install Haploflow
RUN apt-get install -y build-essential sudo git ronn cmake && \
    cd /opt/ && \
    git clone https://github.com/hzi-bifo/Haploflow && \
    cd Haploflow && \
    git checkout 9a5a0ff6c3a0435e723e41f98fe82ec2ad19cf50 && \
    yes | sh build.sh && \
    ln -s /opt/Haploflow/build/haploflow /bin/haploflow

## Install uv
ENV HOME=/opt/uv-home
ENV UV_PROJECT_ENVIRONMENT=/opt/venv
ENV PATH="/opt/uv-home/.local/bin:/opt/venv/bin:${PATH}"

RUN apt-get install -qy tar git
RUN wget -q https://astral.sh/uv/install.sh -O /tmp/uv-install.sh
RUN sh /tmp/uv-install.sh

## Install dependencies for genetracks/drawsvg
RUN apt-get install -q -y libcairo2-dev

## Install just the dependencies of MiCall (for faster build times in development).
COPY pyproject.toml README.md uv.lock /opt/micall/
RUN uv sync --frozen --managed-python --project /opt/micall --extra basespace --no-editable

## Trigger matplotlib to build its font cache
RUN python -c 'import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot'

COPY . /opt/micall/

RUN uv sync --frozen --managed-python --project /opt/micall --extra basespace --no-editable
RUN micall make_blast_db

## Sometimes BaseSpace will crash if python tries to create threads. We prevent thread creation here:
ENV OPENBLAS_NUM_THREADS=1
ENV OMP_NUM_THREADS=1
ENV MKL_NUM_THREADS=1
ENV NUMEXPR_NUM_THREADS=1

WORKDIR /data
ENTRYPOINT ["micall", "analyze"]
