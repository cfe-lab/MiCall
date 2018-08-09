# Generate the Singularity container to run MiCall on Kive.
Bootstrap: docker
From: python:3.4

%help
    MiCall maps all the reads from a sample against a set of reference
    sequences, then stitches all the reads into consensus sequences and
    coverage maps.

    This Singularity container can be run on Kive: http://cfe-lab.github.io/Kive

%label
    MAINTAINER BC CfE in HIV/AIDS https://github.com/cfe-lab/MiCall

%setup
    # Unneeded once Singularity creates parent dirs:
    # https://github.com/singularityware/singularity/issues/1549
    mkdir ${SINGULARITY_ROOTFS}/opt/micall
    mkdir ${SINGULARITY_ROOTFS}/opt/micall/micall

%files
    ## MiCall
    micall_basespace.py /opt/micall/
    micall_kive.py /opt/micall/
    micall_kive_resistance.py /opt/micall/
    micall/__init__.py /opt/micall/micall/
    micall/project* /opt/micall/micall/

    micall/core    /opt/micall/micall/core/
    micall/drivers    /opt/micall/micall/drivers/
    micall/g2p     /opt/micall/micall/g2p/
    micall/resistance   /opt/micall/micall/resistance/
    micall/monitor /opt/micall/micall/monitor/
    micall/utils   /opt/micall/micall/utils/

    ## Gotoh
    micall/alignment /opt/micall/alignment
    requirements.txt /opt/micall/
    requirements-basespace.txt /opt/micall/

%post
    ## Prerequisites
    apt-get update -qq --fix-missing
    apt-get install -qq -y unzip wget
    rm -rf /var/lib/apt/lists/*

    ## bowtie2
    wget -q -O bowtie2.zip https://github.com/BenLangmead/bowtie2/releases/download/v2.2.8/bowtie2-2.2.8-linux-x86_64.zip
    unzip bowtie2.zip -d /opt/
    ln -s /opt/bowtie2-2.2.8/ /opt/bowtie2
    rm bowtie2.zip

    ## Python packages, plus trigger matplotlib to build its font cache
    cd /opt
    pip install -r /opt/micall/requirements-basespace.txt
    ln -s /usr/local/bin/cutadapt /usr/local/bin/cutadapt-1.11
    /usr/local/bin/python -c 'import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot'

%environment
    export PATH=$PATH:/opt/bowtie2