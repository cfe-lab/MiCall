# Generate the Singularity container to run MiCall on Kive.
Bootstrap: docker
From: centos:7

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

    ## HCV genotyping database
    micall/utils/hcv_geno /opt/hcv_geno/

%post
    ## Prerequisites
    yum update -q -y

    yum groupinstall -q -y 'development tools'
    yum install -q -y epel-release
    yum install -q -y python34 python34-devel unzip wget fontconfig

    # Install blast and build the HCV genotype database
    #yum install blast
    #makeblastdb -in /opt/hcv.fasta -parse_seqids -dbtype nucl


    ## Miniconda (Python 2) (Don't use this)
    #wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O miniconda.sh
    #bash miniconda.sh -b -p /opt/miniconda

    ## bowtie2
    wget -q -O bowtie2.zip https://github.com/BenLangmead/bowtie2/releases/download/v2.2.8/bowtie2-2.2.8-linux-x86_64.zip
    unzip bowtie2.zip -d /opt/
    ln -s /opt/bowtie2-2.2.8/ /opt/bowtie2
    rm bowtie2.zip

    ## Python packages, plus trigger matplotlib to build its font cache
    wget -q https://bootstrap.pypa.io/get-pip.py
    python3 get-pip.py
    rm get-pip.py
    cd /opt
    pip install -r /opt/micall/requirements-basespace.txt
    ln -s /usr/bin/cutadapt /usr/bin/cutadapt-1.11
    python3 -c 'import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot'

    # IVA assembler
    pip install iva

    yum groupremove -q -y 'development tools'
    yum remove -q -y epel-release wget python34-devel unzip
    yum autoremove -q -y
    yum clean all

    pip install biopython

    rm -rf /var/cache/yum

    ## CAUTION! This changes the default python command to python3!
    ## This breaks many things, including yum!
    ## To switch back to python2, use this command:
    # sudo alternatives --set python /usr/bin/python2
    alternatives --install /usr/bin/python python /usr/bin/python2 50
    alternatives --install /usr/bin/python python /usr/bin/python3 60

    ## Savage assembler
    #export PATH="/opt/miniconda/bin:$PATH"
    #source /opt/miniconda/bin/activate
    #conda config --add channels r
    #conda config --add channels defaults
    #conda config --add channels conda-forge
    #conda config --add channels bioconda
    #conda install savage
    #ls /opt/miniconda/bin/

%environment
    export PATH=/bin:/opt/bowtie2
