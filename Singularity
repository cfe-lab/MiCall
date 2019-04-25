# Generate the Singularity container to run MiCall on Kive.
Bootstrap: docker
From: centos:7

%help
    MiCall maps all the reads from a sample against a set of reference
    sequences, then stitches all the reads into consensus sequences and
    coverage maps.

    This Singularity container can be run on Kive: http://cfe-lab.github.io/Kive

%labels
    MAINTAINER BC CfE in HIV/AIDS https://github.com/cfe-lab/MiCall
    KIVE_INPUTS fastq1 fastq2 bad_cycles_csv
    KIVE_OUTPUTS g2p_csv g2p_summary_csv remap_counts_csv \
        remap_conseq_csv unmapped1_fastq unmapped2_fastq conseq_ins_csv \
        failed_csv cascade_csv nuc_csv amino_csv coord_ins_csv conseq_csv \
        conseq_region_csv failed_align_csv coverage_scores_csv \
        coverage_maps_tar aligned_csv g2p_aligned_csv
    KIVE_THREADS 1
    KIVE_MEMORY 6000

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
    yum update -q -y

    yum groupinstall -q -y 'development tools'
    yum install -q -y epel-release
    yum install -q -y python36 python36-devel unzip wget fontconfig

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
    ln -s /usr/local/bin/cutadapt /usr/local/bin/cutadapt-1.11
    python3 -c 'import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot'

    yum groupremove -q -y 'development tools'
    yum remove -q -y epel-release wget python34-devel unzip
    yum autoremove -q -y
    yum clean all

    rm -rf /var/cache/yum

    ## CAUTION! This changes the default python command to python3!
    ## This breaks many things, including yum!
    ## To switch back to python2, use this command:
    # sudo alternatives --set python /usr/bin/python2
    alternatives --install /usr/bin/python python /usr/bin/python2 50
    alternatives --install /usr/bin/python python /usr/bin/python3 60

%environment
    export PATH=/opt/bowtie2:/bin:/usr/local/bin
    export LANG=en_US.UTF-8

%runscript
    python /opt/micall/micall_kive.py "$@"

%apphelp filter_quality
    Post-processing of short-read alignments.

%applabels filter_quality
    KIVE_INPUTS quality_csv
    KIVE_OUTPUTS bad_cycles_csv
    KIVE_THREADS 1
    KIVE_MEMORY 200

%apprun filter_quality
    PYTHONPATH=/opt/micall python -m micall.core.filter_quality "$@"

%apphelp resistance
    Combine HCV results with HCV-Midi results, and generate resistance
    interpretation.

%applabels resistance
    KIVE_INPUTS main_amino_csv midi_amino_csv
    KIVE_OUTPUTS resistance_csv mutations_csv resistance_fail_csv \
        resistance_pdf resistance_consensus_csv
    KIVE_THREADS 1
    KIVE_MEMORY 200

%apprun resistance
    python /opt/micall/micall_kive_resistance.py "$@"
