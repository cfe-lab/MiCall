# Generate the Singularity container to run MiCall on Kive.
Bootstrap: docker
From: centos:7

%help
    MiCall maps all the reads from a sample against a set of reference
    sequences, then stitches all the reads into consensus sequences and
    coverage maps.

    This Singularity container can be run on Kive: http://cfe-lab.github.io/Kive

    Change Notes: SARS-CoV-2 added, HIVdb upgraded to 8.9-1, filter primer
    dimers.

%labels
    MAINTAINER BC CfE in HIV/AIDS https://github.com/cfe-lab/MiCall
    KIVE_INPUTS sample_info_csv fastq1 fastq2 bad_cycles_csv
    KIVE_OUTPUTS g2p_csv g2p_summary_csv remap_counts_csv \
        remap_conseq_csv unmapped1_fastq unmapped2_fastq insertions_csv \
        failed_csv cascade_csv nuc_csv amino_csv conseq_csv \
        conseq_all_csv conseq_region_csv conseq_stitched_csv failed_align_csv \
        coverage_scores_csv coverage_maps_tar aligned_csv g2p_aligned_csv \
        genome_coverage_csv genome_coverage_svg
    KIVE_THREADS 1
    KIVE_MEMORY 6000

%setup
    # Unneeded once Singularity creates parent dirs:
    # https://github.com/singularityware/singularity/issues/1549
    mkdir ${SINGULARITY_ROOTFS}/opt/micall
    mkdir ${SINGULARITY_ROOTFS}/opt/micall/micall

%files
    ## MiCall
    micall_docker.py /opt/micall/
    micall_kive.py /opt/micall/
    micall_kive_resistance.py /opt/micall/
    micall/__init__.py /opt/micall/micall/
    micall/project* /opt/micall/micall/

    micall/core    /opt/micall/micall/core
    micall/data    /opt/micall/micall/data
    micall/drivers    /opt/micall/micall/drivers
    micall/g2p     /opt/micall/micall/g2p
    micall/resistance   /opt/micall/micall/resistance
    micall/monitor /opt/micall/micall/monitor
    micall/utils   /opt/micall/micall/utils

    ## Gotoh
    micall/alignment /opt/micall/alignment
    requirements.txt /opt/micall/
    requirements-basespace.txt /opt/micall/

    ## HCV genotyping database
    micall/blast_db /opt/micall/micall/blast_db

%post
    echo ===== Installing Prerequisites ===== >/dev/null
    yum update -q -y

    yum groupinstall -q -y 'development tools'
    yum install -q -y epel-release
    yum install -q -y unzip wget fontconfig bzip2-devel xz-devel openssl-devel \
        libffi-devel sqlite-devel

    echo ===== Installing Python ===== >/dev/null
    wget -q https://www.python.org/ftp/python/3.8.3/Python-3.8.3.tar.xz
    tar xJf Python*
    rm Python*.xz
    cd Python*
    ./configure --enable-optimizations
    make altinstall
    cd ..
    rm -rf Python*
    ln -s /usr/local/bin/python3.8 /usr/local/bin/python3

    echo ===== Installing blast ===== >/dev/null
    cd /root
    # Saved our own copy, because download was slow from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-1.x86_64.rpm
    wget -q https://github.com/cfe-lab/MiCall/releases/download/v7.12.dev28/ncbi-blast-2.6.0+-1.x86_64.rpm
    yum install -q -y ncbi-blast-2.6.0+-1.x86_64.rpm
    rm ncbi-blast-2.6.0+-1.x86_64.rpm
    python3 /opt/micall/micall/blast_db/make_blast_db.py

    echo ===== Installing Rust and merge-mates ===== >/dev/null
    yum install -q -y rust cargo
    cargo install --quiet --root / --git https://github.com/jeff-k/merge-mates.git

    ## Miniconda (Python 2) (Don't use this)
    #wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O miniconda.sh
    #bash miniconda.sh -b -p /opt/miniconda

    echo ===== Installing bowtie2 ===== >/dev/null
    wget -q -O bowtie2.zip https://github.com/BenLangmead/bowtie2/releases/download/v2.2.8/bowtie2-2.2.8-linux-x86_64.zip
    unzip -qq bowtie2.zip -d /opt/
    ln -s /opt/bowtie2-2.2.8/ /opt/bowtie2
    rm bowtie2.zip

    echo ===== Installing IVA dependencies ===== >/dev/null
    yum install -q -y tcsh ncurses-devel zlib-devel
    cd /bin
    wget -q http://sun.aei.polsl.pl/kmc/download-2.1.1/linux/kmc
    wget -q http://sun.aei.polsl.pl/kmc/download-2.1.1/linux/kmc_dump
    chmod +x kmc kmc_dump
    cd /opt
    wget -q https://sourceforge.net/projects/mummer/files/mummer/3.23/MUMmer3.23.tar.gz
    tar -xzf MUMmer3.23.tar.gz --no-same-owner
    cd MUMmer3.23
    make --quiet install
    rm -r docs src ../MUMmer3.23.tar.gz
    ln -s /opt/MUMmer3.23/nucmer \
        /opt/MUMmer3.23/delta-filter \
        /opt/MUMmer3.23/show-coords \
        /bin
    cd /opt
    wget -q https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
    tar -xf samtools-1.3.1.tar.bz2 --no-same-owner --bzip2
    cd samtools-1.3.1
    ./configure --quiet --prefix=/
    make --quiet
    make --quiet install
    cd /opt
    rm -rf samtools-1.3.1*
    wget -q http://downloads.sourceforge.net/project/smalt/smalt-0.7.6-bin.tar.gz
    tar -xzf smalt-0.7.6-bin.tar.gz --no-same-owner
    ln -s /opt/smalt-0.7.6-bin/smalt_x86_64 /bin/smalt

    echo ===== Installing Python packages ===== >/dev/null
    # Also trigger matplotlib to build its font cache.
    wget -q https://bootstrap.pypa.io/get-pip.py
    python3 get-pip.py
    rm get-pip.py
    cd /opt
    pip install --quiet -r /opt/micall/requirements.txt
    ln -s /usr/local/bin/cutadapt /usr/local/bin/cutadapt-1.11
    python3 -c 'import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot'

    # Install dependencies for genetracks/drawSvg
    yum install -q -y cairo-devel cairo cairo-tools zlib-devel

    yum groupremove -q -y 'development tools'
    yum remove -q -y epel-release wget unzip
    yum autoremove -q -y
    yum clean all

    rm -rf /var/cache/yum

    ## CAUTION! This changes the default python command to python3!
    ## This breaks many things, including yum!
    ## To switch back to python2, use this command:
    # sudo alternatives --set python /usr/bin/python2
    alternatives --install /usr/bin/python python /usr/bin/python2 50
    alternatives --install /usr/bin/python python /usr/local/bin/python3 60

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
    KIVE_INPUTS main_amino_csv midi_amino_csv main_nuc_csv
    KIVE_OUTPUTS resistance_csv mutations_csv nuc_mutations_csv \
        resistance_fail_csv resistance_pdf resistance_consensus_csv
    KIVE_THREADS 1
    KIVE_MEMORY 200

%apprun resistance
    python /opt/micall/micall_kive_resistance.py "$@"

%apprun denovo
    python /opt/micall/micall_kive.py --denovo "$@"

%applabels denovo
    KIVE_INPUTS sample_info_csv fastq1 fastq2 bad_cycles_csv
    KIVE_OUTPUTS g2p_csv g2p_summary_csv remap_counts_csv \
        remap_conseq_csv unmapped1_fastq unmapped2_fastq insertions_csv \
        failed_csv cascade_csv nuc_csv amino_csv conseq_csv \
        conseq_all_csv conseq_region_csv conseq_stitched_csv failed_align_csv \
        coverage_scores_csv coverage_maps_tar aligned_csv g2p_aligned_csv \
        genome_coverage_csv genome_coverage_svg contigs_csv read_entropy_csv
    KIVE_THREADS 2
    KIVE_MEMORY 6000

%apphelp denovo
    Standard pipeline with de novo assembly instead of mapping to reference
    sequences.
