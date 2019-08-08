# Generate the Singularity container to run MiCall on Kive.
Bootstrap: docker
From: centos:7

%help
    MiCall maps all the reads from a sample against a set of reference
    sequences, then stitches all the reads into consensus sequences and
    coverage maps.

    This Singularity container can be run on Kive: http://cfe-lab.github.io/Kive

    Change Notes: Includes denovo assembly apps.

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
    micall/data    /opt/micall/micall/data/
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
    micall/blast_db /opt/micall/micall/blast_db

    ## H77 reference
    micall/utils/hcv_geno /opt/micall/micall/utils/hcv_geno

%post
    echo ===== Installing Prerequisites ===== >/dev/null
    yum update -q -y

    yum groupinstall -q -y 'development tools'
    yum install -q -y epel-release
    yum install -q -y https://centos7.iuscommunity.org/ius-release.rpm
    yum install -q -y python36 python36-devel unzip wget fontconfig

    echo ===== Installing Rust and merge-mates ===== >/dev/null
    yum install -q -y rust cargo
    cargo install --root / --git https://github.com/jeff-k/merge-mates.git

    echo ===== Installing Savage ===== >/dev/null
    yum install -q -y zlib-devel boost-timer boost-program-options boost-devel
    cargo install --root / --git https://github.com/jbaaijens/rust-overlaps.git --tag v0.1.1

    wget -q https://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.17.tar.bz2
    tar xjf bwa-0.7.17.tar.bz2 --no-same-owner
    cd bwa-0.7.17
    make
    mv bwa /bin
    cd ..
    rm -rf bwa-0.7.17 bwa-0.7.17.tar.bz2

    wget -q https://github.com/pachterlab/kallisto/releases/download/v0.44.0/kallisto_linux-v0.44.0.tar.gz
    tar xzf kallisto_linux-v0.44.0.tar.gz --no-same-owner
    mv kallisto_linux-v0.44.0/kallisto /bin
    rm -rf kallisto_linux-v0.44.0 kallisto_linux-v0.44.0.tar.gz

    cd /opt
    git clone https://bitbucket.org/jbaaijens/savage.git
    cd savage
    git checkout tags/v0.4.0
    make
    echo \#\!/usr/bin/env sh > /bin/savage
    echo /opt/savage/savage.py \$@ >> /bin/savage
    chmod +x /bin/savage
    cd /

    yum remove -q -y zlib-devel boost-devel

    echo ===== Installing blast ===== >/dev/null
    cd /root
    wget -q ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-1.x86_64.rpm
    yum install -q -y ncbi-blast-2.6.0+-1.x86_64.rpm
    rm ncbi-blast-2.6.0+-1.x86_64.rpm
    python3.6 /opt/micall/micall/blast_db/make_blast_db.py


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
    python3.6 get-pip.py
    rm get-pip.py
    cd /opt
    pip install -r /opt/micall/requirements-basespace.txt
    ln -s /usr/local/bin/cutadapt /usr/local/bin/cutadapt-1.11
    python3 -c 'import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot'

    # build h77 reference
    cd /opt/micall/micall/utils/hcv_geno
    bwa index hxb2.fasta
    bwa index h77.fasta
    bwa index hcv.fasta

    # Install genetracks
    yum install -q -y cairo-devel cairo cairo-tools zlib-devel
    pip install cairosvg
    pip install cairocffi
    pip install pysam
    cd /opt
    git clone https://github.com/jeff-k/genetracks.git
    cd genetracks
    pip install .

    yum groupremove -q -y 'development tools'
    yum remove -q -y epel-release wget python36-devel unzip
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

%apprun denovo_hcv
    sh /opt/micall/micall/utils/iva.sh "$@"

%applabels denovo_hcv
    KIVE_INPUTS wg1 wg2 mid1 mid2 sample_name
    KIVE_OUTPUTS wg_fasta mid_fasta alignment_svg subtyping_csv alignment_png

%apphelp denovo_hcv
    Two-part HCV denovo assembly

%apprun proviral_hiv
    sh /opt/micall/micall/utils/proviral.sh "$@"

%applabels proviral_hiv
    KIVE_INPUTS r1 r2 sample_name
    KIVE_OUTPUTS alignment_svg summary_csv assembly_fasta alignment_png

%apphelp proviral_hiv
    Bespoke pipeline for highly mutated HIV proviruses

%apprun denovo
    python /opt/micall/micall_kive.py --denovo "$@"

%applabels denovo
    KIVE_INPUTS fastq1 fastq2 bad_cycles_csv
    KIVE_OUTPUTS g2p_csv g2p_summary_csv remap_counts_csv \
        remap_conseq_csv unmapped1_fastq unmapped2_fastq conseq_ins_csv \
        failed_csv cascade_csv nuc_csv amino_csv coord_ins_csv conseq_csv \
        conseq_region_csv failed_align_csv coverage_scores_csv \
        coverage_maps_tar aligned_csv g2p_aligned_csv contigs_csv \
        contigs_coverage_csv contigs_coverage_svg
    KIVE_THREADS 2
    KIVE_MEMORY 6000

%apphelp denovo
    Standard pipeline with de novo assembly instead of mapping to reference
    sequences.
