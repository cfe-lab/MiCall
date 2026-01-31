# Generate the Singularity container to run MiCall on Kive.
Bootstrap: docker
From: debian:bookworm-slim

%help
    MiCall maps all the reads from a sample against a set of reference
    sequences, then stitches all the reads into consensus sequences and
    coverage maps.

    This Singularity container can be run on Kive: http://cfe-lab.github.io/Kive

    Change Notes: Comprehensive updates to the contig stitcher,
    including bug fixes and visualization enhancements.

%labels
    MAINTAINER BC CfE in HIV/AIDS https://github.com/cfe-lab/MiCall
    KIVE_INPUTS sample_info_csv fastq1 fastq2 bad_cycles_csv
    KIVE_OUTPUTS g2p_csv g2p_summary_csv remap_counts_csv \
        remap_conseq_csv unmapped1_fastq unmapped2_fastq conseq_ins_csv \
        failed_csv cascade_csv nuc_csv amino_csv insertions_csv conseq_csv \
        conseq_all_csv concordance_csv concordance_seed_csv failed_align_csv \
        coverage_scores_csv coverage_maps_tar aligned_csv g2p_aligned_csv \
        genome_coverage_csv genome_coverage_svg genome_concordance_svg stitcher_plot_svg
    KIVE_THREADS 1
    KIVE_MEMORY 6000

%setup
    # Unneeded once Singularity creates parent dirs:
    # https://github.com/singularityware/singularity/issues/1549
    mkdir ${SINGULARITY_ROOTFS}/opt/micall

%files
    ## These files will be deleted after the install.
    . /opt/micall/

%post
    echo ===== Installing Prerequisites ===== >/dev/null
    apt-get update -q
    apt-get install -q -y unzip wget

    echo ===== Installing blast ===== >/dev/null
    apt-get install -q -y ncbi-blast+

    echo ===== Installing Rust and merge-mates ===== >/dev/null
    wget -qO rustup.sh https://sh.rustup.rs
    chmod +x /rustup.sh
    /rustup.sh -y -q
    . /root/.cargo/env
    rm rustup.sh
    cargo install --quiet --root / --git https://github.com/jeff-k/merge-mates.git --rev 2fec61363f645e2008a4adff553d098beae21469

    echo ===== Installing bowtie2 ===== >/dev/null
    wget -q -O bowtie2.zip https://github.com/BenLangmead/bowtie2/releases/download/v2.2.8/bowtie2-2.2.8-linux-x86_64.zip
    unzip bowtie2.zip -d /opt/
    ln -s /opt/bowtie2-2.2.8/ /opt/bowtie2
    rm bowtie2.zip

    echo ===== Installing IVA ===== >/dev/null
    apt-get install -q -y zlib1g-dev libncurses5-dev libncursesw5-dev
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

    echo ===== Install iva via uv ===== >/dev/null
    export HOME=/opt/uv-home
    export UV_PROJECT_ENVIRONMENT=/opt/venv
    export PATH="/opt/uv-home/.local/bin:/opt/venv/bin:${PATH}"

    wget -q https://astral.sh/uv/install.sh -O /tmp/uv-install.sh
    sh /tmp/uv-install.sh

    uv tool install --python=3.11 'git+https://github.com/cfe-lab/iva.git@v1.1.1'

    echo ===== Installing Python packages ===== >/dev/null
    # Install dependencies for genetracks/drawsvg
    apt-get install -q -y libcairo2-dev
    # Install micall main executable.
    uv sync --managed-python --project /opt/micall --extra basespace --no-editable
    micall make_blast_db
    # Also trigger matplotlib to build its font cache.
    python -c 'import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot'

    # Cleanup.
    rm -rf /opt/micall /root /usr/lib/gcc /opt/uv-home/.cache /etc/apt/ /var/apt /etc/dpkg /var/log /var/cache /var/lib/apt /var/lib/dpkg
    mkdir -p /root

%environment
    export PATH=/opt/uv-home/.local/bin:/opt/bowtie2:/bin:/usr/local/bin
    export LANG=en_US.UTF-8

%runscript
    micall micall_kive "$@"

%apphelp filter_quality
    Post-processing of short-read alignments.

%applabels filter_quality
    KIVE_INPUTS quality_csv
    KIVE_OUTPUTS bad_cycles_csv
    KIVE_THREADS 1
    KIVE_MEMORY 200

%apprun filter_quality
    micall filter_quality "$@"

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
    micall micall_kive_resistance "$@"

%apprun denovo
    micall micall_kive --denovo "$@"

%applabels denovo
    KIVE_INPUTS sample_info_csv fastq1 fastq2 bad_cycles_csv
    KIVE_OUTPUTS g2p_csv g2p_summary_csv remap_counts_csv \
        remap_conseq_csv unmapped1_fastq unmapped2_fastq conseq_ins_csv \
        failed_csv cascade_csv nuc_csv amino_csv insertions_csv conseq_csv \
        conseq_all_csv concordance_csv concordance_seed_csv failed_align_csv \
        coverage_scores_csv coverage_maps_tar aligned_csv g2p_aligned_csv \
        genome_coverage_csv genome_coverage_svg genome_concordance_svg \
        stitcher_plot_svg unstitched_cascade_csv unstitched_conseq_csv \
        unstitched_contigs_csv contigs_csv read_entropy_csv \
        conseq_region_csv
    KIVE_THREADS 2
    KIVE_MEMORY 6000

%apphelp denovo
    Standard pipeline with de novo assembly instead of mapping to reference
    sequences.
