# Generate the Singularity container to run MiCall on Kive.
Bootstrap: docker
From: donaim/micallsys:7.17.0-1694-g0fc07d356

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
    pip install /opt/micall
    micall make_blast_db
    python -c 'import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot'

    # Cleanup.
    rm -rf /opt/micall /root /usr/lib/gcc /opt/uv-home/.cache /etc/apt/ /var/apt /etc/dpkg /var/log /var/cache /var/lib/apt /var/lib/dpkg
    mkdir -p /root

%environment
    export PATH=/opt/bowtie2:/bin:/usr/local/bin
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
        conseq_region_csv conseq_stitched_csv
    KIVE_THREADS 2
    KIVE_MEMORY 6000

%apphelp denovo
    Standard pipeline with de novo assembly instead of mapping to reference
    sequences.
