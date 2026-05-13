#!/usr/bin/env python3
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from pathlib import Path
from typing import Sequence
import subprocess
import sys

from micall.utils.docker_build import build


SINGULARITY_TEMPLATE = """

Bootstrap: docker-archive
From: ./simgs/micall-{container_sha}.tar

%help
    MiCall maps all the reads from a sample against a set of reference
    sequences, then stitches all the reads into consensus sequences and
    coverage maps.

    This Singularity container can be run on Kive: http://cfe-lab.github.io/Kive

%labels
    MAINTAINER BC CfE in HIV/AIDS https://github.com/cfe-lab/MiCall
    KIVE_INPUTS sample_info_csv fastq1 fastq2 bad_cycles_csv
    KIVE_OUTPUTS g2p_csv g2p_summary_csv remap_counts_csv \\
        remap_conseq_csv unmapped1_fastq unmapped2_fastq conseq_ins_csv \\
        failed_csv cascade_csv nuc_csv amino_csv insertions_csv conseq_csv \\
        conseq_all_csv concordance_csv concordance_seed_csv failed_align_csv \\
        coverage_scores_csv coverage_maps_tar aligned_csv g2p_aligned_csv \\
        genome_coverage_csv genome_coverage_svg genome_concordance_svg
    KIVE_THREADS 1
    KIVE_MEMORY 6000

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
    KIVE_OUTPUTS resistance_csv mutations_csv nuc_mutations_csv \\
        resistance_fail_csv resistance_pdf resistance_consensus_csv
    KIVE_THREADS 1
    KIVE_MEMORY 200

%apprun resistance
    micall micall_kive_resistance "$@"

%apprun denovo
    micall micall_kive --denovo "$@"

%applabels denovo
    KIVE_INPUTS sample_info_csv fastq1 fastq2 bad_cycles_csv
    KIVE_OUTPUTS g2p_csv g2p_summary_csv remap_counts_csv \\
        remap_conseq_csv unmapped1_fastq unmapped2_fastq conseq_ins_csv \\
        failed_csv cascade_csv nuc_csv amino_csv insertions_csv conseq_csv \\
        conseq_all_csv concordance_csv concordance_seed_csv failed_align_csv \\
        coverage_scores_csv coverage_maps_tar aligned_csv g2p_aligned_csv \\
        genome_coverage_csv genome_coverage_svg genome_concordance_svg \\
        unstitched_cascade_csv unstitched_conseq_csv unstitched_contigs_csv \\
        contigs_csv stitcher_plot_svg read_entropy_csv \\
        conseq_region_csv
    KIVE_THREADS 2
    KIVE_MEMORY 6000

%apphelp denovo
    Standard pipeline with de novo assembly instead of mapping to reference
    sequences.

"""


def get_container_sha(repository_name: str) -> str:
    """Return the short content-addressable SHA of a docker image."""
    full_id = subprocess.check_output(
        ['docker', 'inspect', '--format={{.Id}}', '--', repository_name],
        text=True,
    ).strip()
    # full_id is typically "sha256:<hex>" — strip the algorithm prefix.
    if ':' in full_id:
        full_id = full_id.split(':', 1)[1]
    return full_id[:12]


def save_docker_archive(repository_name: str, container_sha: str) -> Path:
    """Save the docker image as a tar archive under ./simgs/ and return the path."""
    simgs_dir = Path('simgs')
    simgs_dir.mkdir(exist_ok=True)
    archive_path = simgs_dir / f'micall-{container_sha}.tar'
    print(f'Saving docker image to {archive_path} ...')
    subprocess.check_call(
        ['docker', 'save', '--output', str(archive_path), '--', repository_name],
    )
    return archive_path


def generate_singularity_def(container_sha: str) -> str:
    """Return the content of a Singularity definition file for the given image."""
    return SINGULARITY_TEMPLATE.format(container_sha=container_sha)


def get_parser() -> ArgumentParser:
    parser = ArgumentParser(
        description=(
            'Build a Docker image, save it as a tar archive, and generate '
            'a Singularity definition file that references that archive.'
        ),
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        '--output',
        default='Singularity.def',
        help='Path to write the generated Singularity definition file.',
    )
    return parser


def main(argv: Sequence[str]) -> int:
    parser = get_parser()
    args = parser.parse_args(argv)

    repository_name = build()
    container_sha = get_container_sha(repository_name)
    save_docker_archive(repository_name, container_sha)

    def_content = generate_singularity_def(container_sha)
    output_path = Path(args.output)
    output_path.write_text(def_content)
    print(f'Singularity definition written to {output_path}')

    return 0


def entry() -> None:
    sys.exit(main(sys.argv[1:]))


if __name__ == '__main__':
    entry()
