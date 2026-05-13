#!/usr/bin/env python3
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from pathlib import Path
from typing import Sequence
import logging
import subprocess
import sys
import urllib3

from micall.utils.docker_build import build, get_latest_git_tag


logger = logging.getLogger(__name__)

urllib3.connectionpool.log.setLevel(logging.ERROR)


SINGULARITY_TEMPLATE = """

Bootstrap: docker-archive
From: ./simgs/micall-{container_sha}.tar

%help
    {description}

    This Singularity container can be run on Kive: http://cfe-lab.github.io/Kive

%labels
    MAINTAINER BC CfE in HIV/AIDS https://github.com/cfe-lab/MiCall
    KIVE_INPUTS sample_info_csv fastq1 fastq2 bad_cycles_csv
    KIVE_OUTPUTS g2p_csv g2p_summary_csv remap_counts_csv \
        remap_conseq_csv unmapped1_fastq unmapped2_fastq conseq_ins_csv \
        failed_csv cascade_csv nuc_csv amino_csv insertions_csv conseq_csv \
        conseq_all_csv concordance_csv concordance_seed_csv failed_align_csv \
        coverage_scores_csv coverage_maps_tar aligned_csv g2p_aligned_csv \
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
        unstitched_cascade_csv unstitched_conseq_csv unstitched_contigs_csv \
        contigs_csv stitcher_plot_svg read_entropy_csv \
        conseq_region_csv
    KIVE_THREADS 2
    KIVE_MEMORY 6000

%apphelp denovo
    Standard pipeline with de novo assembly instead of mapping to reference
    sequences.
"""


DESCRIPTION = """
    MiCall maps all the reads from a sample against a set of reference
    sequences, then stitches all the reads into consensus sequences and
    coverage maps.
"""


SINGULARITY_IMAGE_DIR = Path('simgs')
SINGULARITY_DEFINITION_PATH = SINGULARITY_IMAGE_DIR / 'Singularity.def'


def get_parser() -> ArgumentParser:
    parser = ArgumentParser(
        description=(
            'Build a Docker image, save it as a tar archive, generate a '
            'Singularity definition file, and build the Singularity image.'
        ),
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    verbosity_group = parser.add_mutually_exclusive_group()
    verbosity_group.add_argument('--verbose', action='store_true', help='Increase output verbosity.')
    verbosity_group.add_argument('--no-verbose', action='store_true', help='Normal output verbosity.', default=True)
    verbosity_group.add_argument('--debug', action='store_true', help='Maximum output verbosity.')
    verbosity_group.add_argument('--quiet', action='store_true', help='Minimize output verbosity.')

    subparsers = parser.add_subparsers(dest='mode')

    push_parser = subparsers.add_parser(
        'push',
        help='Upload the built Singularity image to Kive.',
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    push_parser.add_argument(
        '--family',
        help='Kive container family name or ID for uploading the built image.',
    )
    push_parser.add_argument(
        '--users',
        nargs='*',
        help='Users to grant access to the uploaded container.',
    )
    push_parser.add_argument(
        '--groups',
        nargs='*',
        default=['Everyone'],
        help='Groups to grant access to the uploaded container.',
    )

    subparsers.add_parser(
        'nopush',
        help='Build artifacts but do not upload to Kive.',
        formatter_class=ArgumentDefaultsHelpFormatter,
    )

    # Default mode is nopush when no subcommand is provided.
    parser.set_defaults(mode='nopush')
    return parser


def configure_logging(args) -> None:
    if args.quiet:
        logger.setLevel(logging.ERROR)
    elif args.verbose:
        logger.setLevel(logging.INFO)
    elif args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.WARN)

    logging.basicConfig(
        level=logger.level,
        format='%(asctime)s[%(levelname)s]%(name)s: %(message)s',
    )


def get_container_sha(repository_name: str) -> str:
    """Return the short content-addressable SHA of a docker image."""
    logger.debug('Inspecting Docker image %s for its immutable ID.', repository_name)
    full_id = subprocess.check_output(
        ['docker', 'inspect', '--format={{.Id}}', '--', repository_name],
        text=True,
    ).strip()
    logger.debug('Docker image %s resolved to raw ID %s.', repository_name, full_id)
    if ':' in full_id:
        full_id = full_id.split(':', 1)[1]
    container_sha = full_id[:12]
    logger.info('Docker image %s will be archived as micall-%s.tar.', repository_name, container_sha)
    return container_sha


def save_docker_archive(repository_name: str, container_sha: str) -> Path:
    """Save the docker image as a tar archive under ./simgs/ and return the path."""
    SINGULARITY_IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    archive_path = SINGULARITY_IMAGE_DIR / f'micall-{container_sha}.tar'
    if archive_path.exists():
        logger.info('Skipping docker save because archive already exists at %s.', archive_path)
        return archive_path

    logger.info('Saving Docker image %s to %s.', repository_name, archive_path)
    subprocess.check_call(
        ['docker', 'save', '--output', str(archive_path), '--', repository_name],
    )
    logger.debug('Docker archive %s created successfully.', archive_path)
    return archive_path


def generate_singularity_def(container_sha: str) -> str:
    """Return the content of a Singularity definition file for the given image."""
    logger.debug('Rendering Singularity definition for archive micall-%s.tar.', container_sha)
    return SINGULARITY_TEMPLATE.format(container_sha=container_sha, description=DESCRIPTION)


def write_singularity_definition(container_sha: str) -> Path:
    definition_path = SINGULARITY_DEFINITION_PATH
    definition_content = generate_singularity_def(container_sha)
    definition_path.write_text(definition_content)
    logger.info('Wrote Singularity definition to %s.', definition_path)
    return definition_path


def build_singularity_image(definition_path: Path, container_sha: str) -> Path:
    SINGULARITY_IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    image_path = SINGULARITY_IMAGE_DIR / f'micall-{container_sha}.sif'
    latest_link_path = SINGULARITY_IMAGE_DIR / 'micall-latest.sif'
    if image_path.exists():
        logger.info('Skipping Singularity build because image already exists at %s.', image_path)
    else:
        logger.info('Building Singularity image %s from %s.', image_path, definition_path)
        subprocess.check_call(
            ['singularity', 'build', str(image_path), str(definition_path)],
        )

    if latest_link_path.exists() or latest_link_path.is_symlink():
        latest_link_path.unlink()
    # Keep link target relative to simgs/ so the tree is portable.
    latest_link_path.symlink_to(image_path.name)
    logger.info('Updated latest image symlink %s -> %s.', latest_link_path, image_path.name)
    logger.debug('Singularity image %s built successfully.', image_path)
    return image_path


def push_image_with_kivecli(
    image_path: Path,
    repository_name: str,
    family_name_or_id: str,
    tag: str,
    users: list[str] | None,
    groups: list[str] | None,
) -> int:
    logger.info('Preparing to push Singularity image %s with kivecli.', image_path)
    logger.debug('Source Docker repository is %s.', repository_name)
    logger.debug('Uploading to Kive family=%s tag=%s users=%s groups=%s.',
                 family_name_or_id,
                 tag,
                 users,
                 groups)

    if not users and not groups:
        raise ValueError('At least one of --users or --groups must be provided for kivecli upload.')

    try:
        import kivecli.makecontainer as makecontainer
        import kivecli.logger
    except ImportError as ex:
        raise ImportError(
            'kivecli is required for upload. Install dev dependencies that include kivecli.'
        ) from ex

    kivecli.logger.logger.setLevel(logger.level)
    logger.info('Starting kivecli makecontainer upload for %s.', image_path)
    description = '\n'.join(line.strip() for line in DESCRIPTION.strip().splitlines())
    is_json = False
    result_code = makecontainer.main_typed(
        image_path=image_path,
        family_name_or_id=family_name_or_id,
        tag=tag,
        description=description,
        users=users,
        groups=groups,
        is_json=is_json,
    )
    logger.info('kivecli upload completed with status code %s.', result_code)
    return result_code


def main(argv: Sequence[str]) -> int:
    parser = get_parser()
    args = parser.parse_args(argv)
    configure_logging(args)

    logger.info('Starting Singularity build workflow.')
    logger.debug('Command-line arguments: %s', args)

    logger.info('Building Docker image first.')
    repository_name = build()
    logger.info('Docker build completed: %s', repository_name)

    container_sha = get_container_sha(repository_name)
    archive_path = save_docker_archive(repository_name, container_sha)
    definition_path = write_singularity_definition(container_sha)
    image_path = build_singularity_image(definition_path, container_sha)
    kive_tag = get_latest_git_tag()

    logger.info('Singularity archive ready at %s.', archive_path)
    logger.info('Singularity image ready at %s.', image_path)

    if args.mode == 'nopush':
        logger.info('Skipping upload because mode is nopush (default).')
        return 0

    logger.info('Using Kive tag %s for uploaded container.', kive_tag)
    return push_image_with_kivecli(
        image_path=image_path,
        repository_name=repository_name,
        family_name_or_id=args.family,
        tag=kive_tag,
        users=args.users,
        groups=args.groups,
    )


def entry() -> None:
    sys.exit(main(sys.argv[1:]))


if __name__ == '__main__':
    entry()
