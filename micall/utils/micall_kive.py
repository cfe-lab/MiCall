#! /usr/bin/env python

"""
Entry script that serves as an entry point of MiCall's Singularity image.
This file is run by Kive.
"""

import logging
import shutil
import tarfile
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import os

from micall.drivers.sample import Sample
from micall.g2p.pssm_lib import Pssm

logger = logging.getLogger(__name__)


def parse_args():
    parser = ArgumentParser(description='Map FASTQ files to references.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    # inputs
    parser.add_argument('sample_info_csv',
                        help='sample name and project code')
    parser.add_argument('fastq1',
                        help='FASTQ containing forward reads')
    parser.add_argument('fastq2',
                        help='FASTQ containing reverse reads')
    parser.add_argument('bad_cycles_csv',
                        help='list of tiles and cycles rejected for poor quality')

    # outputs
    parser.add_argument('g2p_csv',
                        help='CSV containing g2p predictions.')
    parser.add_argument('g2p_summary_csv',
                        help='CSV containing overall call for the sample.')
    parser.add_argument('remap_counts_csv',
                        help='CSV containing numbers of mapped reads')
    parser.add_argument('remap_conseq_csv',
                        help='CSV containing mapping consensus sequences')
    parser.add_argument('unmapped1_fastq',
                        help='FASTQ R1 of reads that failed to map to any region')
    parser.add_argument('unmapped2_fastq',
                        help='FASTQ R2 of reads that failed to map to any region')
    parser.add_argument('conseq_ins_csv',
                        help='CSV containing insertions relative to sample consensus')
    parser.add_argument('failed_csv',
                        help='CSV containing reads that failed to merge')
    parser.add_argument('cascade_csv',
                        help='count of reads at each step')
    parser.add_argument('nuc_csv',
                        help='CSV containing nucleotide frequencies')
    parser.add_argument('amino_csv',
                        help='CSV containing amino frequencies')
    parser.add_argument('insertions_csv',
                        help='CSV containing all insertions')
    parser.add_argument('conseq_csv',
                        help='CSV containing consensus sequences')
    parser.add_argument('conseq_all_csv',
                        help='CSV containing consensus sequences with low coverage')
    parser.add_argument('concordance_csv',
                        help='CSV containing coordinate reference concordance measures for each region')
    parser.add_argument('concordance_seed_csv',
                        help='CSV containing seed concordance measures for each region')
    parser.add_argument('failed_align_csv',
                        help='CSV containing any consensus that failed to align')
    parser.add_argument('coverage_scores_csv',
                        help='CSV coverage scores.')
    parser.add_argument('coverage_maps_tar',
                        help='tar file of coverage maps.')
    parser.add_argument('aligned_csv',
                        help='CSV containing individual reads aligned to consensus')
    parser.add_argument('g2p_aligned_csv',
                        help='CSV containing individual reads aligned to V3LOOP')
    parser.add_argument('genome_coverage_csv',
                        nargs='?',
                        help='CSV of coverage levels in full-genome coordinates')
    parser.add_argument('genome_coverage_svg',
                        nargs='?',
                        help='SVG diagram of coverage in full-genome coordinates')
    parser.add_argument('genome_concordance_svg',
                        nargs='?',
                        help='SVG diagram of concordance in full-genome coordinates')
    parser.add_argument('--denovo',
                        action='store_true',
                        help='Use de novo assembly instead of mapping to '
                             'reference sequences.')
    parser.add_argument('unstitched_cascade_csv',
                        nargs='?',
                        help='count of reads at each step')
    parser.add_argument('unstitched_conseq_csv',
                        nargs='?',
                        help='CSV containing mapping unstitched consensus sequences')
    parser.add_argument('unstitched_contigs_csv',
                        nargs='?',
                        help='CSV containing contigs built by de novo assembly')
    parser.add_argument('contigs_csv',
                        nargs='?',
                        help='CSV containing contigs built by de novo assembly and stitched by our stitcher')
    parser.add_argument('stitcher_plot_svg',
                        nargs='?',
                        help='SVG plot showing the contig stitcher coverage and operations')
    parser.add_argument('read_entropy_csv',
                        nargs='?',
                        help='CSV containing read pair length counts')
    parser.add_argument('conseq_region_csv',
                        nargs='?',
                        help='CSV containing consensus sequences, split by region')

    return parser.parse_args()


def load_sample(args):
    """ Load the data from Kive's command-line arguments. """
    scratch_path = os.path.join(os.path.dirname(args.cascade_csv), 'scratch')
    shutil.rmtree(scratch_path, ignore_errors=True)

    sample = Sample(sample_info_csv=args.sample_info_csv,
                    fastq1=args.fastq1,
                    fastq2=args.fastq2,
                    bad_cycles_csv=args.bad_cycles_csv,
                    g2p_csv=args.g2p_csv,
                    g2p_summary_csv=args.g2p_summary_csv,
                    remap_counts_csv=args.remap_counts_csv,
                    remap_conseq_csv=args.remap_conseq_csv,
                    unmapped1_fastq=args.unmapped1_fastq,
                    unmapped2_fastq=args.unmapped2_fastq,
                    insertions_csv=args.insertions_csv,
                    failed_csv=args.failed_csv,
                    cascade_csv=args.cascade_csv,
                    nuc_csv=args.nuc_csv,
                    amino_csv=args.amino_csv,
                    conseq_csv=args.conseq_csv,
                    conseq_all_csv=args.conseq_all_csv,
                    conseq_region_csv=args.conseq_region_csv,
                    failed_align_csv=args.failed_align_csv,
                    coverage_scores_csv=args.coverage_scores_csv,
                    aligned_csv=args.aligned_csv,
                    g2p_aligned_csv=args.g2p_aligned_csv,
                    unstitched_conseq_csv=args.unstitched_conseq_csv,
                    unstitched_contigs_csv=args.unstitched_contigs_csv,
                    contigs_csv=args.contigs_csv,
                    stitcher_plot_svg=args.stitcher_plot_svg,
                    genome_coverage_csv=args.genome_coverage_csv,
                    genome_coverage_svg=args.genome_coverage_svg,
                    genome_concordance_svg=args.genome_concordance_svg,
                    read_entropy_csv=args.read_entropy_csv,
                    concordance_csv=args.concordance_csv,
                    concordance_seed_csv=args.concordance_seed_csv,
                    scratch_path=scratch_path)
    sample.name = None  # Since the file names are messy in Kive.
    return sample


def main():
    logging.basicConfig(level=logging.WARN)
    args = parse_args()
    sample = load_sample(args)

    pssm = Pssm()
    sample.process(pssm,
                   force_gzip=True,  # dataset files change .gz to .raw
                   use_denovo=args.denovo)

    with tarfile.open(args.coverage_maps_tar, mode='w') as tar:
        for image_name in os.listdir(sample.coverage_maps):
            image_path = os.path.join(sample.coverage_maps, image_name)
            archive_path = os.path.join('coverage_maps', image_name)
            tar.add(image_path, archive_path)


if __name__ == '__main__':
    main()
