#! /usr/bin/env python3

"""
Kive-style bowtie2
Run bowtie2 on paired-end FASTQ data sets with user-supplied *.bt2
bowtie2 SAM format output to <stdout> for redirection via subprocess.Popen
Sort outputs by refname.
Convert to CSV format and write to file.
"""

import argparse
import csv
import logging
import os
import sys
from pathlib import Path
from typing import Optional, Sequence, Set

from micall.core import project_config
from micall.utils.cache import cached
from micall.utils.externals import Bowtie2, Bowtie2Build, LineCounter
from micall.utils.stderr import Stderr
from micall.utils.work_dir import WorkDir

BOWTIE_THREADS = 1    # Bowtie performance roughly scales with number of threads
# Read and reference gap open/extension penalties.
READ_GAP_OPEN = 10
READ_GAP_EXTEND = 3
REF_GAP_OPEN = 10
REF_GAP_EXTEND = 3

logger = logging.getLogger(__name__)
line_counter = LineCounter()


@cached("prelim_map", parameters=['nthreads', 'rdgopen', 'rfgopen', 'gzip', 'excluded_seeds'], outputs=['prelim_csv'])
def prelim_map(
    fastq1: Path,
    fastq2: Path,
    prelim_csv: Path,
    nthreads: int = BOWTIE_THREADS,
    rdgopen: int = READ_GAP_OPEN,
    rfgopen: int = REF_GAP_OPEN,
    gzip: bool = False,
    excluded_seeds: Optional[Set[str]] = None,
    ) -> None:

    bowtie2 = Bowtie2()
    bowtie2_build = Bowtie2Build()
    bowtie2_build.set_logger(logger)

    # Get work_dir and stderr from dynamic scope - required to be set by caller
    work_path = WorkDir.get()
    stderr_file = Stderr.get()

    # check that the inputs exist
    fastq1_str = check_fastq(str(fastq1), gzip)
    fastq2_str = check_fastq(str(fastq2), gzip)

    # generate initial reference files
    projects = project_config.ProjectConfig.loadDefault()
    ref_path = work_path / 'micall.fasta'
    all_excluded_seeds = {project_config.G2P_SEED_NAME}
    if excluded_seeds:
        all_excluded_seeds.update(excluded_seeds)
    with open(ref_path, 'w') as ref:
        projects.writeSeedFasta(ref, all_excluded_seeds)
    reffile_template = str(work_path / 'reference')
    bowtie2_build.build(str(ref_path), reffile_template)

    fieldnames = ['qname',
                  'flag',
                  'rname',
                  'pos',
                  'mapq',
                  'cigar',
                  'rnext',
                  'pnext',
                  'tlen',
                  'seq',
                  'qual']

    # Get stderr from dynamic scope (already a file object)
    # No need to open it - it's already open via Stderr.using()
    with open(prelim_csv, 'w') as csv_file:
        writer = csv.writer(csv_file, lineterminator=os.linesep)
        writer.writerow(fieldnames)

        # do preliminary mapping
        read_gap_open_penalty = rdgopen
        ref_gap_open_penalty = rfgopen

        # stream output from bowtie2
        bowtie_args = ['--wrapper', 'micall-0',
                       '--quiet',
                       '-x', reffile_template,
                       '-1', fastq1_str,
                       '-2', fastq2_str,
                       '--rdg', "{},{}".format(read_gap_open_penalty,
                                               READ_GAP_EXTEND),
                       '--rfg', "{},{}".format(ref_gap_open_penalty,
                                               REF_GAP_EXTEND),
                       '--no-hd',  # no header lines (start with @)
                       '-X', '1200',
                       '-p', str(nthreads)]

        for i, line in enumerate(bowtie2.yield_output(bowtie_args, stderr=stderr_file)):
            writer.writerow(line.split('\t')[:11])  # discard optional items


def check_fastq(filename: str, gzip: bool = False) -> str:
    if not os.path.exists(filename):
        sys.exit('No FASTQ found at ' + filename)
    if gzip:
        if not filename.endswith('.gz'):
            new_filename = filename + '.gz'
            try:
                os.symlink(filename, new_filename)
            except FileExistsError:
                pass
            filename = new_filename
    return filename


def main(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        description='Map contents of FASTQ R1 and R2 data sets to references using bowtie2.')

    parser.add_argument('fastq1', help='<input> FASTQ containing forward reads')
    parser.add_argument('fastq2', help='<input> FASTQ containing reverse reads')
    parser.add_argument('prelim_csv',
                        help='<output> CSV containing preliminary mapping from bowtie2 (modified SAM)')
    parser.add_argument("--rdgopen", type=int, default=READ_GAP_OPEN, help="<optional> read gap open penalty")
    parser.add_argument("--rfgopen", type=int, default=REF_GAP_OPEN, help="<optional> reference gap open penalty")
    parser.add_argument("--gzip", action='store_true', help="<optional> FASTQs are compressed")

    args = parser.parse_args(argv)
    prelim_csv_path = Path(args.prelim_csv)
    with WorkDir.using(prelim_csv_path.parent):
        prelim_map(fastq1=Path(args.fastq1),
                   fastq2=Path(args.fastq2),
                   prelim_csv=prelim_csv_path,
                   rdgopen=args.rdgopen,
                   rfgopen=args.rfgopen,
                   gzip=args.gzip)

    return 0


def entry() -> None:
    sys.exit(main(sys.argv[1:]))


if __name__ == '__main__': entry()  # noqa
