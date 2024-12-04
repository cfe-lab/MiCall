import argparse
import logging
from typing import Optional, TextIO, cast, BinaryIO
from datetime import datetime
from glob import glob
from shutil import rmtree, copyfileobj
from subprocess import CalledProcessError
import subprocess
from tempfile import mkdtemp
from pathlib import Path


HAPLOFLOW = "haploflow"
logger = logging.getLogger(__name__)


def count_fasta_sequences(file_path: Path):
    with open(file_path, 'r') as file:
        return sum(1 for line in file if line.startswith('>'))


def denovo(fastq1_path: Path,
           fastq2_path: Path,
           fasta: TextIO,
           work_dir: Path = Path('.'),
           merged_contigs_csv: Optional[TextIO] = None,
           ):
    """ Use de novo assembly to build contigs from reads.

    :param fastq1: FASTQ file for read 1 reads
    :param fastq2: FASTQ file for read 2 reads
    :param fasta: file to write assembled contigs to
    :param work_dir: path for writing temporary files
    :param merged_contigs_csv: open file to read contigs that were merged from
        amplicon reads
    """

    if merged_contigs_csv is not None:
        # TODO: implement this.
        logger.error("Haploflow implementation does not support contig extensions yet.")

    old_tmp_dirs = glob(str(work_dir / 'assembly_*'))
    for old_tmp_dir in old_tmp_dirs:
        rmtree(old_tmp_dir, ignore_errors=True)

    tmp_dir = Path(mkdtemp(dir=work_dir, prefix='assembly_'))

    start_time = datetime.now()
    joined_path = tmp_dir / 'joined.fastq'
    subprocess.run(['merge-mates',
                    str(fastq1_path),
                    str(fastq2_path),
                    '--interleave',
                    '-o', str(joined_path)],
                   check=True)

    haplo_args = {'long': 0,
                  'filter': 80,
                  'thres': -1,
                  'strict': 5,
                  'error': 0.02,
                  'kmer': 41,
                  'merge': False,
                  'scaffold': False,
                  'patch': False,
                  'ref': None,
                  'RP': False,
                  }

    assembly_out_path = tmp_dir / 'haplo_out'
    contigs_fasta_path = assembly_out_path / 'contigs.fa'

    assembly_out_path.mkdir(exist_ok=True, parents=True)
    contigs_fasta_path.touch()

    haplo_cmd = [HAPLOFLOW,
                 '--read-file', str(joined_path),
                 '--out', str(assembly_out_path),
                 '--k', str(haplo_args['kmer']),
                 '--error-rate', str(haplo_args['error']),
                 '--strict', str(haplo_args['strict']),
                 '--filter', str(haplo_args['filter']),
                 '--thres', str(haplo_args['thres']),
                 '--long', str(haplo_args['long'])]
    try:
        subprocess.run(haplo_cmd, check=True)
    except CalledProcessError:
        logger.warning('Haploflow failed to assemble.', exc_info=True)

    with open(contigs_fasta_path) as reader:
        copyfileobj(cast(BinaryIO, reader), fasta)
        fasta.flush()

    duration = datetime.now() - start_time
    contig_count = count_fasta_sequences(contigs_fasta_path)
    logger.info('Assembled %d contigs in %s (%ds) on %s.',
                contig_count,
                duration,
                duration.total_seconds(),
                fastq1_path)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser(
        description="A script to perform de novo assembly of reads to build contigs."
    )
    parser.add_argument(
        'fastq1',
        type=argparse.FileType('r'),
        help="Path to the FASTQ file containing read 1 of paired-end sequencing data."
    )
    parser.add_argument(
        'fastq2',
        type=argparse.FileType('r'),
        help="Path to the FASTQ file containing read 2 of paired-end sequencing data."
    )
    parser.add_argument(
        'fasta',
        type=argparse.FileType('w'),
        help="Path to the output FASTA file where assembled contigs will be written."
    )

    args = parser.parse_args()
    denovo(args.fastq1.name, args.fastq2.name, args.fasta)
