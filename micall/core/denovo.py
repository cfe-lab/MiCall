import argparse
import logging
import os
from typing import Optional, TextIO, cast, BinaryIO
from datetime import datetime
from glob import glob
from shutil import rmtree, copyfileobj
from subprocess import PIPE, CalledProcessError, STDOUT
import subprocess
from tempfile import mkdtemp


HAPLOFLOW = "haploflow"
logger = logging.getLogger(__name__)


def count_fasta_sequences(file_path):
    with open(file_path, 'r') as file:
        return sum(1 for line in file if line.startswith('>'))


def denovo(fastq1_path: str,
           fastq2_path: str,
           fasta: TextIO,
           work_dir: str = '.',
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

    old_tmp_dirs = glob(os.path.join(work_dir, 'assembly_*'))
    for old_tmp_dir in old_tmp_dirs:
        rmtree(old_tmp_dir, ignore_errors=True)

    tmp_dir = mkdtemp(dir=work_dir, prefix='assembly_')

    start_time = datetime.now()
    start_dir = os.getcwd()
    joined_path = os.path.join(tmp_dir, 'joined.fastq')
    subprocess.run(['merge-mates',
                    fastq1_path,
                    fastq2_path,
                    '--interleave',
                    '-o', joined_path],
                   check=True)

    haplo_args = {'long': 0,
                  'filter': 500,
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
    assembly_out_path = os.path.join(tmp_dir, 'haplo_out')
    contigs_fasta_path = os.path.join(assembly_out_path, 'contigs.fa')

    os.makedirs(assembly_out_path, exist_ok=True)
    with open(contigs_fasta_path, 'w'):
        pass

    haplo_cmd = [HAPLOFLOW,
                 '--read-file', joined_path,
                 '--out', assembly_out_path,
                 '--k', str(haplo_args['kmer']),
                 '--error-rate', str(haplo_args['error']),
                 '--strict', str(haplo_args['strict']),
                 '--filter', str(haplo_args['filter']),
                 '--thres', str(haplo_args['thres']),
                 '--long', str(haplo_args['long'])]
    try:
        subprocess.run(haplo_cmd, check=True, stdout=PIPE, stderr=STDOUT)
    except CalledProcessError as ex:
        output = ex.output and ex.output.decode('UTF8')
        if output != 'Failed to make first seed. Cannot continue\n':
            logger.warning('Haploflow failed to assemble.', exc_info=True)
            logger.warning(output)

    with open(contigs_fasta_path) as reader:
        copyfileobj(cast(BinaryIO, reader), fasta)

    os.chdir(start_dir)
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
