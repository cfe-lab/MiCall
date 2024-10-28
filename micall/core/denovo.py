import argparse
import logging
import os
from typing import Optional, TextIO
from csv import DictReader
from datetime import datetime
from glob import glob
from shutil import rmtree, copyfileobj
from subprocess import PIPE, CalledProcessError, STDOUT
import subprocess
from tempfile import mkdtemp

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


IVA = "iva"
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
    iva_out_path = os.path.join(tmp_dir, 'iva_out')
    contigs_fasta_path = os.path.join(iva_out_path, 'contigs.fasta')
    iva_args = [IVA, '--fr', joined_path, '-t', '2']
    if merged_contigs_csv is not None:
        seeds_fasta_path = os.path.join(tmp_dir, 'seeds.fasta')
        with open(seeds_fasta_path, 'w') as seeds_fasta:
            SeqIO.write((SeqRecord(Seq(row['contig']), f'seed-{i}', '', '')
                         for i, row in enumerate(DictReader(merged_contigs_csv))),
                        seeds_fasta,
                        'fasta')
            seeds_size = seeds_fasta.tell()
        if seeds_size > 0:
            iva_args.extend(['--contigs', seeds_fasta_path, '--make_new_seeds'])
    iva_args.append(iva_out_path)
    try:
        subprocess.run(iva_args, check=True, stdout=PIPE, stderr=STDOUT)
    except CalledProcessError as ex:
        output = ex.output and ex.output.decode('UTF8')
        if output != 'Failed to make first seed. Cannot continue\n':
            logger.warning('iva failed to assemble.', exc_info=True)
            logger.warning(output)
        with open(contigs_fasta_path, 'a'):
            pass

    with open(contigs_fasta_path, 'rb') as reader:
        copyfileobj(reader, fasta)  # type: ignore

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
