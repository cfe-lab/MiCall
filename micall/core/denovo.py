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


MEGAHIT = "megahit"
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
        logger.error("Megahit implementation does not support contig extensions yet.")

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
    output_path = os.path.join(tmp_dir, 'assembler_output')
    contigs_fasta_path = os.path.join(output_path, 'final.contigs.fa')
    megahit_cmd = [MEGAHIT, '--12', joined_path, '-o', output_path]
    try:
        subprocess.run(megahit_cmd, check=True, stdout=PIPE, stderr=STDOUT)
    except CalledProcessError as ex:
        output = ex.output and ex.output.decode('UTF8')
        if output != 'Failed to make first seed. Cannot continue\n':
            logger.warning('iva failed to assemble.', exc_info=True)
            logger.warning(output)
        with open(contigs_fasta_path, 'a'):
            pass

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
