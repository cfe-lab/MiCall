import argparse
import logging
import os
import tempfile
from typing import Optional, TextIO, cast
from csv import DictReader
from datetime import datetime
from glob import glob
from shutil import rmtree, copyfileobj
from subprocess import PIPE, CalledProcessError, STDOUT
import subprocess
from tempfile import mkdtemp, NamedTemporaryFile

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import micall.utils.fasta_to_csv as fasta_to_csv
import micall.core.contig_stitcher as stitcher

IVA = "iva"
logger = logging.getLogger(__name__)

def write_contig_refs(contigs_fasta_path: str,
                      unstitched_contigs_csv: Optional[TextIO],
                      contigs_csv: Optional[TextIO],
                      merged_contigs_csv: Optional[TextIO] = None,
                      blast_csv: Optional[TextIO] = None,
                      stitcher_plot_path: Optional[str] = None) -> int:
    """ Run BLAST search to identify contig sequences.

    :param str contigs_fasta_path: path to file to read contig sequences from
        and append merged contigs to
    :param unstitched_contigs_csv: open file to write assembled contigs to
    :param contigs_csv: open file to write stitched contigs to
    :param merged_contigs_csv: open file to read contigs that were merged from
        amplicon reads
    :param blast_csv: open file to write BLAST search results for each contig
    :param stitcher_plot_path: open file to write the visualizer plot to
    """

    with open(contigs_fasta_path, 'a') as contigs_fasta:
        if merged_contigs_csv is not None:
            contig_reader = DictReader(merged_contigs_csv)
            for i, row in enumerate(contig_reader, 1):
                contig_name = f'merged-contig-{i}'
                contigs_fasta.write(f">{contig_name}\n{row['contig']}\n")

    with NamedTemporaryFile(mode='wt') as temporary_unstitched_csv:
        fasta_to_csv.run(contigs_fasta_path,
                         cast(TextIO, temporary_unstitched_csv),
                         blast_csv)

        if unstitched_contigs_csv:
            with open(temporary_unstitched_csv.name) as input_csv:
                copyfileobj(input_csv, unstitched_contigs_csv)

        if contigs_csv:
            output_csv = contigs_csv
        else:
            output_csv = open("/dev/null", "wt")

        with open(temporary_unstitched_csv.name) as input_csv:
            return stitcher.run(input_csv, output_csv, stitcher_plot_path)


def denovo(fastq1_path: str,
           fastq2_path: str,
           unstitched_contigs_csv: Optional[TextIO],
           contigs_csv: Optional[TextIO],
           work_dir: str = '.',
           merged_contigs_csv: Optional[TextIO] = None,
           blast_csv: Optional[TextIO] = None,
           stitcher_plot_path: Optional[str] = None,
           ):
    """ Use de novo assembly to build contigs from reads.

    :param fastq1_path: FASTQ file name for read 1 reads
    :param fastq2_path: FASTQ file name for read 2 reads
    :param unstitched_contigs_csv: open file to write assembled contigs to
    :param contigs_csv: open file to write stitched contigs to
    :param work_dir: path for writing temporary files
    :param merged_contigs_csv: open file to read contigs that were merged from
        amplicon reads
    :param blast_csv: open file to write BLAST search results for each contig
    :param stitcher_plot_path: open file to write the visualizer plot to
    """

    if unstitched_contigs_csv is None and contigs_csv is None:
        raise ValueError("Must specify either contigs_csv or unstitched_contigs_csv")

    old_tmp_dirs = glob(os.path.join(work_dir, 'assembly_*'))
    for old_tmp_dir in old_tmp_dirs:
        rmtree(old_tmp_dir, ignore_errors=True)

    tmp_dir = mkdtemp(dir=work_dir, prefix='assembly_')

    if contigs_csv is None:
        contigs_csv_tmp = tempfile.NamedTemporaryFile("wt")
        contigs_csv = cast(TextIO, contigs_csv_tmp.file)
    else:
        contigs_csv_tmp = None

    start_time = datetime.now()
    start_dir = os.getcwd()
    joined_path = os.path.join(tmp_dir, 'joined.fastq')
    if stitcher_plot_path is None:
        stitcher_plot_path = os.path.join(tmp_dir, "stitcher_plot.svg")
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

    os.chdir(start_dir)
    duration = datetime.now() - start_time
    contig_count = write_contig_refs(contigs_fasta_path,
                                     unstitched_contigs_csv,
                                     contigs_csv,
                                     blast_csv=blast_csv,
                                     stitcher_plot_path=stitcher_plot_path)
    logger.info('Assembled %d contigs in %s (%ds) on %s.',
                contig_count,
                duration,
                duration.total_seconds(),
                fastq1_path)

    if contigs_csv_tmp:
        contigs_csv_tmp.close()


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq1')
    parser.add_argument('fastq2')
    parser.add_argument('--unstitched_contigs', type=argparse.FileType('w'))
    parser.add_argument('--contigs', type=argparse.FileType('w'))
    parser.add_argument('--stitcher_plot')

    args = parser.parse_args()
    denovo(args.fastq1, args.fastq2, args.unstitched_contigs, args.contigs, args.stitcher_plot_path)
