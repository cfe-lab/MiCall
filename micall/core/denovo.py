import argparse
import logging
import os
from csv import DictWriter, DictReader
from datetime import datetime, timedelta
from enum import Enum
from glob import glob
from io import StringIO
from shutil import rmtree
from subprocess import run, PIPE, CalledProcessError
from tempfile import mkdtemp

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline

from micall.utils.externals import LineCounter

PEAR = "/opt/bin/pear"
SAVAGE = "/opt/savage_wrapper.sh"
IVA = "iva"
DEFAULT_DATABASE = os.path.join(os.path.dirname(__file__),
                                '..',
                                'blast_db',
                                'refs.fasta')
ASSEMBLY_TIMEOUT = timedelta(hours=4).total_seconds()
GENOME_WIDTH = 5000  # Wild guess at the genome width.
TARGET_DEPTH = 1000  # Savage recommends 500 to 1000.
logger = logging.getLogger(__name__)

# noinspection PyArgumentList
Assembler = Enum('Assembler', 'IVA SAVAGE')


def write_genotypes(contigs_fasta_path, contigs_csv, merged_contigs_csv=None):
    writer = DictWriter(contigs_csv,
                        ['genotype', 'match', 'contig'],
                        lineterminator=os.linesep)
    writer.writeheader()
    with open(contigs_fasta_path, 'a') as contigs_fasta:
        if merged_contigs_csv is not None:
            contig_reader = DictReader(merged_contigs_csv)
            for i, row in enumerate(contig_reader, 1):
                contig_name = f'merged-contig-{i}'
                contigs_fasta.write(f">{contig_name}\n{row['contig']}\n")
    genotypes = genotype(contigs_fasta_path)
    genotype_count = 0
    for i, record in enumerate(SeqIO.parse(contigs_fasta_path, "fasta")):
        genotype_name, match_fraction = genotypes.get(record.name, ('unknown', 0))
        writer.writerow(dict(genotype=genotype_name,
                             match=match_fraction,
                             contig=record.seq))
        genotype_count += 1
    return genotype_count


def genotype(fasta, db=DEFAULT_DATABASE):
    """ Use Blastn to search for the genotype of a set of reference sequences.

    :param str fasta: file path of the FASTA file containing the query
        sequences
    :param str db: file path of the database to search for matches
    :return: {query_name: (ref_name, matched_fraction)} where query_name is a
        sequence header from the query sequences FASTA file, ref_name is the
        name of the best match from the database, and matched_fraction is the
        fraction of the query that aligned against the reference (matches and
        mismatches).
    """
    cline = NcbiblastnCommandline(query=fasta,
                                  db=db,
                                  outfmt='"10 qaccver saccver pident score qcovhsp"',
                                  evalue=0.0001,
                                  gapopen=5,
                                  gapextend=2,
                                  penalty=-3,
                                  reward=1,
                                  max_target_seqs=5000)
    stdout, _ = cline()
    samples = {}  # {query_name: (subject_name, matched_fraction)}
    matches = sorted(DictReader(StringIO(stdout), ['qaccver',
                                                   'saccver',
                                                   'pident',
                                                   'score',
                                                   'qcovhsp']),
                     key=lambda row: float(row['score']))
    for match in matches:
        matched_fraction = float(match['qcovhsp']) / 100
        samples[match['qaccver']] = (match['saccver'], matched_fraction)
    return samples


def denovo(fastq1_path,
           fastq2_path,
           contigs,
           work_dir='.',
           merged_contigs_csv=None,
           assembler=Assembler.IVA):
    old_tmp_dirs = glob(os.path.join(work_dir, 'assembly_*'))
    for old_tmp_dir in old_tmp_dirs:
        rmtree(old_tmp_dir, ignore_errors=True)

    tmp_dir = mkdtemp(dir=work_dir, prefix='assembly_')
    start_time = datetime.now()
    start_dir = os.getcwd()
    contigs_fasta_path = os.path.join(tmp_dir, 'contigs.fasta')
    if assembler == Assembler.IVA:
        joined_path = os.path.join(tmp_dir, 'joined.fastq')
        run(['merge-mates',
             fastq1_path,
             fastq2_path,
             '--interleave',
             '-o', joined_path],
            check=True)
        iva_out_path = os.path.join(tmp_dir, 'iva_out')
        contigs_fasta_path = os.path.join(iva_out_path, 'contigs.fasta')
        try:
            run([IVA, '--fr', joined_path, '-t', '8', iva_out_path], check=True)
        except CalledProcessError:
            logger.warning('iva failed to assemble.', exc_info=True)
            with open(contigs_fasta_path, 'a'):
                pass
    else:
        counter = LineCounter()
        read_count = counter.count(fastq1_path) / 4
        expected_depth = read_count / GENOME_WIDTH
        split = max(1, expected_depth // TARGET_DEPTH)
        merged_reads_path = os.path.join(tmp_dir, 'merged.fastq')
        run(['merge-mates',
             '--out', merged_reads_path,
             fastq1_path,
             fastq2_path],
            check=True)

        try:
            run(['savage',
                 '--split', str(split),
                 '-s', merged_reads_path,
                 '-t', '1',
                 '--merge_contigs', '0.01',
                 '--overlap_len_stage_c', '100'],
                cwd=tmp_dir,
                check=True,
                stdout=PIPE)
            contigs_fasta_path = os.path.join(tmp_dir, 'contigs_stage_c.fasta')
        except CalledProcessError:
            logger.warning('De novo assembly failed.')

    os.chdir(start_dir)
    duration = datetime.now() - start_time
    contig_count = write_genotypes(contigs_fasta_path, contigs, merged_contigs_csv)
    logger.info('Assembled %d contigs in %s (%ds) on %s.',
                contig_count,
                duration,
                duration.total_seconds(),
                fastq1_path)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq1')
    parser.add_argument('fastq2')
    parser.add_argument('contigs', type=argparse.FileType('w'))

    args = parser.parse_args()
    denovo(args.fastq1, args.fastq2, args.contigs)
