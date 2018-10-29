import argparse
import logging
import os
import sys
from csv import DictWriter, DictReader
from datetime import datetime
from glob import glob
from io import StringIO
from shutil import rmtree
from subprocess import Popen
from tempfile import mkdtemp

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline

from micall.utils.iva_wrapper import assemble, SeedingError

PEAR = "/opt/bin/pear"
SAVAGE = "/opt/savage_wrapper.sh"
IVA = "iva"
IS_SAVAGE_ENABLED = False
DEFAULT_DATABASE = os.path.join(os.path.dirname(__file__),
                                '..',
                                'blast_db',
                                'refs.fasta')
ASSEMBLY_TIMEOUT = 1800  # 30 minutes
logger = logging.getLogger(__name__)


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


def denovo(fastq1_path, fastq2_path, contigs, work_dir='.', merged_contigs_csv=None):
    prefix = "sample"

    old_tmp_dirs = glob(os.path.join(work_dir, 'iva_*'))
    for old_tmp_dir in old_tmp_dirs:
        rmtree(old_tmp_dir, ignore_errors=True)

    tmp_dir = mkdtemp(dir=work_dir, prefix='iva_')
    start_time = datetime.now()
    start_dir = os.getcwd()
    try:
        assemble(os.path.join(tmp_dir, "iva"),
                 str(fastq1_path),
                 str(fastq2_path))
        is_successful = True
    except SeedingError:
        logger.warning('De novo assembly failed.')
        is_successful = False
    os.chdir(start_dir)
    duration = datetime.now() - start_time

    if IS_SAVAGE_ENABLED:
        pear_proc = Popen([PEAR, "-f", fastq1_path, "-r", fastq2_path, "-o", prefix],
                          cwd=work_dir)

        if pear_proc.wait():
            raise Exception

        savage_proc = Popen([
                    SAVAGE, "--split", "1", "-s",
                    "{}/sample.assembled.fastq".format(work_dir), "-p1",
                    "{}/sample.unassembled.forward.fastq".format(work_dir),
                    "-p2",
                    "{}/sample.unassembled.reverse.fastq".format(work_dir),
                    "-t", "1",
                    "--merge_contigs", "0.01", "--overlap_len_stage_c",
                    "100",
            ], cwd=work_dir, shell=True)

        savage_proc.wait(timeout=ASSEMBLY_TIMEOUT)
        if not savage_proc.wait():
            write_genotypes("{}/contigs_stage_c.fasta".format(work_dir), contigs)
        else:
            print("savage exits with error", file=sys.stderr)

    if is_successful:
        contigs_fasta_path = os.path.join(tmp_dir, 'iva', 'contigs.fasta')
    else:
        contigs_fasta_path = os.path.join(tmp_dir, 'contigs.fasta')
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
