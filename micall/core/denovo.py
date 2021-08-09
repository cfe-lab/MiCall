import argparse
import logging
import os
import typing
from collections import Counter
from csv import DictWriter, DictReader
from datetime import datetime
from glob import glob
from io import StringIO
from itertools import groupby
from operator import itemgetter
from shutil import rmtree
from subprocess import run, PIPE, CalledProcessError, STDOUT
from tempfile import mkdtemp

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from micall.core.project_config import ProjectConfig

from tempfile import mkstemp
from shutil import move, copymode
from os import fdopen, remove
from iva.assembly import Assembly
from pyfastaq.tasks import deinterleave

HAPLOFLOW = "haploflow"
DEFAULT_DATABASE = os.path.join(os.path.dirname(__file__),
                                '..',
                                'blast_db',
                                'refs.fasta')
logger = logging.getLogger(__name__)


def write_contig_refs(contigs_fasta_path,
                      contigs_csv,
                      merged_contigs_csv=None,
                      blast_csv=None):
    """ Run BLAST search to identify contig sequences.

    :param str contigs_fasta_path: path to file to read contig sequences from
        and append merged contigs to
    :param contigs_csv: open file to write assembled contigs to
    :param merged_contigs_csv: open file to read contigs that were merged from
        amplicon reads
    :param blast_csv: open file to write BLAST search results for each contig
    """
    writer = DictWriter(contigs_csv,
                        ['ref', 'match', 'group_ref', 'contig'],
                        lineterminator=os.linesep)
    writer.writeheader()
    with open(contigs_fasta_path, 'a') as contigs_fasta:
        if merged_contigs_csv is not None:
            contig_reader = DictReader(merged_contigs_csv)
            for i, row in enumerate(contig_reader, 1):
                contig_name = f'merged-contig-{i}'
                contigs_fasta.write(f">{contig_name}\n{row['contig']}\n")
    group_refs = {}
    genotypes = genotype(contigs_fasta_path,
                         blast_csv=blast_csv,
                         group_refs=group_refs)
    genotype_count = 0
    for i, record in enumerate(SeqIO.parse(contigs_fasta_path, "fasta")):
        (ref_name, match_fraction) = genotypes.get(record.name, ('unknown', 0))
        seq = record.seq
        if match_fraction < 0:
            seq = seq.reverse_complement()
            match_fraction *= -1
        writer.writerow(dict(ref=ref_name,
                             match=match_fraction,
                             group_ref=group_refs.get(ref_name),
                             contig=seq))
        genotype_count += 1
    return genotype_count


def genotype(fasta, db=DEFAULT_DATABASE, blast_csv=None, group_refs=None):
    """ Use Blastn to search for the genotype of a set of reference sequences.

    :param str fasta: file path of the FASTA file containing the query
        sequences
    :param str db: file path of the database to search for matches
    :param blast_csv: open file to write the blast matches to, or None
    :param dict group_refs: {contig_ref: group_ref} or None. The dictionary
        will get filled in with the mapping from each contig's reference name
        to the best matched reference for the whole seed group.
    :return: {query_name: (ref_name, matched_fraction)} where query_name is a
        sequence header from the query sequences FASTA file, ref_name is the
        name of the best match from the database, and matched_fraction is the
        fraction of the query that aligned against the reference (matches and
        mismatches).
    """
    contig_nums = {}  # {contig_name: contig_num}
    with open(fasta) as f:
        for line in f:
            if line.startswith('>'):
                contig_name = line[1:-1]
                contig_nums[contig_name] = len(contig_nums) + 1
    blast_columns = ['qaccver',
                     'saccver',
                     'pident',
                     'score',
                     'qcovhsp',
                     'qstart',
                     'qend',
                     'sstart',
                     'send']
    cline = NcbiblastnCommandline(query=fasta,
                                  db=db,
                                  outfmt=f'"10 {" ".join(blast_columns)}"',
                                  evalue=0.0001,
                                  gapopen=5,
                                  gapextend=2,
                                  penalty=-3,
                                  reward=1,
                                  max_target_seqs=5000)
    stdout, _ = cline()
    samples = {}  # {query_name: (subject_name, matched_fraction)}
    matches = sorted(DictReader(StringIO(stdout), blast_columns),
                     key=lambda row: (row['qaccver'], float(row['score'])))
    if not blast_csv:
        blast_writer = None
    else:
        blast_writer = DictWriter(blast_csv,
                                  ['contig_num',
                                   'ref_name',
                                   'score',
                                   'match',
                                   'pident',
                                   'start',
                                   'end',
                                   'ref_start',
                                   'ref_end'],
                                  lineterminator=os.linesep)
        blast_writer.writeheader()
    contig_top_matches = {match['qaccver']: match['saccver']
                          for match in matches}
    top_refs = set(contig_top_matches.values())
    projects = ProjectConfig.loadDefault()
    match_scores = Counter()
    for contig_name, contig_matches in groupby(matches, itemgetter('qaccver')):
        contig_top_ref = contig_top_matches[contig_name]
        contig_seed_group = projects.getSeedGroup(contig_top_ref)
        for match in contig_matches:
            ref_name = match['saccver']
            if ref_name not in top_refs:
                continue
            match_seed_group = projects.getSeedGroup(ref_name)
            if match_seed_group == contig_seed_group:
                match_scores[ref_name] += float(match['score'])

    if group_refs is not None:
        group_top_refs = {projects.getSeedGroup(ref_name): ref_name
                          for ref_name, count in reversed(match_scores.most_common())}
        for ref_name in contig_top_matches.values():
            group_refs[ref_name] = group_top_refs[projects.getSeedGroup(ref_name)]

    for match in matches:
        matched_fraction = float(match['qcovhsp']) / 100
        if int(match['send']) < int(match['sstart']):
            matched_fraction *= -1
        pident = round(float(match['pident']))
        contig_name = match['qaccver']
        samples[contig_name] = (match['saccver'], matched_fraction)
        if blast_writer:
            blast_writer.writerow(dict(contig_num=contig_nums[contig_name],
                                       ref_name=match['saccver'],
                                       score=match['score'],
                                       match=matched_fraction,
                                       pident=pident,
                                       start=match['qstart'],
                                       end=match['qend'],
                                       ref_start=match['sstart'],
                                       ref_end=match['send']))
    return samples


def denovo(fastq1_path: str,
           fastq2_path: str,
           contigs_csv: typing.TextIO,
           work_dir: str = '.',
           merged_contigs_csv: typing.TextIO = None,
           blast_csv: typing.TextIO = None,
           haplo_args=None):
    """ Use de novo assembly to build contigs from reads.

    :param fastq1_path: FASTQ file name for read 1 reads
    :param fastq2_path: FASTQ file name for read 2 reads
    :param contigs_csv: open file to write assembled contigs to
    :param work_dir: path for writing temporary files
    :param merged_contigs_csv: open file to read contigs that were merged from
        amplicon reads
    :param blast_csv: open file to write BLAST search results for each contig
    """
    old_tmp_dirs = glob(os.path.join(work_dir, 'assembly_*'))
    for old_tmp_dir in old_tmp_dirs:
        rmtree(old_tmp_dir, ignore_errors=True)

    tmp_dir = mkdtemp(dir=work_dir, prefix='assembly_')
    start_time = datetime.now()
    start_dir = os.getcwd()
    joined_path = os.path.join(tmp_dir, 'joined.fastq')
    run(['merge-mates',
         fastq1_path,
         fastq2_path,
         '--interleave',
         '-o', joined_path],
        check=True)

    if haplo_args is None:
        haplo_args = {'long': 0,
                  'filter': 500,
                  'thres': -1,
                  'strict': 5,
                  'error': 0.02,
                  'kmer': 41}
    haplo_out_path = os.path.join(tmp_dir, 'haplo_out')
    contigs_fasta_path = os.path.join(haplo_out_path, 'contigs.fa')
    haplo_args = [HAPLOFLOW,
                  '--read-file', joined_path,
                  '--out', haplo_out_path,
                  '--k', str(haplo_args['kmer']),
                  '--error-rate', str(haplo_args['error']),
                  '--strict', str(haplo_args['strict']),
                  '--filter', str(haplo_args['filter']),
                  '--thres', str(haplo_args['thres']),
                  '--long', str(haplo_args['long']),
                  '--debug', '1']
    try:
        run(haplo_args, check=True, stdout=PIPE, stderr=STDOUT)
    except CalledProcessError as ex:
        output = ex.output and ex.output.decode('UTF8')
        if output != 'Failed to make first seed. Cannot continue\n':
            logger.warning('Haploflow failed to assemble.', exc_info=True)
            logger.warning(output)
        with open(contigs_fasta_path, 'a'):
            pass

    fh, abs_path = mkstemp()
    i = 0
    with fdopen(fh, 'w') as new_file:
        with open(contigs_fasta_path) as old_file:
            for line in old_file:
                if line.startswith('>'):
                    new_file.write(f">contig{i}\n")
                    i += 1
                else:
                    new_file.write(line)
    copymode(contigs_fasta_path, abs_path)
    remove(contigs_fasta_path)
    move(abs_path, contigs_fasta_path)
    print(f"Number of contigs before trimming and joining: {i}")

    haplo_assembly = Assembly(contigs_file=contigs_fasta_path)
    reads_prefix = os.path.join(tmp_dir, 'reads')
    reads_1 = reads_prefix + '_1.fa'
    reads_2 = reads_prefix + '_2.fa'
    deinterleave(joined_path, reads_1, reads_2, fasta_out=True)
    haplo_assembly._trim_strand_biased_ends(reads_prefix, tag_as_trimmed=True)
    haplo_assembly._remove_contained_contigs(list(haplo_assembly.contigs.keys()))
    haplo_assembly._merge_overlapping_contigs(list(haplo_assembly.contigs.keys()))
    contigs_fasta_path = os.path.join(haplo_out_path, 'finalcontigs.fasta')
    haplo_assembly.write_contigs_to_file(contigs_fasta_path)

    os.chdir(start_dir)
    duration = datetime.now() - start_time
    contig_count = write_contig_refs(contigs_fasta_path,
                                     contigs_csv,
                                     blast_csv=blast_csv)
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
