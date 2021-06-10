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
from micall.core.remap import remap, map_to_contigs

IVA = "iva"
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


def run_iva(tmp_dir: str,
            joined_path: str,
            iva_out_path: str,
            contigs_fasta_path: str,
            merged_contigs_csv: typing.TextIO = None,
            is_pessimistic: bool = False):

    iva_args = [IVA, '--fr', joined_path, '-t', '2']
    if is_pessimistic:
        iva_args.append('--pessimistic')
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
        run(iva_args, check=True, stdout=PIPE, stderr=STDOUT)
    except CalledProcessError as ex:
        output = ex.output and ex.output.decode('UTF8')
        if output != 'Failed to make first seed. Cannot continue\n':
            logger.warning('iva failed to assemble.', exc_info=True)
            logger.warning(output)
        with open(contigs_fasta_path, 'a'):
            pass


def denovo(fastq1_path: str,
           fastq2_path: str,
           contigs_csv: typing.TextIO,
           work_dir: str = '.',
           merged_contigs_csv: typing.TextIO = None,
           blast_csv: typing.TextIO = None):
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

    if True: # set to pessimistic option!
        iva_out_path = os.path.join(tmp_dir, 'iva_out_pess')
        contigs_fasta_path = os.path.join(iva_out_path, 'contigs.fasta')
        run_iva(tmp_dir, joined_path, iva_out_path, contigs_fasta_path, merged_contigs_csv, is_pessimistic=True)
        contig_count = write_contig_refs(contigs_fasta_path, contigs_csv, blast_csv=blast_csv)
        # figure out which contigs are worth keeping from blast results!
        new_contigs_csv = contigs_csv #for now
        with open(os.path.join(tmp_dir, 'remap.csv'), 'w') as remap_csv, \
                open(os.path.join(tmp_dir, 'remap_counts.csv'), 'w') as counts_csv, \
                open(os.path.join(tmp_dir, 'remap_conseq_csv'), 'w') as conseq_csv, \
                open(os.path.join(tmp_dir, 'unmapped1.fastq'), 'w') as unmapped1, \
                open(os.path.join(tmp_dir, 'unmapped2_fastq'), 'w') as unmapped2:
            map_to_contigs(fastq1_path,
                           fastq2_path,
                           new_contigs_csv,
                           remap_csv,
                           counts_csv,
                           conseq_csv,
                           unmapped1,
                           unmapped2,
                           tmp_dir,)

    iva_out_path = os.path.join(tmp_dir, 'iva_out')
    contigs_fasta_path = os.path.join(iva_out_path, 'contigs.fasta')
    run_iva(tmp_dir, joined_path, iva_out_path, contigs_fasta_path, merged_contigs_csv, is_pessimistic=False)

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
