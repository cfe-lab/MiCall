import argparse
import logging
import os
import tempfile
import typing
from typing import Optional, TextIO, Iterable, Dict, cast
from collections import Counter
from csv import DictWriter, DictReader
from datetime import datetime
from glob import glob
from io import StringIO
from itertools import groupby
from operator import itemgetter
from shutil import rmtree, copyfileobj
from subprocess import run, PIPE, CalledProcessError, STDOUT
from tempfile import mkdtemp, NamedTemporaryFile

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from micall.core.project_config import ProjectConfig
from micall.utils.contig_stitcher_contigs import GenotypedContig
import micall.core.contig_stitcher as stitcher

IVA = "iva"
DEFAULT_DATABASE = os.path.join(os.path.dirname(__file__),
                                '..',
                                'blast_db',
                                'refs.fasta')
logger = logging.getLogger(__name__)


def read_assembled_contigs(group_refs, genotypes, contigs_fasta_path: str) -> Iterable[GenotypedContig]:
    projects = ProjectConfig.loadDefault()

    for i, record in enumerate(SeqIO.parse(contigs_fasta_path, "fasta")):
        (ref_name, match_fraction) = genotypes.get(record.name, ('unknown', 0))
        seq = record.seq
        if match_fraction < 0:
            seq = seq.reverse_complement()
            match_fraction *= -1

        group_ref = group_refs.get(ref_name)
        try:
            ref_seq = projects.getGenotypeReference(group_ref)
        except KeyError:
            try:
                ref_seq = projects.getReference(group_ref)
            except KeyError:
                ref_seq = None

        yield GenotypedContig(name=record.name,
                              seq=str(seq),
                              ref_name=ref_name,
                              group_ref=group_ref,
                              ref_seq=str(ref_seq) if ref_seq is not None else None,
                              match_fraction=match_fraction)


def init_contigs_refs(contigs_csv: TextIO):
    writer = DictWriter(contigs_csv,
                        ['ref', 'match', 'group_ref', 'contig'],
                        lineterminator=os.linesep)
    writer.writeheader()
    return writer


def write_unstitched_contigs(writer,
                             group_refs,
                             genotypes,
                             contigs_fasta_path
                             ):

    for contig in read_assembled_contigs(group_refs, genotypes, contigs_fasta_path):
        writer.writerow(dict(ref=contig.ref_name,
                             match=contig.match_fraction,
                             group_ref=contig.group_ref,
                             contig=contig.seq
                             ))


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
        unstitched_writer = init_contigs_refs(cast(TextIO, temporary_unstitched_csv))
        stitched_writer = init_contigs_refs(contigs_csv) if contigs_csv else None
        group_refs: Dict[str, str] = {}

        genotypes = genotype(contigs_fasta_path,
                             blast_csv=blast_csv,
                             group_refs=group_refs)

        write_unstitched_contigs(unstitched_writer,
                                 group_refs,
                                 genotypes,
                                 contigs_fasta_path)
        temporary_unstitched_csv.flush()

        if unstitched_contigs_csv:
            with open(temporary_unstitched_csv.name) as input_csv:
                copyfileobj(input_csv, unstitched_contigs_csv)

        with open(temporary_unstitched_csv.name) as input_csv:
            return stitcher.parse_and_run(input_csv, stitched_writer, stitcher_plot_path)


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

    contig_nums: Dict[str, int] = {}  # {contig_name: contig_num}
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
    match_scores: typing.Counter[str] = Counter()
    for contig_name, contig_matches in groupby(matches, itemgetter('qaccver')):
        contig_top_ref = contig_top_matches[contig_name]
        contig_seed_group = projects.getSeedGroup(contig_top_ref)
        for match in contig_matches:
            ref_name = match['saccver']
            if ref_name not in top_refs:
                continue
            match_seed_group = projects.getSeedGroup(ref_name)
            if match_seed_group == contig_seed_group:
                match_scores[ref_name] += float(match['score'])  # type: ignore[assignment]

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
    run(['merge-mates',
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
        run(iva_args, check=True, stdout=PIPE, stderr=STDOUT)
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
