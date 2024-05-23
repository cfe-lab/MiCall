import typing
from typing import Optional, TextIO, Iterable, Dict, cast
from collections import Counter
from csv import DictWriter, DictReader
from io import StringIO
from itertools import groupby
from operator import itemgetter
import os

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline

from micall.core.project_config import ProjectConfig
from micall.utils.contig_stitcher_contigs import GenotypedContig


DEFAULT_DATABASE = os.path.join(os.path.dirname(__file__),
                                '..',
                                'blast_db',
                                'refs.fasta')


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


def run(contigs_fasta_path: str,
        unstitched_contigs_csv: TextIO,
        blast_csv: Optional[TextIO] = None):
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

    unstitched_writer = init_contigs_refs(cast(TextIO, unstitched_contigs_csv))
    group_refs: Dict[str, str] = {}

    genotypes = genotype(contigs_fasta_path,
                         blast_csv=blast_csv,
                         group_refs=group_refs)

    write_unstitched_contigs(unstitched_writer,
                             group_refs,
                             genotypes,
                             contigs_fasta_path)
    unstitched_contigs_csv.flush()
