import argparse
import logging
import os
import typing
from typing import Optional, TextIO, Iterable, Dict, cast, Sequence
from collections import Counter
from csv import DictWriter, DictReader
from itertools import groupby
from operator import itemgetter
from pathlib import Path

from io import StringIO

from Bio import SeqIO

from micall.core.project_config import ProjectConfig
from micall.utils.contig_stitcher_contigs import GenotypedContig
from micall.utils.externals import Blastn, DefaultBlastDatabase


default_database = DefaultBlastDatabase().path


def read_assembled_contigs(group_refs: Dict[str, str],
                           genotypes: Dict[str, typing.Tuple[str, float]],
                           contigs_fasta: Path) -> Iterable[GenotypedContig]:
    """Read assembled contigs and generate GenotypedContig objects.

    Args:
        group_refs (Dict[str, str]): Mapping of reference names to group references.
        genotypes (Dict[str, Tuple[str, float]]): Mapping of contig names to (reference name, match fraction).
        contigs_fasta (Path): Path to the FASTA file containing contig sequences.

    Returns:
        Iterable[GenotypedContig]: An iterable of GenotypedContig objects.
    """
    projects = ProjectConfig.loadDefault()

    for i, record in enumerate(SeqIO.parse(contigs_fasta, "fasta")):
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
                              match_fraction=match_fraction,
                              reads_count=None)


def init_contigs_refs(contigs_csv: TextIO) -> DictWriter:
    """Initialize a CSV writer with header for contig references.

    Args:
        contigs_csv (TextIO): Open file object to write the contig references.

    Returns:
        DictWriter: A CSV DictWriter object initialized with the headers.
    """
    writer = DictWriter(contigs_csv,
                        ['ref', 'match', 'group_ref', 'contig'],
                        lineterminator=os.linesep)
    writer.writeheader()
    return writer


def write_contigs(writer: DictWriter,
                  group_refs: Dict[str, str],
                  genotypes: Dict[str, typing.Tuple[str, float]],
                  contigs_fasta: Path):
    """Write contigs to a CSV file.

    Args:
        writer (DictWriter): CSV writer to write contigs.
        group_refs (Dict[str, str]): Mapping of reference names to group references.
        genotypes (Dict[str, Tuple[str, float]]): Mapping of contig names to (reference name, match fraction).
        contigs_fasta (Path): Path to the FASTA file containing contig sequences.
    """
    for contig in read_assembled_contigs(group_refs, genotypes, contigs_fasta):
        writer.writerow(dict(ref=contig.ref_name,
                             match=contig.match_fraction,
                             group_ref=contig.group_ref,
                             contig=contig.seq))


def genotype(contigs_fasta: Path, db: Optional[Path] = None,
             blast_csv: Optional[TextIO] = None,
             group_refs: Optional[Dict[str, str]] = None) -> Dict[str, typing.Tuple[str, float]]:
    """Use Blastn to search for the genotype of a set of reference sequences.

    Args:
        fasta (Path): File path of the FASTA file containing the query sequences.
        db (Optional[str]): File path of the database to search for matches.
        blast_csv (Optional[TextIO]): Open file to write the blast matches to, or None.
        group_refs (Optional[Dict[str, str]]): Dictionary to fill with the mapping from
           each contig's reference name to the best matched reference for the whole seed group.

    Returns:
        Dict[str, Tuple[str, float]]: Mapping of query name to (reference name, matched fraction).
    """

    contig_nums: Dict[str, int] = {}  # {contig_name: contig_num}

    for record in SeqIO.parse(contigs_fasta, "fasta"):
        contig_name = record.name
        contig_nums[contig_name] = len(contig_nums) + 1

    def invoke_blast(db: Path) -> str:
        return Blastn().genotype(
            contigs_fasta=contigs_fasta,
            database=db,
        )

    if db is None:
        with default_database() as db:
            stdout = invoke_blast(db)
    else:
        stdout = invoke_blast(db)

    samples = {}  # {query_name: (subject_name, matched_fraction)}
    matches = sorted(DictReader(StringIO(stdout), Blastn.BLAST_COLUMNS),
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


def fasta_to_csv(contigs_fasta: Path,
                 contigs_csv: TextIO,
                 blast_csv: Optional[TextIO] = None) -> None:
    """Run BLAST search to identify contig sequences and write them to CSV.

    Args:
        contigs_fasta (Path): Path to the FASTA file containing contig sequences.
        contigs_csv (TextIO): Open file to write assembled contigs to.
        blast_csv (Optional[TextIO]): Open file to write BLAST search results for each contig.
    """

    contigs_fasta.touch()
    writer = init_contigs_refs(cast(TextIO, contigs_csv))
    group_refs: Dict[str, str] = {}

    genotypes = genotype(contigs_fasta, blast_csv=blast_csv, group_refs=group_refs)

    write_contigs(writer, group_refs, genotypes, contigs_fasta)
    contigs_csv.flush()


def main(argv: Sequence[str]):
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser(description="Convert contigs from FASTA to CSV format with BLAST annotations.")
    parser.add_argument('contigs_fasta', help="Input FASTA file with contig sequences.")
    parser.add_argument('contigs_csv', type=argparse.FileType('w'),
                        help="Output CSV file to write assembled contigs.")
    parser.add_argument('--blast_csv', type=argparse.FileType('w'),
                        help="Optional CSV file to write BLAST search results.")
    args = parser.parse_args(argv)
    fasta_to_csv(Path(args.contigs_fasta), args.contigs_csv, args.blast_csv)


if __name__ == "__main__":
    import sys
    main(sys.argv[1:])
