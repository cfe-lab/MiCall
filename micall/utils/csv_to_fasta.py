#! /usr/bin/env python

import sys
import argparse
import csv
from typing import TextIO, Sequence, Iterator
from pathlib import Path
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, Seq


class NoContigsInCSV(ValueError):
    pass


COLUMNS = ['contig', 'sequence']


def csv_to_fasta(contigs_csv: TextIO, contigs_fasta: Path) -> None:
    reader = csv.DictReader(contigs_csv)
    seqcolumns = [col for col in COLUMNS if col in (reader.fieldnames or [])]
    if len(seqcolumns) == 0:
        raise NoContigsInCSV("Input CSV does not contain contigs.")
    if len(seqcolumns) > 1:
        raise NoContigsInCSV("Input CSV contains multiple possible contig columns.")

    seqcolumn = seqcolumns[0]

    def records() -> Iterator[SeqRecord]:
        for i, row in enumerate(reader):
            seq = row[seqcolumn]
            name = str(i + 1)
            yield SeqRecord(Seq.Seq(seq),
                            description='',
                            id=name,
                            name=name)

    contigs_fasta.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records(), contigs_fasta, "fasta")


def main(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(description="Convert contigs from CSV to FASTA."
                                     " This converter assumes the MiCall's conventions"
                                     " for CSV field names.")
    parser.add_argument('contigs_csv', type=argparse.FileType('r'),
                        help="Input CSV file to read contigs from.")
    parser.add_argument('contigs_fasta', type=Path,
                        help="Output FASTA file to write contigs to.")
    args = parser.parse_args(argv)
    csv_to_fasta(args.contigs_csv, args.contigs_fasta)
    return 0


def cli() -> None:
    sys.exit(main(sys.argv[1:]))


if __name__ == "__main__": cli()  # noqa
