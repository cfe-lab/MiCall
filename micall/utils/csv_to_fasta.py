#! /usr/bin/env python

import argparse
import csv
from typing import TextIO, Sequence, Iterator
from pathlib import Path
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, Seq


def csv_to_fasta(contigs_csv: TextIO, contigs_fasta: Path) -> None:
    reader = csv.DictReader(contigs_csv)
    if reader.fieldnames is None or \
       'contig' not in reader.fieldnames:
        raise ValueError("Input CSV does not contain contigs.")

    def records() -> Iterator[SeqRecord]:
        for i, row in enumerate(reader):
            seq = row['contig']
            name = str(i + 1)
            yield SeqRecord(Seq.Seq(seq),
                            description='',
                            id=name,
                            name=name)

    SeqIO.write(records(), contigs_fasta, "fasta")


def main(argv: Sequence[str]):
    parser = argparse.ArgumentParser(description="Convert contigs from CSV to FASTA."
                                     " This converter assumes the MiCall's conventions"
                                     " for CSV field names.")
    parser.add_argument('contigs_csv', type=argparse.FileType('r'),
                        help="Input CSV file to read contigs from.")
    parser.add_argument('contigs_fasta', type=Path,
                        help="Output FASTA file to write contigs to.")
    args = parser.parse_args(argv)
    csv_to_fasta(args.contigs_csv, args.contigs_fasta)


if __name__ == "__main__":
    import sys
    main(sys.argv[1:])
