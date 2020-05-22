import argparse
from genetracks import Figure, Track, Multitrack, Coverage
import micall.core.plot_contigs as plot_contigs
from pathlib import Path
from csv import DictReader


def main(args):
    with open(args.blast_csv) as f:
        plot_blast(f)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'blast_csv',
        type=Path
    )
    args = parser.parse_args()
    return args


def plot_blast(blast_csv):
    reader = DictReader(blast_csv)
    visited = set()
    figure = Figure()
    for row in reader:
        print(row.contig_num)


if __name__ == '__main__':
    args = parse_args()
    main(args)