from argparse import ArgumentParser, FileType, ArgumentDefaultsHelpFormatter
from csv import DictReader
from io import StringIO
from itertools import groupby
from operator import itemgetter
from tempfile import NamedTemporaryFile
from pathlib import Path

from micall.utils.fasta_to_csv import fasta_to_csv


def parse_args():
    parser = ArgumentParser(
        description='Run a set of contigs through BLAST again.',
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('contigs_csv',
                        type=FileType(),
                        nargs='?',
                        default='contigs.csv',
                        help='contigs to search for')
    return parser.parse_args()


def main():
    args = parse_args()
    fasta_file = NamedTemporaryFile(mode='w', prefix='contigs', suffix='.fasta')
    contig_sources = []  # [(sample_name, contig_num, ref_name, contig_size)]
    ref_name = None
    for sample_name, sample_rows in groupby(DictReader(args.contigs_csv),
                                            itemgetter('sample')):
        for contig_num, row in enumerate(sample_rows, 1):
            ref_name = row['ref']
            header = f'>{sample_name}_{contig_num}-{ref_name}\n'
            fasta_file.write(header)
            fasta_file.write(row['contig'])
            fasta_file.write('\n')
            contig_size = len(row['contig'])
            contig_sources.append((sample_name,
                                   contig_num,
                                   ref_name,
                                   contig_size))
        if __name__ == '__live_coding__' and ref_name != 'unknown':
            break

    fasta_file.flush()
    new_contigs_csv = StringIO()
    blast_csv = StringIO()
    fasta_to_csv(Path(fasta_file.name), new_contigs_csv, blast_csv=blast_csv)
    blast_csv.seek(0)
    for source_contig_num, contig_rows in groupby(DictReader(blast_csv),
                                                  itemgetter('contig_num')):
        contig_rows = sorted(contig_rows, key=lambda r: int(r['score']))
        sample_name, contig_num, ref_name, contig_size = contig_sources[
            int(source_contig_num)-1]
        best_blast_hits = [None] * contig_size
        for row in contig_rows:
            if row['ref_name'] == 'HIV1-CON-XX-Consensus-seed':
                # Doesn't tell us about HIV subtype, so skip it.
                continue
            start = int(row['start'])
            end = int(row['end'])
            best_blast_hits[start-1:end] = [row['ref_name']] * (end-start+1)
        best_subtypes = set()
        matches = []
        for blast_ref, ref_positions in groupby(best_blast_hits):
            match_size = len(list(ref_positions))
            if match_size > 100 and blast_ref is not None:
                subtype = '-'.join(blast_ref.split('-')[:2])
                best_subtypes.add(subtype)
            else:
                blast_ref = 'other'
            matches.append(f'{blast_ref} x {match_size}')
        if len(best_subtypes) == 1:
            summary, = best_subtypes
        elif not best_subtypes:
            summary = None
        else:
            summary = ', '.join(matches)
        print(f'{sample_name}, {contig_num}-{ref_name}: {summary}')


if __name__ == '__main__':
    main()
