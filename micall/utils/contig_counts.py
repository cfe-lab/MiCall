import typing
from argparse import ArgumentParser, OPTIONAL, FileType, ArgumentDefaultsHelpFormatter
from collections import defaultdict, Counter
from csv import DictReader
from io import StringIO
from itertools import groupby
from operator import itemgetter
from pathlib import Path


def create_parser():
    parser = ArgumentParser(
        description='Report counts on all contigs at some reference positions',
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('start',
                        type=int,
                        help='Start position in reference (one-based)')
    parser.add_argument('end',
                        type=int,
                        help='End position in reference (one-based, included)')
    parser.add_argument('scratch',
                        nargs=OPTIONAL,
                        default='.',
                        type=FileSet(['genome_coverage.csv', 'aligned.csv']),
                        help='Scratch folder to search for contigs')
    parser.add_argument('-r', '--ref',
                        help='Prefix of references to include')
    return parser


class FileSet:
    def __init__(self, file_names: typing.List[str]):
        self.file_names = file_names

    def __call__(self, string):
        files = {}
        start_path = Path(string)
        file_type = FileType()
        for file_name in self.file_names:
            file_path = str(start_path / file_name)
            files[file_name] = file_type(file_path)
        return files


class ContigCounts:
    def __init__(self, start: int, end: int):
        self.start = start
        self.end = end
        self.reference_prefix = ''

        # {contig_name: {contig_pos: ref_pos}}
        self.position_map: typing.Dict[str, typing.Dict[int, int]]
        self.position_map = defaultdict(dict)

        self.counts: typing.Dict[str, typing.Dict[int, typing.Dict[str, int]]]
        self.counts = defaultdict(lambda: defaultdict(Counter))

    def read_genome_coverage(self, genome_coverage_csv: typing.IO):
        for row in DictReader(genome_coverage_csv):
            if not row['coordinates'].startswith(self.reference_prefix):
                continue
            ref_pos_text = row['refseq_nuc_pos']
            if not ref_pos_text:
                continue
            ref_pos = int(ref_pos_text)
            if self.start <= ref_pos <= self.end:
                contig_name = row['contig']
                contig_map = self.position_map[contig_name]
                contig_pos_text = row['query_nuc_pos']
                if not contig_pos_text:
                    continue
                contig_pos = int(contig_pos_text)
                contig_map[contig_pos] = ref_pos

    def read_aligned(self, aligned_csv: typing.IO):
        for contig_name, contig_rows in groupby(DictReader(aligned_csv),
                                                itemgetter('refname')):
            contig_positions = self.position_map[contig_name]
            if not contig_positions:
                continue
            contig_counts = self.counts[contig_name]
            start = min(contig_positions.keys())
            end = max(contig_positions.keys())
            for row in contig_rows:
                sequence = row['seq']
                offset = int(row['offset'])
                count = int(row['count'])
                for contig_pos in range(start, end+1):
                    nuc_index = contig_pos - 1 - offset
                    if nuc_index < 0:
                        continue
                    if len(sequence) <= nuc_index:
                        continue
                    nuc = sequence[nuc_index]
                    if nuc in 'nN':
                        # Middle of read pair or low quality, don't count it.
                        continue
                    ref_pos = contig_positions.get(contig_pos)
                    if ref_pos is not None:
                        contig_counts[ref_pos][nuc] += count

    def display(self) -> str:
        report = StringIO()
        contig_groups = sorted(self.counts.items())
        combined = defaultdict(Counter)
        for contig_name, contig_counts in contig_groups:
            for position, position_counts in list(contig_counts.items()):
                position_counts = Counter(position_counts)
                contig_counts[position] = position_counts
                combined[position] += position_counts
        if 1 < len(contig_groups):
            contig_groups.append(('Combined', combined))
        for contig_name, contig_counts in contig_groups:
            print(contig_name + ':', file=report)
            for position, position_counts in sorted(contig_counts.items()):
                position_counts = Counter(position_counts)
                report.write(f'{position}: ')
                coverage = sum(position_counts.values())
                position_display = []
                for nuc, nuc_count in position_counts.most_common():
                    prevalence = f'{nuc_count/coverage:0.2f}'
                    if prevalence != '0.00':
                        position_display.append(f'{nuc} ({prevalence})')
                print(*position_display, sep=', ', file=report)
        return report.getvalue()


def main(command_args: typing.List[str] = None):
    parser = create_parser()
    args = parser.parse_args(command_args)
    scratch_files = args.scratch
    genome_coverage_csv = scratch_files['genome_coverage.csv']
    aligned_csv = scratch_files['aligned.csv']
    counts = ContigCounts(args.start, args.end)
    if args.ref:
        counts.reference_prefix = args.ref
    counts.read_genome_coverage(genome_coverage_csv)
    counts.read_aligned(aligned_csv)
    print(counts.display(), end='')


if __name__ == '__main__':
    main()
