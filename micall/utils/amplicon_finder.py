from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import defaultdict, Counter
from csv import DictReader
from pathlib import Path


AMPLICON_REGIONS = ('HLA-B-exon2', 'HLA-B-exon3', 'V3LOOP')


def parse_args():
    parser = ArgumentParser(description='Summarize amplicon sizes.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    # noinspection PyTypeChecker
    parser.add_argument('raw_data',
                        default='~/data/RAW_DATA',
                        nargs='?',
                        type=Path)
    return parser.parse_args()


def find_good_coverage(version_path):
    old_coverage_scores = version_path / 'coverage_scores.csv'
    with old_coverage_scores.open() as f:
        reader = DictReader(f)
        found_regions = set()  # {(sample, region)}
        for row in reader:
            if row['on.score'] != '4':
                continue
            if row['region'] not in AMPLICON_REGIONS:
                if 'HLA' in row['region']:
                    print(row['region'], '!!!')
                continue
            found_regions.add((row['sample'], row['seed']))
    return found_regions


def count_contig_lengths(version_path: Path, size_counts: dict, good_amplicons: set):
    contigs_path = version_path / 'denovo' / 'contigs.csv'
    with contigs_path.open() as f:
        reader = DictReader(f)
        for row in reader:
            sample = row['sample']
            ref = row['ref']
            if ref.startswith('HIV1-'):
                old_ref = 'HIV1-CON-XX-Consensus-seed'
            else:
                old_ref = ref
            if (sample, old_ref) not in good_amplicons:
                continue
            seq = row['contig']
            contig_length = len(seq)
            if contig_length == 321:
                print(sample, ref)
            ref_size_counts = size_counts[old_ref]
            ref_size_counts[contig_length] += 1


def count_conseq_lengths(version_path: Path, size_counts: dict, good_amplicons: set):
    conseqs_path = version_path / 'conseq.csv'
    with conseqs_path.open() as f:
        reader = DictReader(f)
        for row in reader:
            if row['consensus-percent-cutoff'] != 'MAX':
                continue
            sample = row['sample']
            region = row['region']
            if (sample, region) not in good_amplicons:
                continue
            seed_size_counts = size_counts[region]
            amplicon_chunks = row['sequence'].split('x')
            for chunk in filter(None, amplicon_chunks):
                amplicon_length = len(chunk)
                seed_size_counts[amplicon_length] += 1


def main():
    args = parse_args()
    raw_data: Path = args.raw_data
    runs_path = raw_data.expanduser() / 'MiSeq' / 'runs'
    if not runs_path.exists():
        raise FileNotFoundError(runs_path)
    conseq_size_counts = defaultdict(Counter)  # {seed: {length: count}}
    contig_size_counts = defaultdict(Counter)  # {seed: {length: count}}
    for run_path in sorted(runs_path.iterdir()):
        if run_path.name == 'suspended':
            continue
        results_path = run_path / 'Results'
        version_paths = sorted(results_path.iterdir())
        if not version_paths:
            print(run_path.name, 'no results')
            continue
        version_path = version_paths[-1]
        print(run_path.name, version_path.name)
        good_amplicons = find_good_coverage(version_path)
        if not good_amplicons:
            continue
        count_conseq_lengths(version_path, conseq_size_counts, good_amplicons)
        count_contig_lengths(version_path, contig_size_counts, good_amplicons)
    print(sorted(conseq_size_counts.items()))
    print(sorted(contig_size_counts.items()))


if __name__ == '__main__':
    main()
