import csv
import os
from argparse import ArgumentParser, FileType
from collections import Counter, defaultdict
from itertools import groupby
from operator import itemgetter

from matplotlib import pyplot as plt

from micall.utils.samples_from_454 import format_sample_name


def parse_args():
    parser = ArgumentParser(description='Compare results with 454 data.')
    parser.add_argument('source',
                        type=FileType('r'),
                        help='454 data in CSV')
    parser.add_argument('result_path', help='Micall results folder')
    return parser.parse_args()


def format_percent(numerator, denominator):
    return '{:.0f}'.format(100.0 * numerator / denominator) if denominator else ''


def plot(prelim_percents):
    levels = sorted(prelim_percents.keys())
    levels.append(levels.pop(0))  # move None to end
    labels = ['<= {} seqs'.format(limit) for limit in levels]
    labels[-1] = 'more'
    level_counts = [prelim_percents[limit] for limit in levels]
    plt.hist(level_counts, stacked=True, label=labels)
    plt.xlabel('% mapped (prelim)')
    plt.ylabel('# samples')
    plt.ylim(0, 10000)
    plt.title('Preliminary mapping success (Consensus B ref, soft clip)')
    plt.legend()
    plt.savefig('compare_454_samples.png')


def main():
    args = parse_args()
    with args.source:
        source_reader = csv.DictReader(args.source)
        source_counts = {}
        for (run, sample), rows in groupby(source_reader,
                                           itemgetter('run', 'enum')):
            sample_name = format_sample_name(run, sample) + '_R2'
            source_counts[sample_name] = counts = Counter()
            for row in rows:
                counts['raw'] += 1
                fpr = float(row['g2p_score'])  # mislabelled in file
                if fpr > 3.5:
                    counts['R5'] += 1
                else:
                    counts['X4'] += 1

    remap_counts_path = os.path.join(args.result_path, 'remap_counts.csv')
    with open(remap_counts_path, 'r') as remap_counts_file:
        raw_counts = Counter()
        prelim_counts = Counter()
        remap_counts = Counter()
        remap_counts_reader = csv.DictReader(remap_counts_file)
        for row in remap_counts_reader:
            sample = row['sample']
            count_type = row['type']
            row_count = int(row['count']) // 6
            if count_type == 'raw':
                raw_counts[sample] = row_count
            elif count_type == 'unmapped':
                pass
            else:
                category, ref_name = count_type.split(' ')
                if ref_name == '*':
                    pass
                elif category == 'prelim':
                    prelim_counts[sample] += row_count
                elif category == 'remap-final':
                    remap_counts[sample] += row_count

    coverage_scores = Counter()
    coverage_path = os.path.join(args.result_path, 'coverage_scores.csv')
    with open(coverage_path, 'r') as coverage_file:
        for row in csv.DictReader(coverage_file):
            if row['region'] == 'V3LOOP':
                coverage_scores[row['sample']] = int(row['on.score'])

    g2p_path = os.path.join(args.result_path, 'g2p.csv')
    with open(g2p_path, 'r') as g2p_results:
        result_reader = csv.DictReader(g2p_results)
        result_counts = defaultdict(Counter)
        for row in result_reader:
            row_count = int(row['count']) // 3
            counts = result_counts[row['sample']]
            category = row['call'] or row['error']
            counts[category] += row_count

    with open('compare_454_samples.csv', 'w') as result_file:
        writer = csv.DictWriter(result_file,
                                ['sample',
                                 'raw',
                                 'prelim',
                                 'prelim_pct',
                                 'remap',
                                 'remap_pct',
                                 'cv_score',
                                 'R5_old',
                                 'R5_new',
                                 'R5_diff_pct',
                                 'X4_old',
                                 'X4_new',
                                 'X4_diff_pct',
                                 'cysteines',
                                 'notdiv3',
                                 'length',
                                 'ambig',
                                 'unmapped',
                                 'stops',
                                 'low_qual',
                                 'err_unmap_pct'])
        writer.writeheader()
        read_limits = (10, 50, 100, None)
        prelim_percents = {limit: [] for limit in read_limits}
        for sample in sorted(raw_counts.keys()):
            results = result_counts[sample]
            source_sample_counts = source_counts.pop(sample)
            raw_count = source_sample_counts['raw']
            if raw_count != raw_counts[sample]:
                raise RuntimeError('Expected {} raw, but found {} for {}.'.format(
                    raw_count,
                    raw_counts[sample],
                    sample))
            unmapped_count = raw_count - sum(results.values())
            r5_old = source_sample_counts['R5']
            r5_new = results.pop('R5', 0)
            r5_diff = format_percent(r5_new - r5_old, r5_old)
            x4_old = source_sample_counts['X4']
            x4_new = results.pop('X4', 0)
            x4_diff = format_percent(x4_new - x4_old, x4_old)
            cysteines = results.pop('cysteines', 0)
            notdiv3 = results.pop('notdiv3', 0)
            length = results.pop('length', 0)
            ambiguous = results.pop('> 2 ambiguous', 0)
            stops = results.pop('stop codons', 0)
            low_quality = results.pop('low quality', 0)
            assert not results, results
            err_unmapped_pct = format_percent(cysteines +
                                              notdiv3 +
                                              length +
                                              ambiguous +
                                              unmapped_count +
                                              low_quality,
                                              raw_count)
            prelim_count = prelim_counts[sample]
            prelim_pct = format_percent(prelim_count, raw_count)
            for limit in read_limits:
                if limit is None or raw_count <= limit:
                    prelim_percents[limit].append(float(prelim_pct))
                    break
            writer.writerow(dict(sample=sample,
                                 raw=raw_count,
                                 prelim=prelim_count,
                                 prelim_pct=prelim_pct,
                                 remap=remap_counts[sample],
                                 remap_pct=format_percent(remap_counts[sample],
                                                          raw_count),
                                 cv_score=coverage_scores[sample],
                                 R5_old=r5_old,
                                 R5_new=r5_new,
                                 R5_diff_pct=r5_diff,
                                 X4_old=x4_old,
                                 X4_new=x4_new,
                                 X4_diff_pct=x4_diff,
                                 cysteines=cysteines,
                                 notdiv3=notdiv3,
                                 length=length,
                                 ambig=ambiguous,
                                 stops=stops,
                                 unmapped=unmapped_count,
                                 low_qual=low_quality,
                                 err_unmap_pct=err_unmapped_pct))

    if source_counts:
        print('{} missing sources:'.format(len(source_counts)))
        print(', '.join(sorted(source_counts.keys())))

    plot(prelim_percents)
    print('Done.')

if __name__ == '__main__':
    main()
