import os
from collections import Counter
from csv import DictWriter, DictReader
from pathlib import Path
from subprocess import run, PIPE

import numpy as np
import matplotlib
from scipy.stats import entropy

matplotlib.use('Agg')
from matplotlib import pyplot as plt  # noqa


def merge_reads(fastq1_path, fastq2_path, joined_fastq_path, merge_lengths_csv, use_gzip=False):
    args = ['merge-mates',
            fastq1_path,
            fastq2_path,
            '-o', joined_fastq_path,
            '--csv']
    if use_gzip:
        args.append('--gzip')
    result = run(args, check=True, stdout=PIPE, encoding='UTF8')
    count_summary = result.stdout.strip().split(',')
    bin_counts = count_summary[2:]
    if bin_counts[0].startswith('['):
        bin_counts[0] = bin_counts[0][1:]
    if bin_counts[-1].endswith(']'):
        bin_counts[-1] = bin_counts[-1][:-1]
    writer = DictWriter(merge_lengths_csv,
                        ['merge_length', 'count'],
                        lineterminator=os.linesep)
    writer.writeheader()
    for merge_length, count_text in enumerate(bin_counts):
        count = int(count_text)
        if count:
            writer.writerow(dict(merge_length=merge_length, count=count))


def plot_merge_lengths(read_entropy_path):
    plots = []
    fig, ax1 = plt.subplots()
    # noinspection PyTypeChecker
    merge_entropy = np.genfromtxt(read_entropy_path,
                                  delimiter=',',
                                  names=True,
                                  dtype=float)
    merge_lengths = merge_entropy['merge_length'].astype(int)
    merge_counts = merge_entropy['count'].astype(int)
    plots.append(plt.bar(merge_lengths,
                         merge_counts,
                         width=1.0,
                         label='read count'))
    plt.title('Count of merged reads at each length')
    ax1.set_ylabel('merged read count')
    ax1.set_xlabel('length of merged reads')

    ax2 = ax1.twinx()
    plots.append(ax2.plot(merge_entropy['merge_length'],
                          merge_entropy['entropy'],
                          'r',
                          label='entropy')[0])

    max_length = merge_lengths.max()
    min_buffer = 3
    max_buffer = 7
    length_counts = np.zeros(max_length + 1)
    length_counts[merge_lengths] = merge_counts
    buffered_length_counts = np.zeros(max_length + 1 + 2*max_buffer)
    buffered_length_counts[max_buffer:-max_buffer] = length_counts
    for buffer_size in range(min_buffer, max_buffer+1):
        averages = np.zeros(max_length+1)
        for neighbour_offset in range(-buffer_size, buffer_size+1):
            if neighbour_offset:
                start_index = max_buffer + neighbour_offset
                averages += buffered_length_counts[
                            start_index:start_index+max_length+1]
        averages /= buffer_size * 2
        # noinspection PyArgumentList
        ratios = np.divide(length_counts,
                           averages,
                           out=np.zeros_like(length_counts),
                           where=averages != 0)
        ratios[ratios < 1.0] = np.nan
        # plots.extend(ax2.plot(ratios,
        #                       'o',
        #                       alpha=0.3,
        #                       label=f'+/-{buffer_size} nbrs'))
    ax2.set_ylabel('entropy level')
    ax2.set_ylim(0)
    plt.legend(plots, [p.get_label() for p in plots])
    return fig


def write_merge_lengths_plot(read_entropy_path, plot_path):
    fig = plot_merge_lengths(read_entropy_path)
    fig.savefig(plot_path)
    plt.close(fig)


def count_kmers(sequence, counts, k=12):
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        counts[kmer] += 1


def calculate_entropy_from_counts(counts):
    return entropy(list(counts.values()), base=2.0)


def calculate_entropy(merged_fastq, merge_lengths):
    entropy_levels = {}
    for merge_length in merge_lengths:
        merged_fastq.seek(0)
        kmer_counts = Counter()
        for header, sequence, divider, quality in zip(*(merged_fastq,)*4):
            sequence = sequence.rstrip()  # trim line feed
            if len(sequence) != merge_length:
                continue
            count_kmers(sequence, kmer_counts)
        entropy_levels[merge_length] = calculate_entropy_from_counts(kmer_counts)
    return entropy_levels


def merge_for_entropy(fastq1_path,
                      fastq2_path,
                      read_entropy_csv,
                      work_path,
                      use_gzip=False):
    work_path = Path(work_path)
    writer = DictWriter(read_entropy_csv,
                        ['merge_length', 'count', 'entropy'],
                        lineterminator=os.linesep)
    writer.writeheader()

    merged_fastq_path = work_path / 'merged.fastq'
    merge_lengths_path = work_path / 'merge_lengths.csv'
    with merge_lengths_path.open('w') as merge_lengths_csv:
        merge_reads(fastq1_path,
                    fastq2_path,
                    merged_fastq_path,
                    merge_lengths_csv,
                    use_gzip)

    with merge_lengths_path.open() as merge_lengths_csv:
        merge_lengths = {int(row['merge_length']): int(row['count'])
                         for row in DictReader(merge_lengths_csv)}
    with merged_fastq_path.open() as merged_fastq:
        entropy_levels = calculate_entropy(merged_fastq, merge_lengths)

    for merge_length, entropy_level in sorted(entropy_levels.items()):
        read_count = merge_lengths[merge_length]
        writer.writerow(dict(merge_length=merge_length,
                             count=read_count,
                             entropy=entropy_level))


def main():
    print('Starting.')
    scratch_paths = dict(
        hla=Path('/home/don/git/MiCall/micall/tests/working') /
        'basespace_linked_data_hla_160729_M01841/scratch/71258A-HLA-B_S16',
        hla2=Path('/home/don/git/MiCall/micall/tests/working') /
        'basespace_linked_data/scratch/71258A-HLA-B_S16',
        hcv=Path('/home/don/git/MiCall/micall/tests/working/') /
        'basespace_linked_data_v3loop_hiv_hcv_190621_M04401/scratch' /
        '90317A-HCV_S100',
        hiv=Path('/home/don/git/MiCall/micall/tests/working/') /
        'basespace_linked_data_v3loop_hiv_hcv_190621_M04401/scratch' /
        '2569P1Y02945-1E8-HIV_S32',
        tcr=Path('/home/don/git/MiCall/micall/tests/working/') /
        'basespace_linked_data_tcr_190501_M05995/scratch/15PA2-HLA-B_S2',
        v3=Path('/home/don/git/MiCall/micall/tests/working/') /
        'basespace_linked_data/scratch' /
        '84681A-V3-WYD-T2-1-V3LOOP_S40')
    scratch_path = scratch_paths['v3']
    read_entropy_path = scratch_path / 'read_entropy.csv'
    plot_merge_lengths(read_entropy_path)
    plt.show()
    print('Done.')


if __name__ in ('__main__', '__live_coding__'):
    main()
