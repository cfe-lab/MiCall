from argparse import ArgumentParser
from collections import defaultdict
from csv import DictReader
from itertools import groupby
from matplotlib import pyplot as plt
from operator import itemgetter
import os

import gotoh

from micall.core.project_config import ProjectConfig


def parse_args():
    parser = ArgumentParser(
        description='Search for the join point of chimeric sequences.')
    parser.add_argument('results', help='path to results folder')
    parser.add_argument('sample', help='sample name', nargs='+')
    parser.add_argument('--plot_path',
                        '-p',
                        help='path to write plots',
                        default='.')
    return parser.parse_args()


def plot(counts, filename, window_size=31):
    window_counts = defaultdict(dict)  # {(seed, consensus): {mid_pos: (agree, disagree)}}
    for names, agreement_counts in counts.iteritems():
        for mid_pos in agreement_counts.iterkeys():
            hit_count = 0
            window_start = mid_pos - window_size/2
            window_end = window_start + window_size
            agree_sum = disagree_sum = 0
            for pos in range(window_start, window_end):
                agree_count, disagree_count = agreement_counts.get(pos,
                                                                   (None, None))
                if agree_count is not None:
                    agree_sum += agree_count
                    disagree_sum += disagree_count
                    hit_count += 1
            if hit_count == window_size:
                window_counts[names][mid_pos] = (agree_sum/hit_count,
                                                 disagree_sum/hit_count)

    _fig, axes = plt.subplots(nrows=2, ncols=2)
    for axis_index, names in enumerate(sorted(window_counts.iterkeys())):
        axis = axes[axis_index % 2][axis_index / 2]
        seed_name, consensus_name = names
#         axis.set_yscale('log')
        axis.set_ylim((0.5, 20000))
        axis.set_xlim((2500, 5500))
        axis.set_xlabel(consensus_name + ' nuc pos')
        axis.set_ylabel(seed_name + ' match counts')
        name_counts = window_counts[names]
        x = sorted(name_counts.iterkeys())
        y1 = [name_counts[pos][0] for pos in x]
        y2 = [sum(name_counts[pos]) for pos in x]
        axis.step(x, y1, 'r', linewidth=2, where='mid', label='agree')
        axis.step(x, y2, 'k', linewidth=2, where='mid', label='total')
        axis.legend()
        plt.tight_layout()

    plt.savefig(filename)  # write image to file
    for i in range(2):
        for j in range(2):
            axes[i][j].set_xlim((3000, 3500))
    plt.savefig(filename + '_zoom.png')


def main():
    args = parse_args()
    projects = ProjectConfig.loadDefault()
    for sample_name in args.sample:
        process_file(sample_name, projects, args)
    print('Done.')


def find_consensus(rows, seed_name, region_name):
    max_nucs = {}  # {pos: max_nuc}
    matching_rows = (row
                     for row in rows
                     if row['seed'] == seed_name and
                     row['region'] == region_name)
    for row in matching_rows:
        max_count = 0
        pos = int(row['refseq.nuc.pos'])
        for nuc in 'ACGT':
            nuc_count = int(row[nuc])
            if nuc_count > max_count:
                max_nucs[pos] = nuc
                max_count = nuc_count
    max_pos = max_nucs and max(max_nucs.iterkeys()) or 0
    consensus = ''.join(max_nucs.get(i+1, 'N') for i in xrange(max_pos))
    return consensus


def process_file(sample_name, projects, args):
    print('Starting {}.'.format(sample_name))

    nuc_counts = defaultdict(dict)  # {(source, dest): {pos: (agree, disagree)}}
    nucleotide_path = os.path.join(args.results, 'nuc.csv')
    with open(nucleotide_path, 'rU') as nuc_csv:
        reader = DictReader(nuc_csv)
        sample_rows = (row for row in reader
                       if row['sample'] == sample_name)
        region_rows = (row for row in sample_rows
                       if row['region'].endswith('-NS2') or
                       row['region'].endswith('-NS3'))
        combo_rows = {combo: list(group)
                      for combo, group in groupby(region_rows,
                                                  itemgetter('seed', 'region'))}
    seed_names = {combo[0] for combo in combo_rows.iterkeys()}
    seed_map = {seed_name: projects.getReference(seed_name)
                for seed_name in seed_names}
    consensus_map = {combo: find_consensus(rows, *combo)
                     for combo, rows in combo_rows.iteritems()}
    nucleotides = 'ACGT'
    for combo, rows in combo_rows.iteritems():
        seed_name, _region_name = combo
        consensus = consensus_map[combo]
        for dest_seed_name in seed_names:
            dest_seed = seed_map[dest_seed_name]
            # print seed_name, region_name, dest_seed_name
            seed_positions = map_sequence(dest_seed, consensus)
            for row in rows:
                pos = int(row['refseq.nuc.pos'])
                agree_count = disagree_count = 0
                if pos > len(seed_positions):
                    continue
                seed_pos = seed_positions[pos-1]
                if seed_pos is None:
                    continue
                expected_nuc = dest_seed[seed_pos-1]
                for nuc in nucleotides:
                    count = int(row[nuc])
                    if nuc == expected_nuc:
                        agree_count += count
                    else:
                        disagree_count += count
                nuc_counts[(dest_seed_name, seed_name)][seed_pos] = (agree_count,
                                                                     disagree_count)

    plot_path_pattern = os.path.join(args.plot_path,
                                     '{prefix}{sample}_agreement.png')
    # amino_plot_path = plot_path_pattern.format(sample=sample_name, prefix='aa_')
    nuc_plot_path = plot_path_pattern.format(sample=sample_name, prefix='')
    # plot(amino_counts, amino_plot_path)
    plot(nuc_counts, nuc_plot_path)


def map_sequence(source_seq, dest_seq):
    """ Find the portion of source_seq that dest_seq maps to.

    :return: a list of 1-based positions in source_seq that it mapped to.
    """
    gap_open = 15
    gap_extend = 5
    use_terminal_gap_penalty = 1
    aligned_source, aligned_dest, _score = gotoh.align_it(source_seq,
                                                          dest_seq,
                                                          gap_open,
                                                          gap_extend,
                                                          use_terminal_gap_penalty)
    positions = []
    source_pos = 1
    for source_nuc, dest_nuc in zip(aligned_source, aligned_dest):
        if dest_nuc != '-':
            positions.append(source_pos if source_nuc != '-' else None)
        if source_nuc != '-':
            source_pos += 1

    return positions


if __name__ == '__main__':
    main()
