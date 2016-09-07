from argparse import ArgumentParser
from collections import defaultdict
from csv import DictReader
from matplotlib import pyplot as plt
import os

import gotoh

from micall.core.project_config import ProjectConfig
# from micall.core.aln2counts import AMINO_ALPHABET
# from micall.utils.translation import translate


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


def plot(counts, filename, window_size=51):
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


def main():
    args = parse_args()
    projects = ProjectConfig.loadDefault()
    for sample_name in args.sample:
        process_file(sample_name, projects, args)
    print('Done.')


def find_region_rows(f, sample_name):
    reader = DictReader(f)
    sample_rows = (row for row in reader
                   if row['sample'] == sample_name)
    region_rows = (row for row in sample_rows
                   if row['region'].endswith('-NS2') or
                   row['region'].endswith('-NS3'))
    return region_rows


def process_file(sample_name, projects, args):
    print('Starting {}.'.format(sample_name))

    # amino_counts = defaultdict(dict)  # {(seed, consensus): {pos: (agree, disagree
    # amino_path = os.path.join(result_path, 'amino.csv')
    # with open(amino_path, 'rU') as amino_csv:
    #     region_rows = find_region_rows(amino_csv, sample_name)
    #     seed_map = {row['seed']: projects.getReference(row['seed'])
    #                 for row in region_rows}
    #     amino_csv.seek(0)
    #     region_rows = find_region_rows(amino_csv, sample_name)
    #     for row in region_rows:
    #         amino_pos = int(row['query.aa.pos'])
    #         codon_end = amino_pos*3 - 1
    #         codon_start = codon_end-3
    #         for seed_name in seed_map.keys():
    #             consensus_name = row['seed']
    #             expected_nucs = position_map[seed_name][consensus_name]
    #             expected_codon = expected_nucs[codon_start:codon_end]
    #             expected_amino = translate(expected_codon)
    #             agree_count = disagree_count = 0
    #             for amino in AMINO_ALPHABET:
    #                 count = int(row[amino])
    #                 if amino == expected_amino:
    #                     agree_count += count
    #                 else:
    #                     disagree_count += count
    #             amino_counts[(seed_name, consensus_name)][amino_pos] = (
    #                 agree_count,
    #                 disagree_count)

    nuc_counts = defaultdict(dict)  # {(seed, consensus): {pos: (agree, disagree)}}
    nucleotide_path = os.path.join(args.results, 'nuc.csv')
    with open(nucleotide_path, 'rU') as nuc_csv:
        region_rows = find_region_rows(nuc_csv, sample_name)
        seed_names = {row['seed'] for row in region_rows}
        seed_map = {seed_name: projects.getReference(seed_name)
                    for seed_name in seed_names}
        # seed_map['HCV-2k'] = projects.getReference('HCV-3a')
        position_map = build_position_map(seed_map)
        nuc_csv.seek(0)
        region_rows = find_region_rows(nuc_csv, sample_name)
        nucleotides = 'ACGT'
        report = []
        for row in region_rows:
            pos = int(row['query.nuc.pos'])
            for seed_name in seed_map.iterkeys():
                consensus_name = row['seed']
                expected_nucs = position_map[seed_name][consensus_name]
                if pos-2 >= len(expected_nucs):
                    expected_nuc = '-'
                else:
                    expected_nuc = expected_nucs[pos-2]
                agree_count = disagree_count = 0
                max_nuc = None
                max_count = 0
                for nuc in nucleotides:
                    count = int(row[nuc])
                    if count > max_count:
                        max_count = count
                        max_nuc = nuc
                    if nuc == expected_nuc:
                        agree_count += count
                    else:
                        disagree_count += count
                nuc_counts[(seed_name, consensus_name)][pos] = (agree_count,
                                                                disagree_count)
                if 2900 <= pos <= 3000 or 4000 <= pos <= 4100:
                    report.append(' '.join(map(str, (seed_name,
                                                     consensus_name,
                                                     pos,
                                                     expected_nuc,
                                                     max_nuc,
                                                     agree_count,
                                                     disagree_count))))
    if 'HCV-2k' in seed_map:
        report.sort()
        for line in report:
            print(line)
        print(seed_map['HCV-2k'])

    plot_path_pattern = os.path.join(args.plot_path,
                                     '{prefix}{sample}_agreement.png')
    # amino_plot_path = plot_path_pattern.format(sample=sample_name, prefix='aa_')
    nuc_plot_path = plot_path_pattern.format(sample=sample_name, prefix='')
    # plot(amino_counts, amino_plot_path)
    plot(nuc_counts, nuc_plot_path)


def map_sequence(source_seq, dest_seq):
    gap_open = 15
    gap_extend = 5
    use_terminal_gap_penalty = 1
    aligned_source, aligned_dest, _score = gotoh.align_it(source_seq,
                                                          dest_seq,
                                                          gap_open,
                                                          gap_extend,
                                                          use_terminal_gap_penalty)
    seed_nucs = []
    dest_index = 0
    for source_nuc, dest_nuc in zip(aligned_source, aligned_dest):
        if dest_nuc != '-' or dest_nuc == dest_seq[dest_index]:
            seed_nucs.append(source_nuc)
            dest_index += 1
        if dest_index >= len(dest_seq):
            break

    mapped_seq = ''.join(seed_nucs)
    return mapped_seq


def build_position_map(seed_map):
    position_map = defaultdict(dict)  # {seed_name: {consensus_name: seed_seq}}
    for source_region, source_seq in seed_map.iteritems():
        for dest_region, dest_seq in seed_map.iteritems():
            mapped_seq = map_sequence(source_seq, dest_seq)

            position_map[source_region][dest_region] = mapped_seq
            # if source_region == 'HCV-1a' and dest_region == 'HCV-1a':
            #     print source_seq
            #     print dest_seq
            #     print position_map[source_region][dest_region]
    return position_map

if __name__ == '__main__':
    main()
