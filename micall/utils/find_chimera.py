from collections import defaultdict
from csv import DictReader
from matplotlib import pyplot as plt
import os

import gotoh

from micall.core.project_config import ProjectConfig
from micall.core.aln2counts import AMINO_ALPHABET
from micall.utils.translation import translate

result_path = (
    '/media/raw_data/MiSeq/runs/160715_M01841_0246_000000000-ARF87/'
    'Results/version_7.4')
sample_names = ['73087A-HCV_S2',
                '73103A-HCV_S4',
                '73117A-HCV_S82',
                '73111A-HCV_S5',
                '73090A-HCV_S41']
plot_path_pattern = '/home/don/data/{prefix}{sample}_agreement.png'


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
    projects = ProjectConfig.loadDefault()
    for sample_name in sample_names:
        process_file(sample_name, projects)
    print('Done.')


def process_file(sample_name, projects):
    print('Starting {}.'.format(sample_name))
    consensus_path = os.path.join(result_path, 'conseq.csv')
    with open(consensus_path, 'rU') as conseq_csv:
        reader = DictReader(conseq_csv)
        sample_rows = (row for row in reader
                       if row['sample'] == sample_name and
                       row['consensus-percent-cutoff'] == 'MAX')
        consensus_map = {}
        seed_map = {}
        for row in sample_rows:
            region = row['region']
            offset = int(row['offset'])
            consensus = row['sequence']
            consensus_map[region] = '-'*offset + consensus
            seed = projects.getReference(region)
            seed_map[region] = seed
    position_map = build_position_map(consensus_map, seed_map)

    amino_counts = defaultdict(dict)  # {(seed, consensus): {pos: (agree, disagree
    amino_path = os.path.join(result_path, 'amino.csv')
    with open(amino_path, 'rU') as amino_csv:
        reader = DictReader(amino_csv)
        sample_rows = (row for row in reader
                       if row['sample'] == sample_name)
        region_rows = (row for row in sample_rows
                       if row['region'].endswith('-NS2') or
                       row['region'].endswith('-NS3'))
        for row in region_rows:
            amino_pos = int(row['query.aa.pos'])
            codon_end = amino_pos*3 - 1
            codon_start = codon_end-3
            for seed_name in seed_map.keys():
                consensus_name = row['seed']
                expected_nucs = position_map[seed_name][consensus_name]
                expected_codon = expected_nucs[codon_start:codon_end]
                expected_amino = translate(expected_codon)
                agree_count = disagree_count = 0
                for amino in AMINO_ALPHABET:
                    count = int(row[amino])
                    if amino == expected_amino:
                        agree_count += count
                    else:
                        disagree_count += count
                amino_counts[(seed_name, consensus_name)][amino_pos] = (
                    agree_count,
                    disagree_count)

    nuc_counts = defaultdict(dict)  # {(seed, consensus): {pos: (agree, disagree)}}
    nucleotide_path = os.path.join(result_path, 'nuc.csv')
    with open(nucleotide_path, 'rU') as nuc_csv:
        reader = DictReader(nuc_csv)
        sample_rows = (row for row in reader
                       if row['sample'] == sample_name)
        region_rows = (row for row in sample_rows
                       if row['region'].endswith('-NS2') or
                       row['region'].endswith('-NS3'))
        nucleotides = 'ACGT'
        for row in region_rows:
            pos = int(row['query.nuc.pos'])
            for seed_name in seed_map.keys():
                consensus_name = row['seed']
                expected_nucs = position_map[seed_name][consensus_name]
                if pos-2 >= len(expected_nucs):
                    expected_nuc = '-'
                else:
                    expected_nuc = expected_nucs[pos-2]
                agree_count = disagree_count = 0
                for nuc in nucleotides:
                    count = int(row[nuc])
                    if nuc == expected_nuc:
                        agree_count += count
                    else:
                        disagree_count += count
                nuc_counts[(seed_name, consensus_name)][pos] = (agree_count,
                                                                disagree_count)
    amino_plot_path = plot_path_pattern.format(sample=sample_name, prefix='aa_')
    nuc_plot_path = plot_path_pattern.format(sample=sample_name, prefix='')
    plot(amino_counts, amino_plot_path)
    plot(nuc_counts, nuc_plot_path)
    # for (seed_name, consensus_name), counts in sorted(nuc_counts.iteritems()):
    #     for pos in sorted(counts.iterkeys()):
    #         print seed_name, consensus_name, pos, counts[pos]


def build_position_map(consensus_map, seed_map):
    gap_open = 15
    gap_extend = 5
    use_terminal_gap_penalty = 1
    position_map = defaultdict(dict)  # {seed_name: {consensus_name: seed_seq}}
    for source_region, source_seq in seed_map.iteritems():
        for dest_region, dest_seq in consensus_map.iteritems():
            aligned_source, aligned_dest, _score = gotoh.align_it_aa(
                source_seq,
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

            position_map[source_region][dest_region] = ''.join(seed_nucs)

    return position_map

main()
