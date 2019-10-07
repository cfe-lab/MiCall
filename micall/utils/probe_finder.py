import re
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from csv import DictReader, DictWriter
from itertools import groupby
from operator import itemgetter

import Levenshtein
from gotoh import align_it

from micall.core.project_config import ProjectConfig

TARGET_SEQUENCES = dict(
    gag_probe='TTTTGGCGTACTCACCAGT',
    env_probe='CCTTGGGTTGTTGGGA',
    round2_fwd_primer='GCGCCCGAACAGGGACYTGAAARCGAAAG',
    round2_rev_primer='TAAGCCTCAATAAAGCTTGCCTTGAGTGC')


def parse_args():
    parser = ArgumentParser(
        description='Search contigs.csv for known probe and primer sequences.',
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('contigs_csv',
                        help='CSV file with contig sequences',
                        nargs='?',
                        default='contigs.csv',
                        type=FileType())
    parser.add_argument('probes_csv',
                        help='CSV file to write probe sequences in',
                        nargs='?',
                        default='probes.csv',
                        type=FileType('w'))
    return parser.parse_args()


def find_probes(contigs_csv, probes_csv):
    reader = DictReader(contigs_csv)
    columns = ['sample', 'contig']
    for target_name in TARGET_SEQUENCES:
        for column_type in ['start',
                            'size',
                            'hxb2_start',
                            'hxb2_size',
                            'dist',
                            'score',
                            'seq']:
            columns.append(target_name + '_' + column_type)
    writer = DictWriter(probes_csv, columns)
    writer.writeheader()
    projects = ProjectConfig.loadDefault()
    hxb2 = projects.getReference('HIV1-B-FR-K03455-seed')
    gap_open_penalty = 15
    gap_extend_penalty = 3
    use_terminal_gap_penalty = 1
    for sample_name, sample_rows in groupby(reader, itemgetter('sample')):
        contig_num = 0
        for row in sample_rows:
            seed_name = row['genotype']
            contig_num += 1
            contig_name = f'{contig_num}-{seed_name}'
            contig_seq: str = row['contig']
            aligned_hxb2, aligned_contig_to_hxb2, _ = align_it(
                hxb2,
                contig_seq,
                gap_open_penalty,
                gap_extend_penalty,
                use_terminal_gap_penalty)
            new_row = dict(sample=sample_name, contig=contig_name)
            for target_name, target_seq in TARGET_SEQUENCES.items():
                aligned_contig, aligned_target, score = align_it(contig_seq,
                                                                 target_seq,
                                                                 gap_open_penalty,
                                                                 gap_extend_penalty,
                                                                 use_terminal_gap_penalty)
                match = re.match('-*([^-](.*[^-])?)', aligned_target)
                start = match.start(1)
                end = match.end(1)
                contig_match = aligned_contig[start:end].replace('-', '')
                size = len(contig_match)
                dist = Levenshtein.distance(target_seq, contig_match)

                start_pos = start + 1
                end_pos = start + size
                hxb2_pos = contig_pos = 0
                hxb2_start = hxb2_size = None
                for hxb2_nuc, contig_nuc in zip(aligned_hxb2,
                                                aligned_contig_to_hxb2):
                    if hxb2_nuc != '-':
                        hxb2_pos += 1
                    if contig_nuc != '-':
                        contig_pos += 1
                        if contig_pos == start_pos:
                            hxb2_start = hxb2_pos
                        if contig_pos == end_pos:
                            hxb2_size = hxb2_pos - hxb2_start + 1
                            break

                prefix = target_name + '_'
                new_row[prefix + 'start'] = start_pos
                new_row[prefix + 'size'] = size
                new_row[prefix + 'hxb2_start'] = hxb2_start
                new_row[prefix + 'hxb2_size'] = hxb2_size
                new_row[prefix + 'dist'] = dist
                new_row[prefix + 'score'] = score
                new_row[prefix + 'seq'] = contig_match
            writer.writerow(new_row)


def main():
    args = parse_args()
    find_probes(args.contigs_csv, args.probes_csv)


main()
