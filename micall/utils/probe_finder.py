import re
import typing
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from csv import DictReader, DictWriter
from itertools import groupby
from operator import itemgetter
import logging

import Levenshtein
from gotoh import align_it

from micall.core.project_config import ProjectConfig
from micall.utils.translation import mixture_dict, reverse_and_complement

logger = logging.getLogger('micall')

TARGET_SEQUENCES = dict(
    gag_probe='ACTGGTGAGTACGCCAAAA',
    env_probe='CCTTGGGTTGTTGGGA',
    round2_fwd_primer='GCGCCCGAACAGGGACYTGAAARCGAAAG',
    round2_rev_primer='TAAGCCTCAATAAAGCTTGCCTTGAGTGC')


def parse_args():
    parser = ArgumentParser(
        description='Search contigs.csv for known probe and primer sequences.',
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('contigs_csv',
                        help='CSV file with contig sequences or consensus sequences',
                        nargs='?',
                        default='contigs.csv',
                        type=FileType())
    parser.add_argument('probes_csv',
                        help='CSV file to write probe sequences in',
                        nargs='?',
                        default='probes.csv',
                        type=FileType('w'))
    parser.add_argument('-t',
                        '--test_find_primers',
                        action='store_true'
    )
    return parser.parse_args()


def find_primers(contigs_csv, probes_csv):
    # Note these are 1-based indicies
    primer_targets = {
        'fwd_primer2': {
            'sequence': 'GCGCCCGAACAGGGACYTGAAARCGAAAG',
            # convert to 0-base index
            'hxb2_start': 638 - 1,
            'hxb2_end': 666
        },
        'rev_primer2': {
            'sequence': 'TAAGCCTCAATAAAGCTTGCCTTGAGTGC',
            # convert to 0-base index
            'hxb2_start': 9604 - 1,
            'hxb2_end': 9632
        }
    }
    columns = ['sample', 'contig']
    for target_name in primer_targets:
        for column_type in [
            # 'contig_probe_seq',
            # 'contig_match',
            'contig_probe_hxb2_start',
            'in_contig_start',
            'in_contig_size',
            'in_hxb2_start',
            'in_hxb2_size',
            'is_reversed',
            'seq',
            'actual_primer_seq',
            'overhang',
            'dist'
        ]:
            columns.append(target_name + '_' + column_type)
    writer = DictWriter(probes_csv, columns)
    writer.writeheader()
    results = {}
    reader = DictReader(contigs_csv)
    projects = ProjectConfig.loadDefault()
    hxb2 = projects.getReference('HIV1-B-FR-K03455-seed')
    for sample_name, sample_rows in groupby(reader, itemgetter('sample')):
        contig_num = 0
        for row in sample_rows:
            seed_name = row.get('genotype') or row.get('ref') or row['region']
            conseq_cutoff = row.get('consensus-percent-cutoff')
            if conseq_cutoff and conseq_cutoff != 'MAX':
                continue
            contig_num += 1
            contig_name = f'{contig_num}-{seed_name}'
            contig_seq: str = row.get('contig') or row['sequence']
            prime5_seq = contig_seq[:100]
            prime3_seq = contig_seq[-100:]
            gap_open_penalty = 15
            gap_extend_penalty = 3
            use_terminal_gap_penalty = 1
            new_row = dict(sample=sample_name, contig=contig_name)
            skip = False
            for end, seq in [(5, prime5_seq), (3, prime3_seq)]:
                primer = None
                if 'X' in seq.upper():
                    skip = True
                    break
                finder = ProbeFinder(hxb2, seq)
                if end == 5:
                    name = 'fwd_primer2'
                else:
                    name = 'rev_primer2'
                # If the segment overlaps the primer
                if finder.start <= primer_targets[name]['hxb2_end']:
                    primer = validate_primer(finder, seq, primer_targets[name])
                if not primer:
                    skip = True
                    break
                prefix = name + '_'
                # new_row[prefix + 'contig_probe_seq'] = primer['finder_seq']
                new_row[prefix + 'contig_probe_hxb2_start'] = primer['contig_hxb2_start']
                new_row[prefix + 'in_contig_start'] = primer['seq_start']
                new_row[prefix + 'in_contig_size'] = primer['seq_end'] - primer['seq_start']
                new_row[prefix + 'in_hxb2_start'] = primer['hxb2_start']
                new_row[prefix + 'in_hxb2_size'] = primer['hxb2_end'] - primer['hxb2_start']
                new_row[prefix + 'is_reversed'] = ('Y'
                                                    if finder.is_reversed
                                                    else 'N')
                new_row[prefix + 'seq'] = primer['target_seq']
                new_row[prefix + 'overhang'] = primer['overhang']
                new_row[prefix + 'dist'] = primer['dist']
                new_row[prefix + 'actual_primer_seq'] = primer['real_primer']
                # new_row[prefix + 'contig_match'] = primer['contig_match']
            if not skip:
                writer.writerow(new_row)
    return results


def validate_primer(finder, finder_seq, target, tolerance=1):
    projects = ProjectConfig.loadDefault()
    hxb2 = projects.getReference('HIV1-B-FR-K03455-seed')
    matched_finder_size = len(finder.contig_match)
    finder_seqsize = len(finder_seq)
    overhang = matched_finder_size - finder_seqsize
    primer_in_finder_hxb2_start_coord = max(finder.start, target['hxb2_start'])
    primer_in_finder_hxb2_end_coord = min(
        finder.start + finder_seqsize,
        target['hxb2_end']
    )
    # Get the coordinates of the primer relative to the finder sequence
    primer_in_finder_start_coord = primer_in_finder_hxb2_start_coord - finder.start
    primer_in_finder_end_coord = primer_in_finder_hxb2_end_coord - finder.start
    if target['hxb2_start'] == 9603:
        primer_in_finder_start_coord -= overhang
        primer_in_finder_end_coord -= overhang
    primer_in_finder = finder_seq[primer_in_finder_start_coord:primer_in_finder_end_coord]

    # Get the sequence of the true primer that overlaps the finder sequence
    primer_start_coord = max(
        0,
        primer_in_finder_hxb2_start_coord - target['hxb2_start']
    )

    primer_end_coord = primer_start_coord + primer_in_finder_hxb2_end_coord - primer_in_finder_hxb2_start_coord

    real_primer = target['sequence'][primer_start_coord:primer_end_coord]
    mismatches = 0
    if len(real_primer) != len(primer_in_finder):
        logger.warning('Real primer did not match length of primer in finder!')
        logger.warning(f'Contig seq: {finder.contig_match}')
        logger.warning(f'Contig start in hxb2: {finder.start}')
    if not real_primer or not primer_in_finder:
        return None
    for i in range(len(real_primer)):
        our_nuc = primer_in_finder[i]
        real_nuc = real_primer[i]
        if our_nuc != real_nuc:
            if real_nuc in mixture_dict and our_nuc in mixture_dict[real_nuc]:
                pass
            else:
                mismatches += 1
        if mismatches > tolerance:
            return None
    return {
        'finder_seq': finder_seq,
        'target_seq': primer_in_finder,
        'contig_hxb2_start': finder.start,
        'contig_hxb2_end': matched_finder_size,
        'seq_start': primer_in_finder_start_coord,
        'seq_end': primer_in_finder_end_coord,
        'hxb2_start': primer_in_finder_hxb2_start_coord,
        'hxb2_end': primer_in_finder_hxb2_end_coord,
        'overhang': overhang,
        'contig_match': finder.contig_match,
        # 'dist': Levenshtein.distance(primer_in_finder, real_primer),
        'dist': mismatches,
        'real_primer': real_primer
    }


def find_probes(contigs_csv, probes_csv):
    reader = DictReader(contigs_csv)
    columns = ['sample', 'contig']
    for target_name in TARGET_SEQUENCES:
        for column_type in ['in_contig_start',
                            'in_contig_size',
                            'in_hxb2_start',
                            'in_hxb2_size',
                            'merged_hxb2_start',
                            'merged_hxb2_size',
                            'dist',
                            'end_dist',
                            'score',
                            'is_reversed',
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
            seed_name = row.get('genotype') or row.get('ref') or row['region']
            conseq_cutoff = row.get('consensus-percent-cutoff')
            if conseq_cutoff and conseq_cutoff != 'MAX':
                continue
            contig_num += 1
            contig_name = f'{contig_num}-{seed_name}'
            contig_seq: str = row.get('contig') or row['sequence']
            aligned_hxb2, aligned_contig_to_hxb2, _ = align_it(
                hxb2,
                contig_seq,
                gap_open_penalty,
                gap_extend_penalty,
                use_terminal_gap_penalty)
            new_row = dict(sample=sample_name, contig=contig_name)
            for target_name, target_seq in TARGET_SEQUENCES.items():
                finder = ProbeFinder(contig_seq, target_seq)

                size = len(finder.contig_match)
                start_pos = finder.start + 1
                end_pos = finder.start + size
                hxb2_pos = contig_pos = 0
                merged_hxb2_start = merged_hxb2_size = None
                for hxb2_nuc, contig_nuc in zip(aligned_hxb2,
                                                aligned_contig_to_hxb2):
                    if hxb2_nuc != '-':
                        hxb2_pos += 1
                    if contig_nuc != '-':
                        contig_pos += 1
                        if contig_pos == start_pos:
                            merged_hxb2_start = hxb2_pos
                        if contig_pos == end_pos:
                            merged_hxb2_size = hxb2_pos - merged_hxb2_start + 1
                            break

                aligned_ref, aligned_match, _ = align_it(
                    hxb2,
                    finder.contig_match,
                    gap_open_penalty,
                    gap_extend_penalty,
                    use_terminal_gap_penalty)
                lstripped_match = aligned_match.lstrip('-')
                in_hxb2_start = len(aligned_match) - len(lstripped_match)
                tail_len = len(lstripped_match) - len(lstripped_match.rstrip('-'))
                ref_match = aligned_ref[in_hxb2_start:-tail_len or None]
                in_hxb2_size = len(ref_match.replace('-', ''))

                prefix = target_name + '_'
                new_row[prefix + 'in_contig_start'] = start_pos
                new_row[prefix + 'in_contig_size'] = size
                new_row[prefix + 'in_hxb2_start'] = in_hxb2_start
                new_row[prefix + 'in_hxb2_size'] = in_hxb2_size
                new_row[prefix + 'merged_hxb2_start'] = merged_hxb2_start
                new_row[prefix + 'merged_hxb2_size'] = merged_hxb2_size
                new_row[prefix + 'dist'] = finder.dist
                new_row[prefix + 'end_dist'] = finder.end_dist
                new_row[prefix + 'score'] = finder.score
                new_row[prefix + 'is_reversed'] = ('Y'
                                                   if finder.is_reversed
                                                   else 'N')
                new_row[prefix + 'seq'] = finder.contig_match
            writer.writerow(new_row)


def unpack_mixtures_and_reverse(seq: str) -> typing.Set[typing.Tuple[str, bool]]:
    """ Unpack mixture nucleotide codes, and add reverse complements.

    :param seq: nucleotide sequence, possibly including mixture codes
    :return: unpacked and reversed sequences, along with is_reversed flag
    """
    old_mixtures = {''}
    for mixture in seq:
        new_mixtures = set()
        for nuc in mixture_dict.get(mixture, mixture):
            for old_mixture in old_mixtures:
                new_mixtures.add(old_mixture + nuc)
        old_mixtures = new_mixtures
    forward_results = {(mixture, False)
                       for mixture in old_mixtures}
    reversed_results = {(reverse_and_complement(mixture), True)
                        for mixture in old_mixtures}
    return forward_results | reversed_results


class ProbeFinder:
    def __init__(self, contig_seq: str, target_seq: str):
        gap_open_penalty = 15
        gap_extend_penalty = 3
        use_terminal_gap_penalty = 1
        best_acontig = best_atarget = best_target = best_score = None
        best_reversed = None
        for target_nucs, is_reversed in unpack_mixtures_and_reverse(
                target_seq):
            aligned_contig, aligned_target, score = align_it(
                contig_seq,
                target_nucs,
                gap_open_penalty,
                gap_extend_penalty,
                use_terminal_gap_penalty)
            if best_score is None or score > best_score:
                best_acontig = aligned_contig
                best_atarget = aligned_target
                best_target = target_nucs
                best_score = score
                best_reversed = is_reversed
        aligned_contig = best_acontig
        aligned_target = best_atarget
        target_nucs = best_target
        self.score = best_score
        self.is_reversed = best_reversed
        if self.is_reversed:
            aligned_contig = reverse_and_complement(aligned_contig)
            aligned_target = reverse_and_complement(aligned_target)
        match = re.match('-*([^-](.*[^-])?)', aligned_target)
        self.start = match.start(1)
        end = match.end(1)
        self.contig_match = aligned_contig[self.start:end].replace('-', '')
        self.dist = Levenshtein.distance(target_nucs, self.contig_match)
        stripped_contig = aligned_contig.lstrip('-')
        overhang = len(aligned_contig) - len(stripped_contig)
        if overhang > 0:
            stripped_target = target_nucs[overhang:]
            self.end_dist = Levenshtein.distance(stripped_target,
                                                 self.contig_match)
        else:
            stripped_contig = aligned_contig.rstrip('-')
            overhang = len(aligned_contig) - len(stripped_contig)
            if overhang == 0:
                self.end_dist = self.dist
            else:
                stripped_target = target_nucs[:-overhang]
                self.end_dist = Levenshtein.distance(stripped_target,
                                                     self.contig_match)


def main():
    args = parse_args()
    if args.test_find_primers:
        find_primers(args.contigs_csv, args.probes_csv)
    else:
        find_probes(args.contigs_csv, args.probes_csv)


if __name__ in ('__main__', '__live_coding__'):
    main()
