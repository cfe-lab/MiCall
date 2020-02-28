import statistics
import re
import typing
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from csv import DictReader, DictWriter
from itertools import groupby
from operator import itemgetter
import logging
from collections import Counter

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
    columns = ['sample', 'contig', 'error', 'sequence', 'seqlen']
    for target_name in primer_targets:
        for column_type in [
            'probe_hxb2_start',
            'full_real_primer_seq',
            'in_probe_start',
            'in_probe_size',
            'in_hxb2_start',
            'in_hxb2_size',
            'is_reversed',
            'seq',
            'actual_primer_seq',
            'overhang',
            'dist',
            'error'
        ]:
            columns.append(target_name + '_' + column_type)
    non_tcga = re.compile(r'[^TCGA-]+')
    writer = DictWriter(probes_csv, columns)
    writer.writeheader()
    results = {}
    reader = DictReader(contigs_csv)
    projects = ProjectConfig.loadDefault()
    hxb2 = projects.getReference('HIV1-B-FR-K03455-seed')
    skipped = {}
    total = 0
    viable = 0
    unique_samples = 0
    for sample_name, sample_rows in groupby(reader, itemgetter('sample')):
        contig_num = 0
        unique_samples += 1
        for row in sample_rows:
            total += 1
            seed_name = row.get('genotype') or row.get('ref') or row['region']
            conseq_cutoff = row.get('consensus-percent-cutoff')
            contig_num += 1
            contig_name = f'{contig_num}-{seed_name}'
            uname = f'{sample_name}_{contig_name}_{contig_num}'
            new_row = dict(sample=sample_name, contig=contig_name)
            contig_seq: str = row.get('contig') or row['sequence']
            contig_seq = contig_seq.upper()
            new_row['seqlen'] = len(contig_seq)
            new_row['sequence'] = contig_seq
            if conseq_cutoff and conseq_cutoff != 'MAX':
                skipped[uname] = 'contig not MAX'
                new_row['error'] = skipped[uname]
                writer.writerow(new_row)
                continue
            found_non_tcga = re.findall(non_tcga, contig_seq)
            mixtures = len([x for x in found_non_tcga if x[0].upper() != 'X'])
            if (
                mixtures > 1
            ):
                skipped[uname] = 'contig sequence contained non-TCGA/gap'
                new_row['error'] = skipped[uname]
                writer.writerow(new_row)
                continue
            probelen = 100
            prime5_seq = contig_seq[:probelen]
            prime3_seq = contig_seq[-probelen:]
            gap_open_penalty = 15
            gap_extend_penalty = 3
            use_terminal_gap_penalty = 1
            for key in columns:
                if key not in ['sample', 'contig', 'seqlen', 'error', 'sequence']:
                    new_row[key] = None
            for end, seq in [(5, prime5_seq), (3, prime3_seq)]:
                if end == 5:
                    name = 'fwd_primer2'
                    hxb2_target_start = primer_targets[name]['hxb2_start']
                    hxb2_target_end = primer_targets[name]['hxb2_end'] + 100
                    hxb2_target_seq = hxb2[hxb2_target_start:hxb2_target_end]
                else:
                    name = 'rev_primer2'
                    hxb2_target_start = primer_targets[name]['hxb2_start'] - 100
                    hxb2_target_end = primer_targets[name]['hxb2_end']
                    hxb2_target_seq = hxb2[hxb2_target_start:hxb2_target_end]
                primer = None

                # DEBUG INTERESTING SAMPLE
                # interesting_sample = 'JLAT4CP-6-HIV_S105'
                # if (
                #     (sample_name == interesting_sample)
                #     & (seed_name == '1-1-HIV1-B-FR-K03455-seed')
                # ):
                #     import pdb; pdb.set_trace()

                if 'X' in seq:
                    seqlen = len(seq)
                    seq = handle_x(seq)
                    if not seq or len(seq) < seqlen / 6:
                        skipped[uname] = 'too many X in sequence'
                        new_row[prefix + 'error'] = skipped[uname]
                        continue
                finder = ProbeFinder(hxb2_target_seq, seq)
                finder.start += hxb2_target_start
                prefix = name + '_'
                new_row[prefix + 'probe_hxb2_start'] = finder.start
                new_row[prefix + 'full_real_primer_seq'] = primer_targets[name]['sequence']

                # If the segment overlaps the primer
                if primer_targets[name]['hxb2_start'] - probelen <= finder.start <= primer_targets[name]['hxb2_end']:
                    primer = validate_primer(finder, seq, primer_targets[name], hxb2_target_seq)
                    if primer['error']:
                        skipped[uname] = primer['error']
                        new_row[prefix + 'error'] = skipped[uname]
                else:
                    if primer_targets[name]['hxb2_start'] - probelen > finder.start:
                        skipped[uname] = f'{end} contig probe ends before hxb2 primer start'
                    elif finder.start > primer_targets[name]['hxb2_end']:
                        skipped[uname] = f'{end} contig probe starts after hxb2 primer end'
                    new_row[prefix + 'error'] = skipped[uname]
                    primer = validate_primer(finder, seq, primer_targets[name], hxb2_target_seq)
                    if primer['error']:
                        skipped[uname] = primer['error']
                        new_row[prefix + 'error'] += f' primer_error: {skipped[uname]}'
                new_row[prefix + 'in_probe_start'] = primer['seq_start']
                try:
                    new_row[prefix + 'in_probe_size'] = primer['seq_end'] - primer['seq_start']
                except TypeError:
                    new_row[prefix + 'in_probe_size'] = None
                new_row[prefix + 'in_hxb2_start'] = primer['hxb2_start']
                try:
                    new_row[prefix + 'in_hxb2_size'] = primer['hxb2_end'] - primer['hxb2_start']
                except TypeError:
                    new_row[prefix + 'in_hxb2_size'] = None
                new_row[prefix + 'is_reversed'] = ('Y'
                                                    if finder.is_reversed
                                                    else 'N')
                new_row[prefix + 'seq'] = primer['target_seq']
                new_row[prefix + 'overhang'] = primer['overhang']
                new_row[prefix + 'dist'] = primer['dist']
                new_row[prefix + 'actual_primer_seq'] = primer['real_primer']
            writer.writerow(new_row)
            viable += 1
    with open(f'{probes_csv.name}.counts.csv', 'w') as o:
        counts = Counter(skipped.values())
        counts['total'] = total
        counts['viable'] = viable
        counts['unique_samples'] = unique_samples
        order = [
            'total',
            'contig not MAX',
            'unique_samples',
            'contig sequence contained non-TCGA/gap',
            'too many X in sequence',
            'no 5 primer sequence found',
            'no 3 primer sequence found',
            'real primer not found at expected coordinates',
            'primer in contig sequence not found at expected coordinates',
            'mismatches in primer > tolerance',
            '5 contig probe ends before hxb2 primer start',
            '3 contig probe ends before hxb2 primer start',
            '5 contig probe starts after hxb2 primer end',
            '3 contig probe starts after hxb2 primer end',
            'viable',
            'error'
        ]
        o.write(','.join(order) + '\n')
        o.write(','.join([str(counts[k]) for k in order]))
        print(counts['viable']/counts['unique_samples'])
    return results

def handle_x(sequence):
    length = len(sequence)
    midpoint = length / 2
    x_positions = [i for i,j in enumerate(sequence) if j == 'X']
    try:
        rightmost_x = max([x for x in x_positions if x <= midpoint])
    except ValueError:
        rightmost_x = None
    try:
        leftmost_x = min([x for x in x_positions if x > midpoint])
    except ValueError:
        leftmost_x = None
    if rightmost_x and leftmost_x:
        sequence = sequence[rightmost_x+1:leftmost_x]
    elif rightmost_x:
        sequence = sequence[rightmost_x+1:]
    else:
        sequence = sequence[:leftmost_x]
    return sequence

def validate_primer(finder, finder_seq, target, hxb2_target, tolerance=1):
    if finder.is_reversed:
        finder_seq = reverse_and_complement(finder_seq)
    error = None
    primer_in_finder = None
    finder_seqsize = len(finder_seq)
    matched_finder_size = len(finder.contig_match)
    overhang = matched_finder_size - finder_seqsize
    primer_in_finder_hxb2_start_coord = max(finder.start, target['hxb2_start'])
    primer_in_finder_hxb2_end_coord = min(
        finder.start + finder_seqsize,
        target['hxb2_end']
    )
    # Get the coordinates of the primer relative to the finder sequence
    primer_in_finder_start_coord = primer_in_finder_hxb2_start_coord - finder.start
    primer_in_finder_end_coord = primer_in_finder_hxb2_end_coord - finder.start
    if (
        target['hxb2_start'] == 9603
        & primer_in_finder_end_coord > matched_finder_size
    ):
        primer_in_finder_start_coord -= overhang
        primer_in_finder_end_coord -= overhang

    # Get the primer sequence of the finder sequence
    primer_in_finder = finder_seq[primer_in_finder_start_coord:primer_in_finder_end_coord]

    # Get the sequence of the true primer that overlaps the finder sequence
    primer_start_coord = max(
        0,
        primer_in_finder_hxb2_start_coord - target['hxb2_start']
    )
    primer_end_coord = primer_start_coord + len(primer_in_finder)
    real_primer = target['sequence'][primer_start_coord:primer_end_coord]
    result = {
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
        'dist': 0,
        'real_primer': real_primer,
        'error': error
    }

    if len(real_primer) != len(primer_in_finder):
        result['error'] = 'Real primer did not match length of primer in finder'
        return result
    if not real_primer:
        result['error'] = 'real primer not found at expected coordinates'
        return result
    elif not primer_in_finder:
        result['error'] = 'primer in contig sequence not found at expected coordinates'
        return result
    for i in range(len(real_primer)):
        try:
            our_nuc = primer_in_finder[i]
        except IndexError:
            pass
        real_nuc = real_primer[i]
        if our_nuc != real_nuc:
            if real_nuc in mixture_dict and our_nuc in mixture_dict[real_nuc]:
                pass
            else:
                result['dist'] += 1
        if result['dist'] > tolerance:
            result['error'] = 'mismatches in primer > tolerance'
    return result


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
    for mixture in seq.upper():
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
