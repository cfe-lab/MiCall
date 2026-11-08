"""Read-pair-level phase evidence from exact truth-unique nucleotide motifs.

This procedure does not use MiCall's mapper. It detects a truth-specific
21-mer around each SNP private to a single verified clone, in either mate
and either strand, and counts joint observations within one original pair.

Note: the resulting counts are not unbiased clone frequencies: library
preparation, fragment coverage, errors, PCR duplicates, and omitted
unsequenced haplotypes can all distort them.
"""
import argparse
import gzip
import json
from pathlib import Path

from Bio import SeqIO


def reverse_complement(seq):
    return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def find_private_positions(references):
    lengths = {len(sequence) for sequence in references.values()}
    if len(lengths) != 1:
        raise ValueError('Clone references must be on an aligned coordinate axis')
    return {
        ref: [
            position for position, base in enumerate(seq)
            if all(base != other[position] for name, other in
                   references.items() if name != ref)
        ] for ref, seq in references.items()
    }


def find_private_motifs(references, width=21):
    if width % 2 != 1 or width < 3:
        raise ValueError('motif width must be odd and >= 3')
    radius = width // 2
    candidates = find_private_positions(references)
    motifs = {}
    for ref, positions in candidates.items():
        sequence = references[ref]
        result = []
        for p in positions:
            if p < radius or p + radius >= len(sequence):
                continue
            motif = sequence[p-radius:p+radius+1]
            opposite = reverse_complement(motif)
            if any(motif in alt or opposite in alt for name, alt in
                   references.items() if name != ref):
                continue
            result.append((p + 1, motif, opposite))
        motifs[ref] = result
    return motifs


def open_fastq(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    return opener(path, 'rt')


def sequences(read_file):
    while True:
        header = read_file.readline()
        if not header:
            return
        seq = read_file.readline().strip().upper()
        plus = read_file.readline()
        quality = read_file.readline()
        if (not header.startswith('@') or not plus.startswith('+')
                or len(seq) != len(quality.rstrip('\r\n'))):
            raise ValueError('Invalid paired FASTQ record')
        yield header.split()[0].removeprefix('@').removesuffix('/1').removesuffix('/2'), seq


def evaluate(references, r1_path, r2_path):
    motifs = find_private_motifs(references)
    results = {
        ref: {
            'private_positions': [pos for pos, _, _ in patterns],
            'at_least_one': 0, 'at_least_two': 0, 'at_least_three': 0,
            'all_marker_motifs': 0,
            'per_position_pairs': {str(pos): 0 for pos, _, _ in patterns},
        } for ref, patterns in motifs.items()
    }
    original_pairs = 0
    with open_fastq(r1_path) as fh1, open_fastq(r2_path) as fh2:
        from itertools import zip_longest
        for left, right in zip_longest(sequences(fh1), sequences(fh2)):
            if left is None or right is None or left[0] != right[0]:
                raise ValueError('Paired FASTQ streams are not synchronized')
            original_pairs += 1
            a, b = left[1], right[1]
            for ref, patterns in motifs.items():
                if not patterns:
                    continue
                hits = 0
                for pos, motif, rc in patterns:
                    if (motif in a or motif in b or rc in a or rc in b):
                        hits += 1
                        results[ref]['per_position_pairs'][str(pos)] += 1
                results[ref]['at_least_one'] += hits >= 1
                results[ref]['at_least_two'] += hits >= 2
                results[ref]['at_least_three'] += hits >= 3
                results[ref]['all_marker_motifs'] += hits == len(patterns)
    return {'original_pairs': original_pairs, 'motif_width': 21,
            'per_clone': results,
            'interpretation': ('Exact private-marker cooccurrence in original '
                               'read pairs; not calibrated clone abundance.')}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--truth', type=Path, required=True)
    parser.add_argument('--r1', type=Path, required=True)
    parser.add_argument('--r2', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    references = {record.id: str(record.seq).upper()
                  for record in SeqIO.parse(args.truth, 'fasta')}
    result = evaluate(references, args.r1, args.r2)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print('Examined', result['original_pairs'], 'paired reads')
    for ref, data in result['per_clone'].items():
        print(ref, data['at_least_one'], data['at_least_two'],
              data['at_least_three'], flush=True)


if __name__ == '__main__':
    main()
