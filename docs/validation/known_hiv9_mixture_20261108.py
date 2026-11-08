"""Known-source HIV9 minority challenge from independently Sanger-verified clones.

Construct paired-end reads from clone references with EXACT input proportions.
The sequences are simulated and error-free, not clinical library molecules.
Do not confuse demonstrated synthetic behavior with clinical sensitivity.
"""
import argparse
import gzip
import json
from pathlib import Path
import random

from Bio import SeqIO

RECORD_LENGTH = 250
INSERT_LENGTH = 450
WEIGHTS = {
    '0': [502, 250, 125, 63, 32, 16, 8, 4, 0],
    '0p2': [500, 250, 125, 63, 32, 16, 8, 4, 2],
    '1': [492, 250, 125, 63, 32, 16, 8, 4, 10],
}


def reverse_complement(text):
    return text.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def generate(truth, folder, percentage, count, seed):
    if count % 1000:
        raise ValueError('Total input pairs must be divisible by 1000')
    weights = WEIGHTS[percentage]
    if sum(weights) != 1000:
        raise ValueError('Input distribution must sum to 100%')
    originals = list(SeqIO.parse(truth, 'fasta'))
    if len(originals) != 9:
        raise ValueError('Expected nine externally validated reference clones')
    lengths = [len(record.seq) for record in originals]
    if min(lengths) < INSERT_LENGTH:
        raise ValueError('Reference sequence too short for fragment simulation')
    rng = random.Random(seed)
    clone_of_pair = [idx for idx, weight in enumerate(weights)
                     for _ in range(weight * count // 1000)]
    rng.shuffle(clone_of_pair)
    folder.mkdir(parents=True, exist_ok=False)
    first = folder / f'{folder.name}_L001_R1_001.fastq.gz'
    second = folder / f'{folder.name}_L001_R2_001.fastq.gz'
    with gzip.open(first, 'wt', compresslevel=3) as a, \
         gzip.open(second, 'wt', compresslevel=3) as b:
        for num, clone in enumerate(clone_of_pair):
            template = str(originals[clone].seq).upper()
            start = rng.randrange(0, len(template)-INSERT_LENGTH+1)
            fragment = template[start:start+INSERT_LENGTH]
            one = fragment[:RECORD_LENGTH]
            two = reverse_complement(fragment[-RECORD_LENGTH:])
            header = f'@clone{clone+1}_pair{num:08d}\n'
            for destination, read in ((a, one), (b, two)):
                destination.write(
                    header + read + '\n+\n' + 'I' * RECORD_LENGTH + '\n')
    summary = {
        'source': 'Sanger-confirmed HIV9 clone reference sequences',
        'method': 'error-free controlled synthetic fragments, no PCR biases',
        'input_pairs': count,
        'clone_pairs': {
            f'clone{i+1}': weight * count // 1000
            for i, weight in enumerate(weights)
        },
        'nominal_absent_or_minority_percentage': percentage,
        'rng_seed': seed,
        'read_length': RECORD_LENGTH,
        'insert_length': INSERT_LENGTH,
        'limitations': [
            'Synthetic, error-free input; not independent clinical validation',
            'Exactly matching reference seeds create a favorable condition',
            'No PCR duplication, cross-species contamination, or real phasing',
        ],
    }
    (folder/'synthetic_manifest.json').write_text(
        json.dumps(summary, indent=2)+'\n')
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--truth', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--minority', choices=list(WEIGHTS), required=True)
    parser.add_argument('--pairs', type=int, default=20000)
    parser.add_argument('--seed', type=int, default=441)
    args = parser.parse_args()
    print(json.dumps(generate(args.truth, args.output,
                              args.minority, args.pairs, args.seed),
                     indent=2))


if __name__ == '__main__':
    main()
