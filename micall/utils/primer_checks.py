from pathlib import Path

# noinspection PyPackageRequirements
from Bio import SeqIO


def print_duplicates(primers, label):
    seen = set()
    for primer in primers:
        name = primer.name
        if name in seen:
            print(f'{label} duplicate: {name}')
        seen.add(name)


def main():
    data_path = Path(__file__).parent.parent / 'data'
    for stem in ('primers_hivgha_left', 'primers_hivgha_right'):
        print(f'###{stem}')
        source_path = data_path / (stem + '.fasta')
        target_path = data_path / (stem + '_end.fasta')
        source_primers = list(SeqIO.parse(source_path, 'fasta'))
        target_primers = list(SeqIO.parse(target_path, 'fasta'))
        targets = {primer.name: str(primer.seq)
                   for primer in target_primers}
        print_duplicates(source_primers, 'Start')
        print_duplicates(target_primers, 'End')
        for source_primer in source_primers:
            expected_target = str(source_primer.reverse_complement().seq)
            target_seq = targets.get(source_primer.name)
            if target_seq is None:
                print(f'>{source_primer.name}\n{expected_target}')
            elif target_seq != expected_target:
                print(f'>{source_primer.name}\n-{expected_target}\n+{target_seq}')




main()
