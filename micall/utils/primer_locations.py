from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, OPTIONAL
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO

from micall.core.project_config import ProjectConfig
from micall.utils.alignment_wrapper import align_nucs

PRIMER_SETS = {'HIVB': 'HIV1-B-FR-K03455-seed',
               'HIV-AD': 'HIV1-B-FR-K03455-seed',
               'HIVGHA': 'HIV1-B-FR-K03455-seed'}


@dataclass(order=True)
class PrimerLocation:
    start_pos: int
    end_pos: int
    name: str

    def __str__(self):
        return f'{self.start_pos}-{self.end_pos}: {self.name}'


def parse_args():
    parser = ArgumentParser(
        description="Find reference locations of a set of primers.",
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('project_code',
                        choices=sorted(PRIMER_SETS),
                        nargs=OPTIONAL,
                        default='HIVGHA',
                        help='Primer set to find.')
    return parser.parse_args()


def find_primer(seq: str, name: str, ref: str) -> PrimerLocation:
    aref, aseq, score = align_nucs(ref, seq)
    ltrimmed = aseq.lstrip('-')
    start_index = len(aseq) - len(ltrimmed)
    rtrimmed = ltrimmed.rstrip('-')
    aligned_length = len(rtrimmed)
    aref_section = aref[start_index:start_index+aligned_length]
    ref_section = aref_section.replace('-', '')
    end_index = start_index + len(ref_section)  # exclusive
    end_pos = end_index  # inclusive
    start_pos = start_index + 1
    return PrimerLocation(start_pos, end_pos, name)


def main():
    args = parse_args()
    project_code = args.project_code.upper()
    projects = ProjectConfig.loadDefault()
    ref = projects.getReference(PRIMER_SETS[project_code])
    print(f'Primers for {project_code}:')

    project_code = project_code.lower()

    primers_path = Path(__file__).parent.parent / 'data'
    left_path = primers_path / f'primers_{project_code}_left.fasta'
    right_path = primers_path / f'primers_{project_code}_right_end.fasta'

    locations = []
    with left_path.open() as left_fasta:
        for primer in SeqIO.parse(left_fasta, 'fasta'):
            assert primer.seq[0] == 'X'
            trimmed_seq = str(primer.seq[1:])
            locations.append(find_primer(trimmed_seq, primer.name, ref))
    with right_path.open() as right_fasta:
        for primer in SeqIO.parse(right_fasta, 'fasta'):
            assert primer.seq[-1] == 'X'
            trimmed_seq = str(primer.seq[:-1])
            locations.append(find_primer(trimmed_seq, primer.name, ref))
    locations.sort()
    print(*locations, sep='\n')

if __name__ == '__main__':
    main()
