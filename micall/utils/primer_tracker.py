import typing
from pathlib import Path

from Bio import SeqIO

from micall.utils.alignment_wrapper import align_nucs


class PrimerTracker:
    def __init__(self, conseq: str, seed_name: str):
        self.conseq = conseq
        self.seed_name = seed_name
        self.ignored_positions: typing.Optional[typing.Set[int]] = None

    def is_ignored(self, pos: int):
        """ Check if a position should be ignored because it's in a primer.

        :param pos: the one-based nucleotide position within the full seed
            sequence
        :return: True if that position should be ignored.
        """
        self.check_positions()
        return pos in self.ignored_positions

    def check_positions(self):
        if self.ignored_positions is not None:
            return
        self.ignored_positions = set()
        if not self.seed_name.startswith('HCV-'):
            return
        cleaned_conseq = self.conseq.replace('-', 'x')
        data_path = Path(__file__).parent.parent / 'data'
        left_primers_path = data_path / f'primers_hcv_left.fasta'
        right_primers_path = data_path / f'primers_hcv_right_end.fasta'
        for primers_path in (left_primers_path, right_primers_path):
            with primers_path.open() as f:
                for primer in SeqIO.parse(f, 'fasta'):
                    primer_name = primer.name
                    if not primer_name.startswith('HCV'):
                        continue
                    if 'dA20' in primer_name or 'TIM' in primer_name:
                        continue
                    primer_seq = str(primer.seq).replace('X', '')
                    aligned_primer, aligned_conseq, score = align_nucs(
                        primer_seq,
                        cleaned_conseq)
                    primer_start = aligned_primer.lstrip('-')
                    start = len(aligned_primer) - len(primer_start)
                    primer_end = aligned_primer.rstrip('-')
                    padded_end = len(primer_end)
                    conseq_match = aligned_conseq[start:padded_end]
                    unpadded_match = conseq_match.replace('-', '')
                    end = start + len(unpadded_match)
                    self.ignored_positions.update(range(start+1, end+1))
