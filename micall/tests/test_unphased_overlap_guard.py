"""Reference-free rejection of unphased long overlaps, using original reads."""
import random
import unittest
from collections import Counter
from io import StringIO

from Bio import SeqIO
from micall.utils.referenceless_contig_stitcher import (
    _cannot_phase_long_overlap, referenceless_contig_stitcher,
)


class LongUnphasedOverlapTests(unittest.TestCase):
    def test_boundary_and_disabled_behavior(self):
        for merged, rejected in ((2000, False), (1700, False),
                                 (1699, True), (1460, True)):
            with self.subTest(merged=merged):
                self.assertEqual(
                    _cannot_phase_long_overlap(1000, 1000, merged, 150, True),
                    rejected)
        self.assertFalse(_cannot_phase_long_overlap(1000, 1000, 1460, 150, False))

    def test_raw_reads_cannot_phase_a_540_base_overlap(self):
        rng = random.Random(472)
        molecule = ''.join(rng.choices('ACGT', k=1800))
        first, second = molecule[:1050], molecule[510:]
        self.assertEqual(first[-540:], second[:540])
        fasta = f'>left\n{first}\n>right\n{second}\n'
        reads = {150: Counter(molecule[i:i+150] for i in range(1651))}
        def assembled_lengths(evidence):
            out = StringIO()
            referenceless_contig_stitcher(StringIO(fasta), out,
                                          read_index=evidence,
                                          minimum_read_depth=1, read_length=150)
            return sorted(len(row.seq) for row in SeqIO.parse(
                StringIO(out.getvalue()), 'fasta'))
        self.assertEqual(assembled_lengths(None), [1800])
        self.assertEqual(assembled_lengths(reads), [1050, 1290])


if __name__ == '__main__':
    unittest.main()
