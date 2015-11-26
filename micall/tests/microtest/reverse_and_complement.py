""" Reverse a nucleotide sequence and replace with complementary nucleotides.
Mixtures are allowed, as well as *, N, and -.
If you want to compare the result to an expected sequence, put the expected
sequence in reverse_compare.
Source: https://github.com/ArtPoon/bioinfo/blob/master/seqUtils.py#L143
"""
from micall.utils.translation import reverse_and_complement

nuc_seq = ''.join([
        "TGTACAAGACCCAACAACAATACAAGAAAAAGTATACATATAGGACCAGGGAGAGCATTTTATGC",
        "AACAGGAGAAATAATAGGAGATATAAGACAAGCACATTGT"
      ])
reverse_compare = ''.join([
        "ACAATGTGCTTGTCTTATATCTCCTATTATTTCTCCTGTTGCATAAAATGCTCTCCCTGGTCCTA",
        "TATGTATACTTTTTCTTGTATTGTTGTTGGGTCTTGTACA"
      ])


reverse_seq = reverse_and_complement(nuc_seq)
pairs = zip(reverse_seq, reverse_compare)
diffs = [' ' if a == b else '*' for a, b in pairs]
print 'result ', reverse_seq
print 'diffs  ', ''.join(diffs) if reverse_seq != reverse_compare else 'no diffs'
print 'compare', reverse_compare
