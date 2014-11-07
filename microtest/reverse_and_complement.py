""" Reverse a nucleotide sequence and replace with complementary nucleotides.
Mixtures are allowed, as well as *, N, and -.
If you want to compare the result to an expected sequence, put the expected
sequence in reverse_compare.
Source: https://github.com/ArtPoon/bioinfo/blob/master/seqUtils.py#L143
"""

nuc_seq = ''.join([
        "TGTACAAGACCCAACAACAATACAAGAAAAAGTATACATATAGGACCAGGGAGAGCATTTTATGC",
        "AACAGGAGAAATAATAGGAGATATAAGACAAGCACATTGT"
      ])
reverse_compare = ''.join([
        "ACAATGTGCTTGTCTTATATCTCCTATTATTTCTCCTGTTGCATAAAATGCTCTCCCTGGTCCTA",
        "TATGTATACTTTTTCTTGTATTGTTGTTGGGTCTTGTACA"
      ])

complement_dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 
                    'W':'S', 'R':'Y', 'K':'M', 'Y':'R', 'S':'W', 'M':'K',
                    'B':'V', 'D':'H', 'H':'D', 'V':'B',
                    '*':'*', 'N':'N', '-':'-'}

def reverse_and_complement(seq):
    rseq = seq[::-1]
    rcseq = ''
    for i in rseq:  # reverse order
        rcseq += complement_dict[i]
    return rcseq

reverse_seq = reverse_and_complement(nuc_seq)
pairs = zip(reverse_seq, reverse_compare)
diffs = [' ' if a == b else '*' for a, b in pairs]
print 'result ', reverse_seq
print 'diffs  ', ''.join(diffs) if reverse_seq != reverse_compare else 'no diffs'
print 'compare', reverse_compare
