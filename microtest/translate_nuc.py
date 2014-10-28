import aln2counts

""" Translate nucleotide sequence to amino acid sequence and compare the result
with an expected sequence.
"""

nuc_seq = ''.join([
        "CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGC",
        "TCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAAC",
        "CAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAA",
        "ATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAG",
        "AAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT"
      ])
aa_compare = ''.join([
        "PQVTLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIE",
        "ICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF"
      ])

aa_seq = aln2counts.translate(nuc_seq)
pairs = zip(aa_seq, aa_compare)
diffs = [' ' if a == b else '*' for a, b in pairs]
print 'result ', aa_seq
print 'diffs  ', ''.join(diffs) if aa_seq != aa_compare else 'no diffs'
print 'compare', aa_compare
