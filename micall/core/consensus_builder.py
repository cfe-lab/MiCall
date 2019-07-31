from collections import Counter, defaultdict

from micall.core.aln2counts import SeedNucleotide, MAX_CUTOFF


class ConsensusBuilder:
    def __init__(self):
        self.length_counts = Counter()
        self.length_nucleotides = defaultdict(
            lambda: defaultdict(SeedNucleotide))
        self.nucleotides = defaultdict(SeedNucleotide)

    def build(self, merged_reads):
        for read in merged_reads:
            yield read
            merged_seq = read[3]
            if merged_seq is None:
                # Did not merge.
                continue
            seq_length = len(merged_seq)
            self.length_counts[seq_length] += 1
            nucleotides = self.length_nucleotides[seq_length]
            for i, nuc in enumerate(merged_seq):
                seed_nucleotide = nucleotides[i]
                seed_nucleotide.count_nucleotides(nuc)

    def get_consensus(self):
        seq_length = self.length_counts.most_common(1)[0][0]
        nucleotides = self.length_nucleotides[seq_length]
        return ''.join(nucleotides[i].get_consensus(MAX_CUTOFF)
                       for i in range(seq_length))
