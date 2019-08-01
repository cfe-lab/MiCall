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
        return self.get_consensus_for_length(seq_length)

    def get_consensus_for_length(self, length):
        nucleotides = self.length_nucleotides[length]
        return ''.join(nucleotides[i].get_consensus(MAX_CUTOFF)
                       for i in range(length))

    def get_consensus_by_lengths(self):
        if not self.length_counts:
            return

        window_radius = 20  # number of neighbours on either side to compare
        min_ratio = 20  # spike must be this many times higher than average
        window_size = window_radius * 2 + 1
        window_totals = []
        running_total = 0
        max_length = max(self.length_counts)
        for length in range(0, max_length+window_radius+1):
            length_count = self.length_counts[length]
            running_total += length_count
            window_totals.append(length_count)
            if len(window_totals) == window_size:
                centre_length = length - window_radius
                centre_count = self.length_counts[centre_length]
                average_count = max(
                    1,
                    (running_total - centre_count) / (window_size-1))
                if centre_count and (centre_count >= average_count * min_ratio):
                    yield self.get_consensus_for_length(centre_length)
                running_total -= window_totals.pop(0)
