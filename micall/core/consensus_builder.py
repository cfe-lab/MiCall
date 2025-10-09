from collections import Counter, defaultdict

from micall.utils.report_amino import SeedNucleotide, FIRST_CUTOFF


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
        # Many assemblers (such as IVA) can't handle seeds with mixtures, so always avoid them.
        return ''.join(nucleotides[i].get_consensus(FIRST_CUTOFF)
                       for i in range(length))

    def get_consensus_by_lengths(self):
        if not self.length_counts:
            return

        window_radius = 20  # number of neighbours on either side to compare
        min_ratio = 20  # spike must be this many times higher than 3rd quartile
        window_size = window_radius * 2 + 1
        window_totals = []
        max_length = max(self.length_counts)
        for length in range(0, max_length+window_radius+1):
            length_count = self.length_counts[length]
            window_totals.append(length_count)
            if len(window_totals) == window_size:
                centre_length = length - window_radius
                centre_count = self.length_counts[centre_length]
                neighbour1_count = self.length_counts[centre_length - 1]
                neighbour2_count = self.length_counts[centre_length + 1]
                is_local_max = neighbour1_count <= centre_count > neighbour2_count
                if is_local_max:
                    sorted_window = sorted(window_totals)
                    quartile = sorted_window[-len(sorted_window)//4]
                    threshold = min_ratio * max(quartile, 0.5)
                    if centre_count >= threshold:
                        yield self.get_consensus_for_length(centre_length)
                window_totals.pop(0)
