from collections import Counter
from typing import Optional

from micall.utils.translation import translate, ambig_dict

AMINO_ALPHABET = 'ACDEFGHIKLMNPQRSTVWY*'
MAX_CUTOFF = 'MAX'
FIRST_CUTOFF = 'FIRST'


class SeedAmino(object):
    """
    Records the frequencies of amino acids at a given position of the
    aligned reads as determined by the consensus sequence.
    """
    def __init__(self, consensus_nuc_index, counts=None):
        self.v3_overlap = 0
        self.consensus_nuc_index = consensus_nuc_index
        self.all_consensus_nuc_indexes = set()
        if consensus_nuc_index is not None:
            self.all_consensus_nuc_indexes.add(consensus_nuc_index)
        self.counts = counts or Counter()  # {amino: count}
        self.codon_counts: Counter = Counter()  # {codon_nucs: count}
        self.nucleotides = []
        for i in range(3):
            seed_nuc = SeedNucleotide()
            if consensus_nuc_index is not None:
                seed_nuc.consensus_index = consensus_nuc_index + i
            self.nucleotides.append(seed_nuc)
        self.low_quality = 0
        self.partial = 0
        self.deletions = 0
        self.read_count = 0
        self.ref_offset = 0
        self.nucleotides_to_skip = 0

    def __repr__(self):
        if self.counts:
            return 'SeedAmino({!r}, {!r})'.format(self.consensus_nuc_index,
                                                  dict(self.counts))
        return 'SeedAmino({})'.format(self.consensus_nuc_index)

    def count_aminos(self, codon_seq, count):
        """ Record a set of reads at this position in the seed reference.
        @param codon_seq: a string of three nucleotides that were read at this
                          position, may be padded with spaces at the start
                          or end of a sequence, or dashes for deletions
        @param count: the number of times they were read
        """
        self.read_count += count
        self.codon_counts[codon_seq] += count
        if 'N' in codon_seq:
            self.low_quality += count
        elif '---' == codon_seq:
            self.deletions += count
        elif '-' in codon_seq:
            self.partial += count  # Partial deletion
        elif ' ' not in codon_seq and 'n' not in codon_seq and len(codon_seq) == 3:
            amino = translate(codon_seq.upper())
            self.counts[amino] += count
        elif 'nnn' == codon_seq:
            # Don't count the gap between forward and reverse reads in a pair.
            self.read_count -= count
        for i, nuc in enumerate(codon_seq):
            if nuc != ' ':
                seed_nucleotide = self.nucleotides[i]
                seed_nucleotide.count_nucleotides(nuc, count)

    def add(self, other: 'SeedAmino', start_nuc: int = 0, end_nuc: int = 2):
        """ Add counts from another SeedAmino to this one.

        :param other: source to copy from
        :param start_nuc: first nucleotide index to copy: 0, 1, or 2.
        :param end_nuc: last nucleotide index to copy: 0, 1, or 2.
        """
        self.all_consensus_nuc_indexes.update(other.all_consensus_nuc_indexes)
        if (self.read_count and other.read_count and
                self.consensus_nuc_index != other.consensus_nuc_index):
            self.consensus_nuc_index = None
        elif other.read_count:
            self.consensus_nuc_index = other.consensus_nuc_index
        if 0 < start_nuc or end_nuc < 2:
            prefix = ' ' * start_nuc
            for nucs, count in other.codon_counts.items():
                self.count_aminos(prefix + nucs[start_nuc:end_nuc+1], count)
            if self.consensus_nuc_index == other.consensus_nuc_index:
                for seed_nuc, other_nuc in zip(self.nucleotides,
                                               other.nucleotides):
                    seed_nuc.consensus_index = other_nuc.consensus_index
        else:
            self.counts += other.counts
            for nuc, other_nuc in zip(self.nucleotides, other.nucleotides):
                nuc.add(other_nuc)
        self.partial += other.partial
        self.deletions += other.deletions
        self.read_count += other.read_count
        self.low_quality += other.low_quality
        self.nucleotides_to_skip = other.nucleotides_to_skip
        self.ref_offset = other.ref_offset

    def get_report(self) -> str:
        """ Build a report string with the counts of each amino acid.

        Report how many times each amino acid was seen in count_aminos().
        @return: comma-separated list of counts in the same order as the
        AMINO_ALPHABET list
        """
        return ','.join([str(self.counts[amino])
                         for amino in AMINO_ALPHABET])

    def apply_repeat(self, repeated_nuc: int) -> 'SeedAmino':
        new_amino = SeedAmino(self.consensus_nuc_index)
        for codon, count in self.codon_counts.items():
            new_codon = codon[:repeated_nuc + 1] + codon[repeated_nuc:2]
            new_amino.count_aminos(new_codon, count)
        return new_amino

    def get_consensus(self) -> str:
        """ Find the amino acid that was seen most often in count_aminos().

        If there is a tie, just pick one of the tied amino acids.
        @return: the letter of the most common amino acid
        """
        consensus = self.counts.most_common(1)
        if consensus:
            return consensus[0][0]
        if self.read_count:
            return '?'
        return '-'

    def count_overlap(self, other):
        for nuc1, nuc2 in zip(self.nucleotides, other.nucleotides):
            nuc1.count_overlap(nuc2)
            self.v3_overlap = max(self.v3_overlap, nuc1.v3_overlap)


class SeedNucleotide(object):
    """
    Records the frequencies of nucleotides at a given position of the
    aligned reads as determined by the consensus sequence.
    """
    COUNTED_NUCS = 'ACTG-'

    def __init__(self, counts=None):
        self.v3_overlap = self.clip_count = self.insertion_count = 0
        self.counts = counts or Counter()
        self.consensus_index = None

    def __repr__(self):
        return 'SeedNucleotide({!r})'.format(dict(self.counts))

    def count_nucleotides(self, nuc_seq, count=1):
        """ Record a set of reads at this position in the seed reference.
        @param nuc_seq: a single nucleotide letter that was read at this
        position
        @param count: the number of times it was read
        """
        if nuc_seq == 'n':
            "Represents gap between forward and reverse read, ignore."
        else:
            self.counts[nuc_seq] += count

    def add(self, other):
        total_count = sum(self.counts.values())
        other_count = sum(other.counts.values())
        if total_count and other_count:
            # Only set to None if consensus indexes are different
            if self.consensus_index != other.consensus_index:
                self.consensus_index = None
        elif other_count:
            self.consensus_index = other.consensus_index
        self.counts += other.counts
        self.clip_count += other.clip_count
        self.insertion_count += other.insertion_count

    def get_report(self):
        """ Build a report string with the counts of each nucleotide.

        Report how many times each nucleotide was seen in count_nucleotides().
        @return: comma-separated list of counts for A, C, G, and T.
        """
        return ','.join(map(str, [self.counts[nuc] for nuc in 'ACGT']))

    def get_coverage(self):
        return sum(self.counts[nuc] for nuc in self.COUNTED_NUCS)

    def get_consensus(self, mixture_cutoff, no_coverage='', discard_deletions=False):
        """ Choose consensus nucleotide or mixture from the counts.

        @param mixture_cutoff: the minimum fraction of reads
            that a nucleotide must be found in for it to be considered,
            or MAX_CUTOFF to consider only the most common nucleotide.
        @param no_coverage: what to return when there are no reads mapped to
            this position.
        @param discard_deletions: whether to return nothing for the deletions (if False: '-')
        @return: The letter for the consensus nucleotide or mixture.
            Nucleotide mixtures are encoded by IUPAC symbols, and the most common
            nucleotide can be a mixture if there is a tie.
        """
        if not self.counts:
            return no_coverage

        coverage = self.get_coverage()
        if mixture_cutoff not in (MAX_CUTOFF, FIRST_CUTOFF):
            min_count = coverage * mixture_cutoff
        else:
            min_count = 0
        mixture = []
        for nuc, count in self.counts.most_common():
            if count < min_count:
                break
            if nuc in self.COUNTED_NUCS:
                mixture.append(nuc)
                if mixture_cutoff in (MAX_CUTOFF, FIRST_CUTOFF):
                    # Catch any ties before breaking out.
                    min_count = count

        has_deletion = '-' in mixture
        if has_deletion:
            mixture.remove('-')
        if len(mixture) > 1:
            mixture.sort()
            if mixture_cutoff == FIRST_CUTOFF:
                consensus = mixture[0]
            else:
                consensus = ambig_dict[''.join(mixture)]
        elif len(mixture) == 1:
            # no ambiguity
            consensus = mixture[0]
        else:
            # Nothing left to go in the mixture.
            consensus = '-' if has_deletion else 'N'
        if has_deletion:
            consensus = consensus.lower()
        if consensus == '-' and discard_deletions:
            consensus = ''
        return consensus

    def count_overlap(self, other):
        for nuc in 'ACGT':
            self.v3_overlap += other.counts[nuc]


class ReportNucleotide:
    def __init__(self, position: int, seed_nucleotide: Optional[SeedNucleotide] = None):
        self.position = position
        if seed_nucleotide is None:
            self.seed_nucleotide = SeedNucleotide()
        else:
            self.seed_nucleotide = seed_nucleotide

    def __repr__(self):
        return f'ReportNucleotide({self.position!r}, {self.seed_nucleotide!r})'


class ReportAmino(object):
    def __init__(self, seed_amino: SeedAmino, position: int):
        """ Create a new instance.

        @param seed_amino: Counts for the
        """
        self.seed_amino = seed_amino
        self.position = position
        self.max_clip_count = 0
        self.insertion_count = 0

    def __repr__(self):
        return 'ReportAmino({!r}, {})'.format(self.seed_amino, self.position)
