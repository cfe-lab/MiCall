from typing import FrozenSet
from dataclasses import dataclass
from micall.utils.referenceless_contig_with_aligner import ContigWithAligner
from micall.utils.referenceless_score import Score, SCORE_NOTHING
from micall.utils.contig_stitcher_contigs import Contig


@dataclass(frozen=True)
class ContigsPath:
    # Contig representing all combined contigs in the path.
    whole: ContigWithAligner

    # Id's of contigs that comprise this path.
    parts_ids: FrozenSet[int]

    # Higher is better. This is an estimated probability that
    # all the components in this path came together by accident.
    probability: Score

    def score(self) -> Score:
        return self.probability

    def __lt__(self, other: 'ContigsPath') -> bool:
        return self.score() < other.score()

    def __le__(self, other: 'ContigsPath') -> bool:
        return self.score() <= other.score()

    def __gt__(self, other: 'ContigsPath') -> bool:
        return self.score() > other.score()

    def __ge__(self, other: 'ContigsPath') -> bool:
        return self.score() >= other.score()

    def has_contig(self, contig: Contig) -> bool:
        return contig.id in self.parts_ids

    @staticmethod
    def singleton(contig: ContigWithAligner) -> 'ContigsPath':
        return ContigsPath(whole=contig,
                           parts_ids=frozenset((contig.id,)),
                           probability=SCORE_NOTHING,
                           )
