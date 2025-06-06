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
    # Does not contain contigs covered
    contigs_ids: FrozenSet[int]

    # Ids of contigs that are contained within this path.
    # Different from `contigs_ids` in that it also contains ones that are covered by this path.
    contains_contigs_ids: FrozenSet[int]

    # Higher is better. The lower the score is, the higher the probability
    # that all the components of `whole` came together on accident.
    score: Score

    def get_score(self) -> Score:
        return self.score

    def __lt__(self, other: 'ContigsPath') -> bool:
        return self.get_score() < other.get_score()

    def __le__(self, other: 'ContigsPath') -> bool:
        return self.get_score() <= other.get_score()

    def __gt__(self, other: 'ContigsPath') -> bool:
        return self.get_score() > other.get_score()

    def __ge__(self, other: 'ContigsPath') -> bool:
        return self.get_score() >= other.get_score()

    def has_contig(self, contig: Contig) -> bool:
        return contig.id in self.contains_contigs_ids

    @staticmethod
    def singleton(contig: ContigWithAligner) -> 'ContigsPath':
        return ContigsPath(whole=contig,
                           contigs_ids=frozenset((contig.id,)),
                           contains_contigs_ids=frozenset((contig.id,)),
                           score=SCORE_NOTHING,
                           )
