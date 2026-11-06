"""Select equally best paired placements for de novo contig remapping.

Only the de novo pipeline uses this module.  A placement is one complete,
concordant R1/R2 Bowtie2 alignment; a secondary placement is not an extra
molecule.  The downstream max-depth consensus intentionally includes each
compatible *placement* as an alternative hypothesis.
"""
from collections import Counter
from itertools import groupby
import re

BAD_FLAGS = 0x4 | 0x8 | 0x400 | 0x800
REFERENCE_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def _score(row):
    for field in row[11:]:
        if field.startswith("AS:i:"):
            return int(field[5:])
    return 0


def _interval(row):
    start = int(row[3])
    consumed = sum(int(n) for n, op in REFERENCE_CIGAR.findall(row[5])
                   if op in "MDN=X")
    return start, start + consumed - 1


def _valid_pair(first, second):
    a, b = int(first[1]), int(second[1])
    return (first[0] == second[0]
            and not ((a | b) & BAD_FLAGS)
            and (a & 0x43) == 0x43 and (b & 0x83) == 0x83
            and not (a & 0x80 or b & 0x40)
            and bool(a & 0x100) == bool(b & 0x100)
            and first[2] == second[2] != "*"
            and first[6] in ("=", first[2])
            and second[6] in ("=", second[2])
            and int(first[7]) == int(second[3])
            and int(second[7]) == int(first[3])
            and first[5] != "*" and second[5] != "*")


def _overlaps(a, b):
    return a[0] <= b[1] and b[0] <= a[1]


def _select(rows, limit):
    # Bowtie2 -k emits the two mates of each reported paired placement
    # together. Fail closed for malformed ordering instead of cross-pairing.
    pairs = [(rows[i], rows[i + 1]) for i in range(0, len(rows) - 1, 2)]
    candidates = [p for p in pairs if _valid_pair(*p)]
    if not candidates:
        # Preserve the original mixed-reference rescue path only when the
        # primary pair actually maps to two different references.
        primary = [row for row in rows if not int(row[1]) & 0x100]
        if (len(primary) == 2
                and primary[0][0] == primary[1][0]
                and primary[0][2] != primary[1][2]
                and all(not (int(row[1]) & BAD_FLAGS) for row in primary)):
            return primary
        return []
    best_score = max(_score(a) + _score(b) for a, b in candidates)
    winners = []
    used = set()
    intervals: dict[str, list[tuple[tuple[int, int], tuple[int, int]]]] = {}
    for a, b in candidates:
        if _score(a) + _score(b) != best_score:
            continue
        key = (a[2], a[3], a[5], b[3], b[5], int(a[1]) & 0x10,
               int(b[1]) & 0x10)
        if key in used:
            continue
        region = a[2]
        positions = (_interval(a), _interval(b))
        prev = intervals.setdefault(region, [])
        if any(any(_overlaps(x, y) for x in positions for y in prior)
               for prior in prev):
            # Never count one molecule twice at the same contig locus.
            continue
        used.add(key)
        prev.append(positions)
        winners.append((a, b))
        if len(winners) == limit:
            break
    return [row for pair in winners for row in pair]


def select_best_pairs(source, destination, limit=8):
    """Write complete best-scoring alternatives, retain original QNAMEs.

    Input is Bowtie2's QNAME-grouped SAM stream. The output can be consumed
    by MiCall's existing two-at-a-time sam2aln matcher unchanged.
    """
    if limit < 1:
        raise ValueError("limit must be positive")
    counts: Counter[str] = Counter()

    def records():
        for line in source:
            if line.startswith("@"):
                destination.write(line)
            elif line.strip():
                yield line.rstrip("\n").split("\t")

    for _, group in groupby(records(), key=lambda row: row[0]):
        for row in _select(list(group), limit):
            destination.write("\t".join(row) + "\n")
            if not int(row[1]) & 0x4 and row[2] != "*":
                counts[row[2]] += 1
    return counts
