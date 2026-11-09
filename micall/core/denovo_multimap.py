"""Select equally best paired placements for de novo contig remapping.

Only the de novo pipeline uses this module.  A placement is one complete,
concordant R1/R2 Bowtie2 alignment; a secondary placement is not an extra
molecule.  The downstream max-depth consensus intentionally includes each
compatible *placement* as an alternative hypothesis.
"""
from collections import Counter
from itertools import groupby
import re

# Bound Bowtie2 output and downstream hypothesis fan-out per original pair.
# This is not an exhaustive search when a read has more than eight placements.
MAX_DENOVO_ALIGNMENTS = 8
BAD_FLAGS = 0x4 | 0x8 | 0x400 | 0x800
REFERENCE_CIGAR = re.compile(r"(\d+)([MIDS])")
VALID_CIGAR = re.compile(r"(?:[1-9][0-9]*[MIDS])+\Z")


def _score(row):
    for field in row[11:]:
        if field.startswith("AS:i:"):
            try:
                return int(field[5:])
            except ValueError:
                return None
    # Older MiCall test fixtures omit the optional Bowtie2 score.
    return 0


def _interval(row):
    start = int(row[3])
    consumed = sum(int(n) for n, op in REFERENCE_CIGAR.findall(row[5])
                   if op in "MD")
    return start, start + consumed - 1


def _valid_mates(first, second):
    """Only pair genuine R1 and R2 SAM records from the same molecule."""
    if len(first) < 11 or len(second) < 11:
        return False
    try:
        a, b = int(first[1]), int(second[1])
        left, right = int(first[3]), int(second[3])
    except ValueError:
        return False
    return (first[0] == second[0]
            and not ((a | b) & BAD_FLAGS)
            and (a & 0x41) == 0x41 and (b & 0x81) == 0x81
            and not (a & 0x80 or b & 0x40)
            and bool(a & 0x100) == bool(b & 0x100)
            and bool(a & 0x20) == bool(b & 0x10)
            and bool(b & 0x20) == bool(a & 0x10)
            and first[2] != "*" and second[2] != "*"
            and left > 0 and right > 0
            and all(VALID_CIGAR.fullmatch(row[5]) for row in (first, second))
            and all(row[9] != "*" and len(row[9]) == len(row[10])
                    for row in (first, second))
            and all(sum(int(n) for n, op in REFERENCE_CIGAR.findall(row[5])
                        if op in "MIS") == len(row[9])
                    for row in (first, second))
            and all("M" in row[5] for row in (first, second)))


def _valid_pair(first, second):
    if not _valid_mates(first, second):
        return False
    try:
        a, b = int(first[1]), int(second[1])
        mate_left, mate_right = int(first[7]), int(second[7])
    except ValueError:
        return False
    return (bool(a & 0x2) and bool(b & 0x2)
            and first[2] == second[2]
            and first[6] in ("=", first[2])
            and second[6] in ("=", second[2])
            and mate_left == int(second[3])
            and mate_right == int(first[3]))


def _overlaps(a, b):
    return a[0] <= b[1] and b[0] <= a[1]


def _select(rows, limit):
    # Bowtie2 -k emits the two mates of each reported paired placement
    # together. Fail closed for malformed ordering instead of cross-pairing.
    pairs = [(rows[i], rows[i + 1]) for i in range(0, len(rows) - 1, 2)]
    candidates = [p for p in pairs if _valid_pair(*p)
                  and _score(p[0]) is not None
                  and _score(p[1]) is not None]
    if not candidates:
        # Preserve MiCall's existing fallback for discordant *primary*
        # alignments. MixedReferenceSplitter can retry these against one
        # reference; properly paired but internally inconsistent rows must
        # never be silently resurrected by the fallback.
        primary = []
        for row in rows:
            try:
                if not int(row[1]) & 0x100:
                    primary.append(row)
            except (ValueError, IndexError):
                continue
        if len(primary) == 2 and _valid_mates(*primary):
            a, b = int(primary[0][1]), int(primary[1][1])
            if (primary[0][2] != primary[1][2]
                    or not (a & 0x2 and b & 0x2)):
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


def select_best_pairs(source, destination, limit=MAX_DENOVO_ALIGNMENTS,
                      diagnostics=None):
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
        rows = list(group)
        if diagnostics is not None and len(rows) // 2 >= limit:
            diagnostics['at_report_limit'] = (
                diagnostics.get('at_report_limit', 0) + 1)
        selected = _select(rows, limit)
        # Distinguish intentionally unmapped pairs from mapped records
        # rejected by SAM integrity or quality constraints.
        if (diagnostics is not None and not selected
                and any(len(row) > 2 and row[2] != '*' for row in rows)):
            diagnostics['discarded_mapped_pairs'] = (
                diagnostics.get('discarded_mapped_pairs', 0) + 1)
        for row in selected:
            destination.write("\t".join(row) + "\n")
            if not int(row[1]) & 0x4 and row[2] != "*":
                counts[row[2]] += 1
    return counts
