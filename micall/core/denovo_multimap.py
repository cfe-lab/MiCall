"""Select equally best paired placements for de novo contig remapping.

Only the de novo pipeline uses this module.  A placement is one complete,
concordant R1/R2 Bowtie2 alignment; a secondary placement is not an extra
molecule.  The downstream max-depth consensus intentionally includes each
compatible *placement* as an alternative hypothesis.
"""
from collections import Counter
from csv import DictReader
from itertools import groupby
import re

# Full Bowtie2 reporting is used for small panels or modest reference
# collections. Limit the contig count AND total reference bases for the
# expanded tier: output can grow dramatically for repetitive inputs.
MAX_DENOVO_ALIGNMENTS = 8
ALL_ALIGNMENT_MAX_CONTIGS = 16
ALL_ALIGNMENT_EXPANDED_MAX_CONTIGS = 32
ALL_ALIGNMENT_EXPANDED_MAX_BASES = 65536
DENOVO_BOWTIE_THREADS = 4


def require_complete_denovo_mapping_for_clinical_report(counts_file):
    """Refuse clinical reporting when alternatives are known incomplete.

    At a finite -k reporting limit, additional equally good placements
    may be absent. Oversized QNAME groups are similarly discarded. Such
    unknown depth may not become an unqualified clinical consensus.
    Capping does not prove omissions, and passing does not guarantee
    exhaustive Bowtie2 heuristic search.
    """
    reader = DictReader(counts_file)
    if not reader.fieldnames or not {'type', 'count'} <= set(reader.fieldnames):
        raise ValueError('Missing de novo remap QC type/count headers')
    failures: Counter[str] = Counter()
    raw_present = False
    for row in reader:
        kind = row['type']
        if kind == 'raw':
            raw_present = True
        if kind in ('multimap-capped-pairs', 'multimap-oversized-pairs'):
            n = int(row['count'])
            if n < 0:
                raise ValueError('Negative de novo remap QC count')
            failures[kind] += n
    if not raw_present:
        raise ValueError('De novo remap QC has no raw input count')
    if any(failures.values()):
        details = ', '.join(f'{kind}={n}' for kind, n in failures.items() if n)
        raise RuntimeError(
            'Clinical de novo consensus withheld: incomplete alignment '
            'search or discarded large read-pair groups (' + details
            + '). Do not interpret partial max depth as genotype evidence.')


def denovo_alignment_limit(contigs):
    """Return None for all reported hits, else a bounded search limit.

    Preserve the existing full-search path for <=16 contigs of any length.
    Also use it on moderate 17–32-contig panels with a small overall
    reference footprint. The expanded tier has a deliberate resource
    ceiling, not a claim of exhaustive Bowtie2 heuristic search.
    """
    count = len(contigs)
    if count <= ALL_ALIGNMENT_MAX_CONTIGS:
        return None
    if (count <= ALL_ALIGNMENT_EXPANDED_MAX_CONTIGS
            and sum(map(len, contigs.values())) <=
            ALL_ALIGNMENT_EXPANDED_MAX_BASES):
        return None
    return MAX_DENOVO_ALIGNMENTS
# Never buffer more than this many SAM records for a single read pair.
# Pathologically repetitive reads are excluded, not partially miscounted.
MAX_DENOVO_RECORDS_PER_PAIR = 8192
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
    return None


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


def _valid_orphan(first, second):
    """Preserve an unambiguously mapped read when its mate was unmapped.

    This is a legacy MiCall behavior distinct from a second mapping
    hypothesis. Only the primary pair may provide singleton evidence.
    """
    if len(first) < 11 or len(second) < 11 or first[0] != second[0]:
        return False
    try:
        a, b = int(first[1]), int(second[1])
        mapped = [row for row in (first, second)
                  if not int(row[1]) & 0x4]
        unmapped = [row for row in (first, second)
                    if int(row[1]) & 0x4]
    except ValueError:
        return False
    if (len(mapped) != 1 or len(unmapped) != 1
            or (a | b) & (0x400 | 0x800 | 0x100)
            or not (a & b & 0x1)
            or not (a & 0x40) or not (b & 0x80)
            or (a & 0x80) or (b & 0x40)):
        return False
    good, missing = mapped[0], unmapped[0]
    try:
        flag = int(good[1])
        mate_flag = int(missing[1])
        start = int(good[3])
    except ValueError:
        return False
    return (
        bool(flag & 0x8) and not (flag & 0x2)
        and not (mate_flag & 0x8 or mate_flag & 0x2)
        and good[2] != '*' and start > 0
        and bool(VALID_CIGAR.fullmatch(good[5]))
        and 'M' in good[5]
        and sum(int(n) for n, op in REFERENCE_CIGAR.findall(good[5])
                if op in 'MIS') == len(good[9])
        and good[9] != '*' and len(good[9]) == len(good[10])
        and missing[5] == '*' and missing[9] != '*'
        and len(missing[9]) == len(missing[10])
        and missing[2] in ('*',good[2]))


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
    candidates = [p for p in pairs if _valid_pair(*p)]
    # Bowtie2 normally emits an AS tag on both mates. If scores exist,
    # require both mates to have parseable scores; otherwise an omitted
    # score (implicitly treated as zero) could defeat a real negative score.
    # Older MiCall test fixtures contain no AS tags, so allow wholly
    # unscored fixtures while still rejecting malformed mixed-score pairs.
    scored = [p for p in candidates
              if all(_score(row) is not None for row in p)]
    if scored:
        candidates = scored
    elif any(any(field.startswith('AS:i:') for field in row[11:])
             for pair in candidates for row in pair):
        candidates = []
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
        if len(primary) == 2:
            # Bowtie2 can report singleton mates in either record order.
            # Do not discard genuine single-mate coverage or pretend it is
            # an additional complete paired placement.
            ordered = sorted(primary,
                             key=lambda row: bool(int(row[1]) & 0x80))
            if _valid_orphan(*ordered):
                return ordered
            if _valid_mates(*ordered):
                a, b = int(ordered[0][1]), int(ordered[1][1])
                if (ordered[0][2] != ordered[1][2]
                        or not (a & 0x2 and b & 0x2)):
                    return ordered
        return []
    best_score = max((_score(a) or 0) + (_score(b) or 0)
                     for a, b in candidates)
    winners = []
    used = set()
    intervals: dict[str, list[tuple[tuple[int, int], tuple[int, int]]]] = {}
    candidates.sort(key=lambda pair: (
        pair[0][2], int(pair[0][3]), int(pair[1][3]),
        pair[0][5], pair[1][5], int(pair[0][1]) & 0x10,
        int(pair[1][1]) & 0x10, bool(int(pair[0][1]) & 0x100)))
    for a, b in candidates:
        if (_score(a) or 0) + (_score(b) or 0) != best_score:
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
        if limit is not None and len(winners) == limit:
            break
    return [row for pair in winners for row in pair]


def select_best_pairs(source, destination, limit=MAX_DENOVO_ALIGNMENTS,
                      diagnostics=None):
    """Write complete best-scoring alternatives, retain original QNAMEs.

    Input is Bowtie2's QNAME-grouped SAM stream. The output can be consumed
    by MiCall's existing two-at-a-time sam2aln matcher unchanged.
    """
    if limit is not None and limit < 1:
        raise ValueError("limit must be positive")
    counts: Counter[str] = Counter()

    def records():
        for line in source:
            if line.startswith("@"):
                destination.write(line)
            elif line.strip():
                yield line.rstrip("\n").split("\t")

    for _, group in groupby(records(), key=lambda row: row[0]):
        rows: list[list[str]] = []
        oversized = False
        for row in group:
            if len(rows) < MAX_DENOVO_RECORDS_PER_PAIR:
                rows.append(row)
            else:
                oversized = True
        if oversized:
            # Never emit a truncated, falsely exhaustive read-pair group.
            if diagnostics is not None:
                diagnostics['oversized_pair_groups'] = (
                    diagnostics.get('oversized_pair_groups', 0) + 1)
            continue
        if (diagnostics is not None and limit is not None
                and len(rows) // 2 >= limit):
            diagnostics['at_report_limit'] = (
                diagnostics.get('at_report_limit', 0) + 1)
        selected = _select(rows, limit)
        # Bowtie2 -a legitimately reports a different number of alternative
        # alignments for R1 and R2 when no concordant pair is found. Such a
        # group can have an odd number of rows yet still contain exactly one
        # safe primary discordant/singleton pair for legacy rescue. Do not
        # misclassify those legitimate outputs as malformed SAM streams.
        primary_rescue = (len(selected) == 2
                          and not _valid_pair(selected[0], selected[1]))
        if (diagnostics is not None and len(rows) % 2
                and not primary_rescue):
            diagnostics['incomplete_pair_groups'] = (
                diagnostics.get('incomplete_pair_groups', 0) + 1)
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
