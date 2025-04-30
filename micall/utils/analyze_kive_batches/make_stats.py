import json
from typing import Iterator, Sequence, \
    Iterable, Optional, Union, MutableMapping
from pathlib import Path
import re
import csv
import math
from datetime import datetime
from itertools import tee
import ast

from micall.utils.new_atomic_file import new_atomic_text_file
from micall.utils.dir_path import DirPath
from .logger import logger
from .find_file import find_file


Row = MutableMapping[str, Optional[Union[str, int, float, Sequence['Row']]]]
Rows = Iterable[Row]


def calc_overlap_pvalue(L: int, M: int, match_prob: float = 1/4) -> float:
    """
    Compute the probability (p-value) of observing at least M matches
    out of L under a binomial model where each position has
    probability `match_prob` of matching.

    NOTE: This is a very simplistic model.

    :param L: Total length of the overlap region (int)
    :param M: Number of matching positions observed (int)
    :param match_prob: Probability of a match at a single position (float),
                       default = 0.25 (uniform base model)
    :return: p-value (float), the probability of seeing at least M
    matches by chance
    """

    pval = 0.0

    # Accounting for the fact that the contigs differ at the ends
    # (left end is different and right end is different).
    L += 2

    # Summation of
    # Binomial(L, x) * match_prob^x * (1-match_prob)^(L-x)
    #   from x = M to L
    for x in range(M, L + 1):
        pval += (math.comb(L, x) *
                 (match_prob ** x) *
                 ((1 - match_prob) ** (L - x)))

    return pval


def get_differences(arr: Iterable[float]) -> Iterator[float]:
    arr, it = tee(arr)
    next(it)
    for x, y in zip(arr, it):
        yield y - x


def mult_arr(arr: Iterable[float]) -> Iterator[float]:
    arr, it = tee(arr)
    next(it)
    for x, y in zip(arr, it):
        yield y * x


def find_number_of_mismatches(arr: Sequence[float]) -> int:
    if len(arr) < 2:
        return 0

    rdiff = get_differences(arr)
    mul = mult_arr(rdiff)
    mismatches = sum(1 for x in mul if x < 0)
    return mismatches


expr = re.compile(r'and full concordance.*\]')
expr2 = re.compile(r'\[.*\]')


def find_overlaps(the_log_path: Path) -> Iterator[Sequence[float]]:
    with open(the_log_path) as reader:
        for line in reader:
            mat = expr.findall(line)
            if mat:
                assert len(mat) == 1
                mat = expr2.findall(mat[0])
                assert len(mat) == 1
                arr: Sequence[float] = ast.literal_eval(mat[0])
                yield arr


def calculate_seconds_between(start_time: str, end_time: str) -> float:
    # Parse the start and end times into datetime objects
    start_dt = datetime.fromisoformat(start_time)
    end_dt = datetime.fromisoformat(end_time)

    # Calculate the difference between the two datetime objects
    duration = end_dt - start_dt

    # Return the total seconds as a float
    return duration.total_seconds()


def read_coverage_rows(path_to_file: Path) -> Rows:
    with open(path_to_file) as fd:
        reader = csv.DictReader(fd)
        for row in reader:
            contig = row["contig"]
            if "unknown" in contig.lower():
                continue

            yield row


def read_contigs_rows(path_to_file: Path) -> Rows:
    with open(path_to_file) as fd:
        reader = csv.DictReader(fd)
        for row in reader:
            ref = row["ref"]
            if "unknown" in ref.lower():
                continue

            yield row


def average(values: Iterable[float]) -> float:
    total = 0.0
    count = 0

    for value in values:
        total += value
        count += 1

    return total / count if count > 0 else float('nan')


def collect_values(column: str, rows: Rows) -> Iterator[float]:
    for row in rows:
        try:
            value = row[column]
            assert isinstance(value, str)
            assert len(value) > 0
            numeric_value = float(value)
        except BaseException:
            continue
        yield numeric_value


def calculate_avg(column: str, rows: Rows) -> Optional[float]:
    num = average(collect_values(column, rows))
    if math.isnan(num):
        return None
    else:
        return num


def count_matches(rows: Rows) -> Optional[int]:
    total = 0
    anything = False

    for row in rows:
        anything = True
        value = row["link"]
        if value == "M":
            total += 1

    if anything:
        return total
    else:
        return None


def count_distinct_matches(rows: Rows) -> Optional[int]:
    ret = set()
    anything = False

    for row in rows:
        anything = True
        value = row["link"]
        if value == "M":
            ret.add(row["refseq_nuc_pos"])

    if anything:
        return len(ret)
    else:
        return None


def count_duplicate_matches(rows: Rows) -> int:
    counter = 0
    ret = set()
    for row in rows:
        value = row["link"]
        if value == "M":
            pos = row["refseq_nuc_pos"]
            if pos in ret:
                counter += 1
            ret.add(pos)
    return counter


def collect_contigs_lengths(rows: Rows) -> Iterator[int]:
    for row in rows:
        ref = row["ref"]
        assert isinstance(ref, str)
        if "unknown" in ref.lower():
            continue
        else:
            contig = row["contig"]
            assert isinstance(contig, str)
            yield len(contig)


def count_contigs(rows: Rows) -> int:
    return sum(1 for x in collect_contigs_lengths(rows))


def avg_contigs_lengths(rows: Rows) -> float:
    lengths = list(collect_contigs_lengths(rows))
    if lengths:
        return sum(lengths) / len(lengths)
    else:
        return 0


def get_stats(info_file: Path) -> Optional[Row]:
    with info_file.open() as reader:
        obj = json.load(reader)

    run_id = obj["id"]
    logger.debug("Processing %r.", run_id)
    o: Row = {}

    #
    # Copying from `info.json`.
    #
    state = obj['state']
    app_name = obj['app_name']
    safe_app = app_name.replace(':', '-') \
                       .replace('/', '-') \
                       .replace(' ', '-') \
                       .replace('--', '-') \
                       .replace('--', '-') \
                       .replace('--', '-') \
                       .replace('--', '-')

    o['app'] = safe_app
    o["run_id"] = run_id
    o["state"] = state

    if state in ['N', 'L', 'R']:
        logger.warning("Run %r is still going.", run_id)
        return None

    if state != 'C':
        logger.warning("Run %r is incomplete.", run_id)
        return o

    start_time = obj['start_time']
    end_time = obj['end_time']
    run_time = calculate_seconds_between(start_time, end_time)
    o['run_time'] = run_time

    directory = DirPath(info_file.parent)

    for subdir in directory.iterdir():
        if subdir.name.endswith("_info.csv"):
            sample = subdir.name[:-len("_info.csv")]
            break
    else:
        logger.warning("Cannot determine sample name for run %r.", run_id)
        sample = None
    o["sample"] = sample

    try:
        the_csv_path = find_file(directory, "genome_coverage.*[.]csv$")
    except ValueError as ex:
        the_csv_path = None
        logger.error("%s", ex)

    try:
        contigs_csv_path = find_file(directory, "^contigs_.*[.]csv$")
    except ValueError as ex:
        contigs_csv_path = None
        logger.error("%s", ex)

    try:
        the_log_path = find_file(directory, ".*stitcher.*[.]log$")
    except ValueError as ex:
        the_log_path = None
        logger.error("%s", ex)

    rc = read_coverage_rows
    ro = read_contigs_rows

    if the_csv_path:
        o["concordance"] = calculate_avg("concordance", rc(the_csv_path))
        o["depth"] = calculate_avg("coverage", rc(the_csv_path))
        o["total_mlen"] = count_matches(rc(the_csv_path))
        o["mlen"] = count_distinct_matches(rc(the_csv_path))
        o["overlap_count"] = count_duplicate_matches(rc(the_csv_path))

    if contigs_csv_path:
        o["number_of_contigs"] = count_contigs(ro(contigs_csv_path))
        o["avg_contigs_size"] = avg_contigs_lengths(ro(contigs_csv_path))

    overlaps: list[Row] = []
    o["overlaps"] = overlaps
    if the_log_path:
        for overlap in find_overlaps(the_log_path):
            overlap_size = len(overlap)
            overlap_mismatches = find_number_of_mismatches(overlap)
            overlap_pvalue = calc_overlap_pvalue(
                L=overlap_size, M=(overlap_size - overlap_mismatches))

            overlaps.append({
                "overlap_size": overlap_size,
                "overlap_mismatches": overlap_mismatches,
                "overlap_pvalue": overlap_pvalue,
            })

    logger.debug("Processed %r.", run_id)
    return o


def make_stats(input: Path, output: Path) -> None:
    result = get_stats(input)
    if result:
        with new_atomic_text_file(output) as stats_file:
            json.dump(result, stats_file, indent='\t')
