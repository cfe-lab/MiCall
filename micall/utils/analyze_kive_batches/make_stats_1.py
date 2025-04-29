import json
from typing import Iterator, \
    Iterable, Optional, Union, Mapping, MutableMapping
from pathlib import Path
import re
import csv
import math

from micall.utils.new_atomic_file import new_atomic_text_file
from micall.utils.dir_path import DirPath
from .logger import logger


Row = Mapping[str, Optional[Union[str, int, float]]]
Rows = Iterable[Row]


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


def calculate_avg(column: str, rows: Rows) -> str:
    num = average(collect_values(column, rows))
    if math.isnan(num):
        return ""
    else:
        return str(num)


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


def get_stats(directory: DirPath) -> Row:
    def find_file(directory: Path, pattern: str) -> Path:
        for subdir in directory.iterdir():
            name = subdir.name
            if re.findall(pattern, name):
                return directory / name
        raise ValueError(f"Cannot find file {pattern!r}"
                         f" in directory {str(directory)!r}.")

    logger.debug("Processing %r.", str(directory))

    o: MutableMapping[str, Optional[Union[int, str, float]]] = {}

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

    logger.debug("Processed %r.", str(directory))
    return o


def make_stats_1(input: DirPath, output: Path) -> None:
    result: Row = get_stats(input)
    with new_atomic_text_file(output) as stats_file:
        json.dump(result, stats_file, indent='\t')
