import json
from typing import Sequence, Iterator, \
    Iterable, Optional, Union, MutableMapping
from pathlib import Path
import re
from itertools import tee
import ast
import math
import subprocess

from .logger import logger
from micall.utils.new_atomic_file import new_atomic_text_file
from micall.utils.dir_path import DirPath


Row = MutableMapping[str, Optional[Union[str, int, float]]]
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


def get_stats(info_file: Path) -> Optional[Row]:
    def find_file(directory: DirPath, pattern: str) -> Path:
        for subdir in directory.iterdir():
            name = subdir.name
            if re.findall(pattern, name):
                return directory / name
        raise ValueError(f"Cannot find file {pattern!r}"
                         f" in directory {str(directory)!r}.")

    with info_file.open() as reader:
        obj = json.load(reader)

    run_id = obj["id"]
    logger.debug("Processing %r.", run_id)

    state = obj['state']
    if state != 'C':
        logger.warning("Run %r is incomplete.", run_id)
        return None

    directory = DirPath(info_file.parent)
    o: Row = {}

    try:
        the_unstitched_contigs_path = find_file(directory, ".*unstitched.*contig.*[.]csv$")
    except ValueError:
        try:
            the_unstitched_contigs_path = find_file(directory, ".*contigs.*[.]csv$")
        except ValueError as ex:
            logger.error("%s", ex)
            return None

    plot = directory / "stitcher_plot.svg"
    log = directory / "stitcher.log"
    if log.exists() and plot.exists():
        logger.debug("Run %r already stitched.", run_id)
    else:
        with open(log, "w") as log_writer:
            subprocess.check_call(
                ["micall", "contig_stitcher",
                 "--debug",
                 "--plot", str(plot),
                 str(the_unstitched_contigs_path),
                 "/dev/null",
                 ],
                stderr=log_writer,
            )

    try:
        the_log_path = find_file(directory, ".*stitcher.*[.]log$")
    except ValueError as ex:
        the_log_path = None
        logger.error("%s", ex)

    #
    # Copying from `info.json`.
    #
    app_name = obj['app_name']
    category = app_name.replace(':', '-') \
                       .replace('/', '-') \
                       .replace(' ', '-') \
                       .replace('--', '-') \
                       .replace('--', '-') \
                       .replace('--', '-') \
                       .replace('--', '-')

    o['assembler'] = category
    o["run_id"] = run_id

    if the_log_path:
        for overlap in find_overlaps(the_log_path):
            overlap_size = len(overlap)
            overlap_mismatches = find_number_of_mismatches(overlap)
            overlap_pvalue = calc_overlap_pvalue(
                L=overlap_size, M=(overlap_size - overlap_mismatches))
            o["overlap_size"] = overlap_size
            o["overlap_mismatches"] = overlap_mismatches
            o["overlap_pvalue"] = "{:.10f}".format(overlap_pvalue)
            return o

    logger.debug("Processed %r.", run_id)
    return o


def calculate_overlaps(input: Path, output: Path) -> None:
    result = get_stats(input)
    if result:
        with new_atomic_text_file(output) as stats_file:
            json.dump(result, stats_file, indent='\t')
