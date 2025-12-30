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
import os
from Bio.Align import PairwiseAligner
from fractions import Fraction

from micall.utils.new_atomic_file import new_atomic_text_file
from micall.utils.dir_path import DirPath
from micall.core.project_config import ProjectConfig
from micall.utils.exact_coverage import calculate_exact_coverage
from micall.utils.list_fastq_files import find_fastq_source_folder
from .logger import logger
from .find_file import find_file


# set up the aligner
ALIGNER = PairwiseAligner()
ALIGNER.mode = "global"
# By default end gaps are free (score 0).
# Make them penalized by picking -1 per gap position at ends.
ALIGNER.end_gap_score = -1


PROJECTS = ProjectConfig.loadDefault()


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
        comb = Fraction(math.comb(L, x))
        pval += (comb *
                 Fraction((match_prob ** x) *
                          ((1 - match_prob) ** (L - x))))

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


def read_conseqs_rows(path_to_file: Path) -> Rows:
    with open(path_to_file) as fd:
        reader = csv.DictReader(fd)
        for row in reader:
            # This field might be called "region" by mistake.
            ref = row.get("seed") or row.get("region")
            if ref is None:
                continue

            if "unknown" in ref.lower():
                continue

            cutoff = row["consensus-percent-cutoff"]
            if cutoff != "MAX":
                continue

            yield row


def read_contigs_or_conseqs_rows(path_to_file: Path) -> Rows:
    with open(path_to_file) as fd:
        reader = csv.DictReader(fd)
        fieldnames = reader.fieldnames
        assert fieldnames is not None, f"Unexpected file {path_to_file} without columns."

    if "consensus-percent-cutoff" in fieldnames:
        yield from read_conseqs_rows(path_to_file)
    else:
        yield from read_contigs_rows(path_to_file)


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


def count_soft_clips_from_nuc(path_to_file: Path) -> Optional[int]:
    """
    Count the total number of soft-clipped reads from the nuc CSV file.

    The nuc.csv file contains per-position nucleotide counts including a 'clip' column.
    Each row represents a position in a region/gene, and the 'clip' value is the count
    of soft-clipped reads at that position.

    To avoid overcounting when the same reference position appears in multiple overlapping regions,
    we track the maximum clip count per unique (seed, refseq.nuc.pos) combination.

    @param path_to_file: Path to the nuc.csv file
    @return: Total count of soft-clipped reads, or None if file is empty/invalid
    """
    # Track max clip count per (seed, refseq.nuc.pos) to avoid overcounting
    # when positions appear in multiple overlapping regions
    position_clips: dict[tuple[str, str], int] = {}
    has_data = False

    try:
        with open(path_to_file) as fd:
            reader = csv.DictReader(fd)
            for row in reader:
                has_data = True
                clip_str = row.get("clip", "0")

                # Skip empty values
                if not clip_str:
                    continue

                try:
                    clip_count = int(clip_str)
                except ValueError:
                    logger.warning("Invalid clip value %r in nuc file at position %s.",
                                 clip_str, row.get("refseq.nuc.pos", "unknown"))
                    continue

                # Skip zero counts
                if clip_count == 0:
                    continue

                # Use (seed, refseq.nuc.pos) as key to deduplicate across overlapping regions
                seed = row.get("seed", "")
                refseq_pos = row.get("refseq.nuc.pos", "")

                # Skip if no reference position (shouldn't happen in valid data)
                if not refseq_pos:
                    continue

                key = (seed, refseq_pos)

                # Keep the maximum clip count for this position
                if key not in position_clips or clip_count > position_clips[key]:
                    position_clips[key] = clip_count

    except (IOError, OSError) as e:
        logger.warning("Could not read nuc file %s: %s", path_to_file, e)
        return None

    if not has_data:
        return None

    # Sum up all the deduplicated clip counts
    total_clips = sum(position_clips.values())
    return total_clips


def calculate_max_ex_cov_dip(
        sample: str,
        run_name: str,
        conseq_stitched_csv_path: Path,
        overlap_size: int = 70) -> Optional[int]:
    """
    Calculate maximum decrease in exact coverage between adjacent positions.

    :param sample: Sample name (e.g., "2130A")
    :param run_name: Run name from sample_info CSV
    :param conseq_stitched_csv_path: Path to conseq_stitched CSV file (final consensus)
    :param overlap_size: Minimum overlap size for exact coverage calculation
    :return: Maximum drop in exact coverage between adjacent positions, or None if files not found
    """
    # Get RAW_DATA path from environment variable
    raw_data_path = os.environ.get('RAW_DATA')
    if not raw_data_path:
        raise ValueError("Environment variable $RAW_DATA not set")

    raw_data_dir = Path(raw_data_path)
    if not raw_data_dir.exists():
        raise ValueError("RAW_DATA path does not exist: {}".format(raw_data_path))

    # Find run directory - try MiSeq/runs structure first
    run_dirs = list((raw_data_dir / 'MiSeq' / 'runs').glob(run_name + '*'))
    if not run_dirs:
        raise ValueError("Run directory not found for run name: {}".format(run_name))
    if len(run_dirs) > 1:
        raise ValueError("Multiple run directories found for run name: {}".format(run_name))

    run_dir = run_dirs[0]
    if not run_dir.exists():
        raise ValueError("Run directory does not exist: {}".format(run_dir))

    # Use find_fastq_source_folder to locate FASTQ files (handles various folder structures)
    fastq_folder = find_fastq_source_folder(run_dir)
    if fastq_folder is None:
        logger.warning("No FASTQ folder found for run %s", run_name)
        return None

    # Find R1 and R2 FASTQ files for this sample
    r1_pattern = f"{sample}_*_R1_*.fastq*"
    r2_pattern = f"{sample}_*_R2_*.fastq*"

    r1_files = list(fastq_folder.glob(r1_pattern))
    r2_files = list(fastq_folder.glob(r2_pattern))

    if not r1_files or not r2_files:
        logger.warning("FASTQ files not found for sample %s in %s (pattern: %s)",
                      sample, fastq_folder, r1_pattern)
        return None

    if len(r1_files) > 1 or len(r2_files) > 1:
        logger.warning("Multiple FASTQ files found for sample %s in %s, using first match",
                      sample, fastq_folder)

    fastq1_path = r1_files[0]
    fastq2_path = r2_files[0]

    # Calculate exact coverage using the exact_coverage utility with conseq_stitched CSV
    try:
        with open(conseq_stitched_csv_path, 'r') as conseq_file:
            coverage_dict, _contigs = calculate_exact_coverage(
                fastq1_path,
                fastq2_path,
                conseq_file,
                overlap_size
            )
    except Exception as ex:
        logger.error("Failed to calculate exact coverage for sample %s: %s", sample, ex)
        return None

    # Calculate maximum drop between adjacent positions across all consensus sequences
    max_drop = 0
    for contig_name, coverage_array in coverage_dict.items():
        for i in range(len(coverage_array) - 1):
            drop = abs(int(coverage_array[i]) - int(coverage_array[i + 1]))
            if drop > max_drop:
                max_drop = drop

    return max_drop if max_drop > 0 else 0


def calculate_alignment_scores(run_id: object, rows: Rows) -> Optional[float]:
    def collect_scores() -> Iterator[float]:
        for row in rows:
            # This field might be called "region" by mistake.
            ref_name_raw = row.get("seed") or row.get("region")
            if not ref_name_raw:
                raise ValueError("Cannot get name")

            ref_name_raw = str(ref_name_raw)

            try:
                #
                # Some ref_name_raws are prefix with "1-" or "2-" or "n-".
                # Strip that here.
                #
                index, _dash, ref_name = ref_name_raw.partition('-')
                int(index, 10)
            except BaseException:
                ref_name = ref_name_raw

            if ref_name.endswith('-partial'):
                ref_name = ref_name[:-len('-partial')]

            try:
                # ref_name = "HIV1-B-ZA-KP109515-seed"
                reference = PROJECTS.getReference(ref_name)
                assert reference is not None
            except BaseException:
                logger.warning("Invalid reference name %r in run %s.", ref_name, run_id)
                continue

            query = str(row["sequence"])
            # do the alignment, grab top hit
            aln = next(iter(ALIGNER.align(reference, query)))
            yield aln.score

    scores = collect_scores()
    return max(scores, default=None)


def get_stats(info_file: Path) -> Optional[Row]:
    with info_file.open() as reader:
        obj = json.load(reader)

    run_id = obj["id"]
    logger.debug("Processing %r.", run_id)
    o: Row = {}

    #
    # Copying from `$DIR.json`.
    #
    state = obj['state']
    app = obj['app_name']

    o['app'] = app
    o["run_id"] = run_id
    o["state"] = state
    if state != 'C':
        logger.error("Should only be processing runs with state 'C', but got state %s.", state)
        return o

    start_time = obj['start_time']
    end_time = obj['end_time']
    run_time = calculate_seconds_between(start_time, end_time)
    o['run_time'] = run_time

    assert info_file.name.endswith(".json")
    directory = DirPath(info_file.with_suffix(""))

    for subdir in directory.iterdir():
        if subdir.name.endswith("_info.csv"):
            sample = subdir.name[:-len("_info.csv")]
            sample_info_path = subdir
            break
    else:
        logger.warning("Cannot determine sample name for run %r.", run_id)
        sample = None
        sample_info_path = None
    o["sample"] = sample

    # Read run_name from sample_info CSV if available
    run_name = None
    if sample_info_path and sample_info_path.exists():
        try:
            with open(sample_info_path, 'r') as info_csv:
                info_reader = csv.DictReader(info_csv)
                info_row = next(info_reader, None)
                if info_row:
                    run_name = info_row.get('run_name')
        except Exception as ex:
            logger.warning("Failed to read run_name from %s: %s", sample_info_path, ex)

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

    try:
        conseqs_csv_path = find_file(directory, "^conseq.*[.]csv$")
    except ValueError as ex:
        conseqs_csv_path = None
        logger.error("%s", ex)

    try:
        conseqs_stitched_csv_path = find_file(directory, "^conseq_stitched.*[.]csv$")
    except ValueError as ex:
        conseqs_stitched_csv_path = None
        logger.error("%s", ex)

    try:
        nuc_csv_path = find_file(directory, "^nuc.*[.]csv$")
    except ValueError as ex:
        nuc_csv_path = None
        logger.error("%s", ex)

    rc = read_coverage_rows
    ro = read_contigs_or_conseqs_rows

    if the_csv_path:
        o["concordance"] = calculate_avg("concordance", rc(the_csv_path))
        o["depth"] = calculate_avg("coverage", rc(the_csv_path))
        o["total_mlen"] = count_matches(rc(the_csv_path))
        o["mlen"] = count_distinct_matches(rc(the_csv_path))
        o["overlap_count"] = count_duplicate_matches(rc(the_csv_path))

    if contigs_csv_path:
        o["number_of_contigs"] = count_contigs(ro(contigs_csv_path))
        o["avg_contigs_size"] = avg_contigs_lengths(ro(contigs_csv_path))

    if conseqs_csv_path:
        o["alignment_score"] = calculate_alignment_scores(run_id, ro(conseqs_csv_path))

    if conseqs_stitched_csv_path:
        o["stitched_alignment_score"] = calculate_alignment_scores(run_id, ro(conseqs_stitched_csv_path))

    if nuc_csv_path:
        o["soft_clips_count"] = count_soft_clips_from_nuc(nuc_csv_path)

    # Calculate maximum exact coverage dip based on final consensus (conseq_stitched)
    if sample and run_name and conseqs_stitched_csv_path:
        o["maximum_ex_cov_dip"] = calculate_max_ex_cov_dip(
            sample=sample,
            run_name=run_name,
            conseq_stitched_csv_path=conseqs_stitched_csv_path
        )
    else:
        if not sample:
            logger.warning("Cannot calculate maximum_ex_cov_dip: sample name not found")
        if not run_name:
            logger.warning("Cannot calculate maximum_ex_cov_dip: run_name not found in sample_info")
        if not conseqs_stitched_csv_path:
            logger.warning("Cannot calculate maximum_ex_cov_dip: conseq_stitched CSV not found")

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
