#!/usr/bin/env python3
"""
Calculate exact coverage for every base in contigs.

Exact coverage is defined as how many reads map around that base exactly,
meaning there are no mutations, deletions, or insertions when aligned.

This tool uses k-mer hashing for fast exact matching of reads to contigs.
"""

import argparse
import csv
import sys
import logging
from collections import defaultdict
from gzip import GzipFile
from io import TextIOWrapper
from pathlib import Path
from typing import Dict, Sequence, Tuple, TextIO, Iterator, cast
from Bio import SeqIO

logger = logging.getLogger(__name__)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="""Calculate exact coverage for every base in contigs.""",
        epilog="""Exact coverage is defined as how many reads map around that base exactly,
        meaning there are no mutations, deletions, or insertions when aligned.
        When the whole contig is covered exactly by reads, we can be confident
        that the contig represents at least one genome from the sample.
        """,
    )

    parser.add_argument(
        "fastq1",
        help="<input> FASTQ file containing forward reads (read 1), can be gzip compressed",
    )
    parser.add_argument(
        "fastq2",
        help="<input> FASTQ file containing reverse reads (read 2), can be gzip compressed",
    )
    parser.add_argument(
        "contigs_file",
        type=argparse.FileType("r"),
        help="<input> FASTA or CSV file containing contig sequences (CSV: 'sequence' or 'contig' column)",
    )
    parser.add_argument(
        "output_csv",
        type=argparse.FileType("w"),
        help="<output> CSV file with exact coverage counts per base",
    )

    verbosity_group = parser.add_mutually_exclusive_group()
    verbosity_group.add_argument(
        "--verbose", action="store_true", help="Increase output verbosity."
    )
    verbosity_group.add_argument(
        "--no-verbose",
        action="store_true",
        help="Normal output verbosity.",
        default=True,
    )
    verbosity_group.add_argument(
        "--debug", action="store_true", help="Maximum output verbosity."
    )
    verbosity_group.add_argument(
        "--quiet", action="store_true", help="Minimize output verbosity."
    )

    return parser


def open_fastq(filename: Path) -> TextIO:
    """
    Open a FASTQ file for reading, automatically handling gzip compression.

    Detects compression by file extension (.gz). Opens compressed files
    using GzipFile and wraps in TextIOWrapper for text mode reading.

    :param filename: Path to FASTQ file (can be .fastq or .fastq.gz)
    :return: Text file handle for reading
    """
    is_gzipped = str(filename).endswith(".gz")

    if is_gzipped:
        # Open in binary mode, wrap with GzipFile, then TextIOWrapper for text mode
        binary_file = open(filename, "rb")
        gzip_file = GzipFile(fileobj=binary_file)
        return TextIOWrapper(gzip_file)
    else:
        # Regular text file
        return open(filename, "r")


def read_fastq_pairs(fastq1: TextIO, fastq2: TextIO) -> Iterator[Tuple[str, str]]:
    """
    Read paired FASTQ files and yield (read1_seq, read2_seq) tuples.

    :param fastq1: Forward reads FASTQ file
    :param fastq2: Reverse reads FASTQ file
    :yield: Tuples of (forward_sequence, reverse_sequence)
    """
    # Read 4 lines at a time from each file (FASTQ format)
    while True:
        # Read forward read
        header1 = fastq1.readline()
        if not header1:
            break
        seq1 = fastq1.readline().strip()
        fastq1.readline()  # plus line
        fastq1.readline()  # quality line

        # Read reverse read
        header2 = fastq2.readline()
        if not header2:
            break
        seq2 = fastq2.readline().strip()
        fastq2.readline()  # plus line
        fastq2.readline()  # quality line

        yield (seq1, seq2)


def read_contigs(contigs_file: TextIO) -> Dict[str, str]:
    """
    Read contigs from either FASTA or CSV format.

    Automatically detects format based on file content:
    - FASTA: starts with '>'
    - CSV: has header line, looks for 'sequence' or 'contig' column for sequences

    For CSV files:
    - Sequence column: prioritizes 'sequence' over 'contig'
    - Name column: uses 'sample' or 'region' (in that order), falls back to position

    :param contigs_file: File handle to read contigs from
    :return: Dictionary mapping contig_name -> sequence
    """
    contigs = {}

    # Peek at first line to detect format
    first_line = contigs_file.readline().strip()
    contigs_file.seek(0)  # Reset to beginning

    if first_line.startswith(">"):
        # FASTA format
        logger.debug("Detected FASTA format for contigs file")
        for record in SeqIO.parse(contigs_file, "fasta"):
            contigs[record.id] = str(record.seq).upper()
    else:
        # CSV format
        logger.debug("Detected CSV format for contigs file")
        reader = csv.DictReader(contigs_file)

        for i, row in enumerate(reader, start=1):
            # Find sequence column: prioritize 'sequence' over 'contig'
            seq_column = None
            if "sequence" in row:
                seq_column = "sequence"
            elif "contig" in row:
                seq_column = "contig"
            else:
                raise ValueError(
                    f"CSV must have either 'sequence' or 'contig' column. Found columns: {list(row.keys())}"
                )

            contig_seq = row[seq_column]
            if not contig_seq:
                continue  # Skip empty sequences

            # Find name column: prioritize 'sample', then 'region'
            # Fall back to position if neither is present
            contig_name = None
            if "sample" in row and row["sample"]:
                contig_name = row["sample"]
            elif "region" in row and row["region"]:
                contig_name = row["region"]
            else:
                contig_name = f"contig{i}"

            # If same name appears multiple times, append index
            if contig_name in contigs:
                contig_name = f"{contig_name}_{i}"

            contigs[contig_name] = contig_seq.upper()

    logger.debug(f"Loaded {len(contigs)} contigs")
    return contigs


def reverse_complement(seq: str) -> str:
    """
    Return the reverse complement of a DNA sequence.

    :param seq: DNA sequence
    :return: Reverse complement
    """
    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(complement.get(base, base) for base in reversed(seq))


def build_kmer_index(
    contigs: Dict[str, str], kmer_size: int
) -> Dict[str, Sequence[Tuple[str, int]]]:
    """
    Build a k-mer index for all contigs at the specified k-mer size.

    :param contigs: Dictionary mapping contig_name -> sequence
    :param kmer_size: Size of k-mers to index
    :return: Dictionary mapping kmer -> list of (contig_name, position) tuples
    """
    kmer_index = defaultdict(list)

    for contig_name, sequence in contigs.items():
        seq_len = len(sequence)
        for i in range(seq_len - kmer_size + 1):
            kmer = sequence[i : i + kmer_size]
            if "N" not in kmer:  # Skip k-mers with N
                kmer_index[kmer].append((contig_name, i))

    ret = cast(Dict[str, Sequence[Tuple[str, int]]], kmer_index)
    return ret


def build_kmer_index_for_size(
    contigs: Dict[str, str], k: int
) -> Dict[str, Sequence[Tuple[str, int]]]:
    """
    Build a k-mer index for a specific k-mer size.
    Used for lazy computation when encountering reads shorter than the default kmer_size.

    :param contigs: Dictionary mapping contig_name -> sequence
    :param k: Size of k-mers to index
    :return: Dictionary mapping kmer -> list of (contig_name, position) tuples
    """
    kmer_index = defaultdict(list)

    for contig_name, sequence in contigs.items():
        seq_len = len(sequence)
        for i in range(seq_len - k + 1):
            kmer = sequence[i : i + k]
            if "N" not in kmer:  # Skip k-mers with N
                kmer_index[kmer].append((contig_name, i))

    ret = cast(Dict[str, Sequence[Tuple[str, int]]], kmer_index)
    return ret


def find_exact_matches(
    read_seq: str,
    kmer_index: Dict[int, Dict[str, Sequence[Tuple[str, int]]]],
    contigs: Dict[str, str],
) -> Iterator[Tuple[str, int, int]]:
    """
    Find exact matches of a read in contigs using k-mer hashing.

    Lazily builds k-mer indices as needed when encountering reads of different lengths.

    :param read_seq: Read sequence
    :param kmer_index: Multi-level k-mer index mapping k-mer size -> (kmer -> [(contig_name, position)])
    :param contigs: Dictionary of contig sequences
    :return: Iterator of (contig_name, start_pos, end_pos) tuples for exact matches
    """

    read_len = len(read_seq)

    if read_len == 0:
        return

    # Check if we already have an index for this read length
    if read_len not in kmer_index:
        # Build it lazily
        logger.debug(f"Building k-mer index for size {read_len}")
        kmer_index[read_len] = build_kmer_index_for_size(contigs, read_len)

    # Use the appropriate index
    effective_index = kmer_index[read_len]

    # Use entire read as search key
    search_kmer = read_seq

    for contig_name, contig_pos in effective_index[search_kmer]:
        yield (contig_name, contig_pos, contig_pos + read_len)


def calculate_exact_coverage(
    fastq1_filename: Path,
    fastq2_filename: Path,
    contigs_file: TextIO,
) -> Tuple[Dict[str, Sequence[int]], Dict[str, str]]:
    """
    Calculate exact coverage for every base in contigs.

    :param fastq1_filename: Path to forward reads FASTQ file (can be gzipped)
    :param fastq2_filename: Path to reverse reads FASTQ file (can be gzipped)
    :param contigs_file: FASTA or CSV file with contigs
    :return: Tuple of (coverage_dict, contigs_dict) where coverage_dict maps
             contig_name -> list of coverage counts and contigs_dict maps
             contig_name -> sequence
    """
    # Read contigs
    logger.debug("Reading contigs...")
    contigs = read_contigs(contigs_file)

    logger.debug(f"Loaded {len(contigs)} contigs")

    # Initialize coverage arrays
    coverage = {}
    for contig_name, sequence in contigs.items():
        coverage[contig_name] = [0] * len(sequence)
        logger.debug(f"Initialized coverage for {contig_name} ({len(sequence)} bases)")

    # Initialize k-mer index structure (multi-level: k-mer size -> index)
    kmer_index: Dict[int, Dict[str, Sequence[Tuple[str, int]]]] = {}

    # Process read pairs - open files with automatic gzip detection
    logger.debug("Processing reads...")
    read_count = 0
    match_count = 0

    with open_fastq(fastq1_filename) as fastq1, open_fastq(fastq2_filename) as fastq2:
        for read1_seq, read2_seq in read_fastq_pairs(fastq1, fastq2):
            read_count += 1
            if read_count % 100000 == 0:
                logger.debug(
                    f"Processed {read_count} read pairs, {match_count} exact matches found"
                )

            # Try forward orientation for read1
            for read_seq in [read1_seq, read2_seq]:
                # Try both forward and reverse complement
                for seq in [read_seq, reverse_complement(read_seq)]:
                    matches = find_exact_matches(seq, kmer_index, contigs)

                    for contig_name, start_pos, end_pos in matches:
                        match_count += 1
                        # Increment coverage for all bases covered by this read
                        for i in range(start_pos, end_pos):
                            coverage[contig_name][i] += 1

    logger.debug(f"Finished processing {read_count} read pairs")
    logger.debug(f"Total exact matches: {match_count}")

    # Report on lazy k-mer indices built
    if kmer_index:
        read_sizes = sorted(kmer_index.keys())
        logger.debug(
            f"Built {len(kmer_index)} k-mer indices for read sizes: {read_sizes}"
        )
    else:
        logger.debug("No k-mer indices built (no reads processed)")

    coverage_ret = cast(Dict[str, Sequence[int]], coverage)
    return coverage_ret, contigs


def write_coverage_csv(
    coverage: Dict[str, Sequence[int]], contigs: Dict[str, str], output_csv: TextIO
) -> None:
    """
    Write coverage data to CSV file.

    :param coverage: Dictionary mapping contig_name -> coverage array
    :param contigs: Dictionary mapping contig_name -> sequence
    :param output_csv: Output CSV file
    """
    writer = csv.DictWriter(
        output_csv,
        ["contig", "position", "exact_coverage"],
        lineterminator="\n",
    )
    writer.writeheader()

    for contig_name in sorted(coverage.keys()):
        contig_coverage = coverage[contig_name]

        for pos, cov in enumerate(contig_coverage, start=1):
            writer.writerow(
                {
                    "contig": contig_name,
                    "position": pos,
                    "exact_coverage": cov,
                }
            )


def main_typed(
    fastq1_filename: Path,
    fastq2_filename: Path,
    contigs_file: TextIO,
    output_csv: TextIO,
) -> None:
    coverage, contigs = calculate_exact_coverage(
        fastq1_filename, fastq2_filename, contigs_file
    )
    write_coverage_csv(coverage, contigs, output_csv)


def main(argv: Sequence[str]) -> int:
    parser = get_parser()
    args = parser.parse_args(argv)

    # Configure logging based on verbosity flags
    if args.quiet:
        logger.setLevel(logging.ERROR)
    elif args.verbose:
        logger.setLevel(logging.INFO)
    elif args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.WARN)

    logging.basicConfig(
        level=logger.level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    main_typed(
        Path(args.fastq1),
        Path(args.fastq2),
        args.contigs_file,
        args.output_csv,
    )

    logger.debug("Exact coverage calculation complete!")
    return 0


def entry() -> None:
    sys.exit(main(sys.argv[1:]))


if __name__ == "__main__": entry()  # noqa
