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
from typing import Dict, List, Sequence, Tuple, TextIO, Iterator
from Bio import SeqIO

logger = logging.getLogger(__name__)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Calculate exact coverage for every base in contigs using k-mer hashing."
    )

    parser.add_argument(
        "fastq1",
        type=argparse.FileType("r"),
        help="<input> FASTQ file containing forward reads (read 1)",
    )
    parser.add_argument(
        "fastq2",
        type=argparse.FileType("r"),
        help="<input> FASTQ file containing reverse reads (read 2)",
    )
    parser.add_argument(
        "contigs_fasta",
        type=argparse.FileType("r"),
        help="<input> FASTA file containing contig sequences",
    )
    parser.add_argument(
        "output_csv",
        type=argparse.FileType("w"),
        help="<output> CSV file with exact coverage counts per base",
    )
    parser.add_argument(
        "--kmer-size", type=int, default=31, help="K-mer size for hashing (default: 31)"
    )
    parser.add_argument(
        "--min-overlap",
        type=int,
        default=20,
        help="Minimum overlap required for exact match (default: 20)",
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
) -> Dict[str, List[Tuple[str, int]]]:
    """
    Build a k-mer index for all contigs.

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

    return kmer_index


def find_exact_matches(
    read_seq: str,
    kmer_index: Dict[str, List[Tuple[str, int]]],
    contigs: Dict[str, str],
    kmer_size: int,
    min_overlap: int,
) -> List[Tuple[str, int, int]]:
    """
    Find exact matches of a read in contigs using k-mer hashing.

    :param read_seq: Read sequence
    :param kmer_index: K-mer index of contigs
    :param contigs: Dictionary of contig sequences
    :param kmer_size: Size of k-mers
    :param min_overlap: Minimum overlap required
    :return: List of (contig_name, start_pos, end_pos) tuples for exact matches
    """
    if len(read_seq) < min_overlap:
        return []

    matches = []
    read_len = len(read_seq)

    # Try to find k-mer matches
    if read_len >= kmer_size:
        first_kmer = read_seq[:kmer_size]
        if first_kmer in kmer_index:
            # Check each potential match location
            for contig_name, contig_pos in kmer_index[first_kmer]:
                contig_seq = contigs[contig_name]
                # Check if the entire read matches exactly at this position
                if contig_pos + read_len <= len(contig_seq):
                    if contig_seq[contig_pos : contig_pos + read_len] == read_seq:
                        matches.append((contig_name, contig_pos, contig_pos + read_len))

    return matches


def calculate_exact_coverage(
    fastq1: TextIO,
    fastq2: TextIO,
    contigs_fasta: TextIO,
    kmer_size: int = 31,
    min_overlap: int = 20,
) -> Tuple[Dict[str, List[int]], Dict[str, str]]:
    """
    Calculate exact coverage for every base in contigs.

    :param fastq1: Forward reads FASTQ file
    :param fastq2: Reverse reads FASTQ file
    :param contigs_fasta: FASTA file with contigs
    :param kmer_size: Size of k-mers for hashing
    :param min_overlap: Minimum overlap required
    :return: Tuple of (coverage_dict, contigs_dict) where coverage_dict maps
             contig_name -> list of coverage counts and contigs_dict maps
             contig_name -> sequence
    """
    # Read contigs
    logger.info("Reading contigs...")
    contigs = {}
    for record in SeqIO.parse(contigs_fasta, "fasta"):
        contigs[record.id] = str(record.seq).upper()

    logger.info(f"Loaded {len(contigs)} contigs")

    # Initialize coverage arrays
    coverage = {}
    for contig_name, sequence in contigs.items():
        coverage[contig_name] = [0] * len(sequence)
        logger.debug(f"Initialized coverage for {contig_name} ({len(sequence)} bases)")

    # Build k-mer index
    logger.info("Building k-mer index...")
    kmer_index = build_kmer_index(contigs, kmer_size)
    logger.info(f"Indexed {len(kmer_index)} unique k-mers")

    # Process read pairs
    logger.info("Processing reads...")
    read_count = 0
    match_count = 0

    for read1_seq, read2_seq in read_fastq_pairs(fastq1, fastq2):
        read_count += 1
        if read_count % 100000 == 0:
            logger.info(
                f"Processed {read_count} read pairs, {match_count} exact matches found"
            )

        # Try forward orientation for read1
        for read_seq in [read1_seq, read2_seq]:
            # Try both forward and reverse complement
            for seq in [read_seq, reverse_complement(read_seq)]:
                matches = find_exact_matches(
                    seq, kmer_index, contigs, kmer_size, min_overlap
                )

                for contig_name, start_pos, end_pos in matches:
                    match_count += 1
                    logger.debug(
                        f"Match: {contig_name}:{start_pos}-{end_pos} (length {end_pos - start_pos})"
                    )
                    # Increment coverage for all bases covered by this read
                    for i in range(start_pos, end_pos):
                        coverage[contig_name][i] += 1

    logger.info(f"Finished processing {read_count} read pairs")
    logger.info(f"Total exact matches: {match_count}")

    return coverage, contigs


def write_coverage_csv(
    coverage: Dict[str, List[int]], contigs: Dict[str, str], output_csv: TextIO
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

    coverage, contigs = calculate_exact_coverage(
        args.fastq1, args.fastq2, args.contigs_fasta, args.kmer_size, args.min_overlap
    )

    write_coverage_csv(coverage, contigs, args.output_csv)

    logger.info("Exact coverage calculation complete!")
    return 0


def entry() -> None:
    sys.exit(main(sys.argv[1:]))


if __name__ == "__main__": entry()  # noqa
