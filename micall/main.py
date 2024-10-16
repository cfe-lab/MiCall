#! /usr/bin/env python

import sys
import argparse
import os
from typing import Sequence
from pathlib import Path
from importlib.metadata import version


EXECUTABLES = [
    "release_test_publish.py",
    "micall_kive.py",
    "micall_watcher.py",
    "release_test_microtest.py",
    "docker_build.py",
    "release_test_setup.py",
    "release_test_compare.py",
    "micall_kive_resistance.py",
    "micall_docker.py",
    "micall/main.py",
    "micall/resistance/genreport.py",
    "micall/resistance/resistance.py",
    "micall/resistance/pdfreport.py",
    "micall/core/filter_quality.py",
    "micall/core/sam2aln.py",
    "micall/core/denovo.py",
    "micall/core/trim_fastqs.py",
    "micall/core/plot_contigs.py",
    "micall/core/cascade_report.py",
    "micall/core/remap.py",
    "micall/core/prelim_map.py",
    "micall/core/aln2counts.py",
    "micall/core/contig_stitcher.py",
    "micall/core/coverage_plots.py",
    "micall/core/plot_simple.py",
    "micall/tests/test_installation.py",
    "micall/tests/test_hcv_rules_import.py",
    "micall/g2p/fastq_g2p.py",
    "micall/blast_db/make_blast_db.py",
    "micall/utils/concordance_evaluation.py",
    "micall/utils/basespace_upload.py",
    "micall/utils/compare_mapping.py",
    "micall/utils/project_seeds_from_compendium.py",
    "micall/utils/fasta_to_csv.py",
    "micall/utils/hcv_rules_import.py",
    "micall/utils/dd.py",
    "micall/utils/find_reads_in_sam.py",
    "micall/utils/hcv_rules_display.py",
    "micall/utils/coverage_data.py",
    "micall/utils/find_by_coverage.py",
    "micall/utils/primer_locations.py",
    "micall/utils/fetch_sequences.py",
    "micall/utils/sam_g2p_simplify.py",
    "micall/utils/contig_summary.py",
    "micall/utils/compare_454_samples.py",
    "micall/utils/genreport_rerun.py",
    "micall/utils/remove_dupe_dirs.py",
    "micall/utils/find_missing_samples.py",
    "micall/utils/denovo_simplify.py",
    "micall/utils/sort_sam.py",
    "micall/utils/sample_fastq.py",
    "micall/utils/sample_sheet_parser.py",
    "micall/utils/projects_upload.py",
    "micall/utils/projects_dump.py",
    "micall/utils/find_chimera.py",
    "micall/utils/probe_finder.py",
    "micall/utils/aln2counts_simplify.py",
    "micall/utils/samples_from_454.py",
    "micall/utils/amplicon_finder.py",
    "micall/utils/driver_utils.py",
    "micall/utils/seed_alignments.py",
    "micall/utils/remap_fastq_simplify.py",
    "micall/utils/contig_counts.py",
    "micall/utils/ref_aligner.py",
    "micall/utils/scan_run_folders.py",
    "micall/utils/contig_blaster.py",
    "micall/utils/hcv_reference_tree.py",
    "micall/utils/sample_project_summary.py",
    "micall/utils/get_list_of_executables.py",
    "micall/monitor/update_qai.py",
    "micall/tcr/igblast.py",
]


def executable_name(path: str) -> str:
    file_name = Path(path).name
    name, extension = os.path.splitext(file_name)
    return name


def executable_module(path: str) -> str:
    noext, extension = os.path.splitext(path)
    pythized = noext.replace(os.path.sep, '.')
    return pythized


EXECUTABLES_MAP = {executable_name(path): path for path in EXECUTABLES}


def get_version() -> str:
    if __package__ is None:
        return "development"
    else:
        return str(version(__package__))


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run MiCall script.", add_help=False)
    parser.add_argument("--version", action="store_true", help="Print version and exit.")
    parser.add_argument('--help', action='store_true', help='Show this help message and exit.')
    parser.add_argument("program", nargs='?', choices=EXECUTABLES_MAP.keys(), help="Program name.")
    parser.add_argument("arguments", nargs=argparse.REMAINDER, help="Program arguments.")
    return parser


def main(argv: Sequence[str]) -> int:
    parser = get_parser()
    args = parser.parse_args(argv)

    if args.version:
        print(get_version())
        return 0

    elif args.help:
        parser.print_help()
        return 0

    elif EXECUTABLES_MAP.get(args.program):
        path = EXECUTABLES_MAP[args.program]
        mod = executable_module(path)
        evalstring = f'__import__({mod!r}, fromlist=[""])'
        evaluated_module = eval(evalstring)
        return evaluated_module.main(args.arguments)

    else:
        parser.print_help()
        return 1


def cli() -> int:
    return main(sys.argv[1:])


if __name__ == '__main__':
    exit(cli())
