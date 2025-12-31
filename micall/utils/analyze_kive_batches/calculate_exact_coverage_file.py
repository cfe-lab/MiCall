from pathlib import Path
import os

from micall.utils.dir_path import DirPath
from micall.utils.exact_coverage import calculate_exact_coverage, NoContigsError
from micall.utils.list_fastq_files import find_fastq_source_folder
from .find_file import find_file
from .logger import logger
import csv


def calculate_exact_coverage_file(info_file: Path, output: Path) -> None:
    """Calculate exact coverage and save to CSV file.

    :param info_file: Path to the run info JSON file
    :param output: Path to output CSV file
    """
    assert info_file.name.endswith(".json")
    directory = DirPath(info_file.with_suffix(""))

    # Find sample info and conseq_stitched files
    sample_info_path = None
    for subdir in directory.iterdir():
        if subdir.name.endswith("_info.csv"):
            sample = subdir.name[: -len("_info.csv")]
            sample_info_path = subdir
            break
    else:
        logger.warning(
            "Cannot determine sample name for run %r, skipping exact coverage",
            info_file.name,
        )
        output.touch()
        return

    # Read run_name from sample_info CSV
    run_name = None
    if sample_info_path and sample_info_path.exists():
        try:
            with open(sample_info_path, "r") as info_csv:
                info_reader = csv.DictReader(info_csv)
                info_row = next(info_reader, None)
                if info_row:
                    run_name = info_row.get("run_name")
        except Exception as ex:
            logger.warning("Failed to read run_name from %s: %s", sample_info_path, ex)

    if not run_name:
        logger.warning(
            "No run_name found in sample_info for %r, skipping exact coverage",
            info_file.name,
        )
        output.touch()
        return

    # Find conseq_stitched CSV
    try:
        conseq_stitched_csv_path = find_file(directory, "^conseq_stitched.*[.]csv$")
    except ValueError:
        logger.debug(
            "No conseq_stitched file found for run %r, skipping exact coverage",
            info_file.name,
        )
        output.touch()
        return

    # Get RAW_DATA path
    raw_data_path = os.environ.get("RAW_DATA")
    if not raw_data_path:
        raise ValueError("Environment variable $RAW_DATA not set")

    raw_data_dir = Path(raw_data_path)
    if not raw_data_dir.exists():
        raise ValueError("RAW_DATA path does not exist: {}".format(raw_data_path))

    # Find run directory
    run_dirs = list((raw_data_dir / "MiSeq" / "runs").glob(run_name + "*"))
    if not run_dirs:
        raise ValueError("Run directory not found for run name: {}".format(run_name))
    if len(run_dirs) > 1:
        raise ValueError(
            "Multiple run directories found for run name: {}".format(run_name)
        )

    run_dir = run_dirs[0]
    if not run_dir.exists():
        raise ValueError("Run directory does not exist: {}".format(run_dir))

    # Find FASTQ folder
    fastq_folder = find_fastq_source_folder(run_dir)
    if fastq_folder is None:
        logger.warning(
            "No FASTQ folder found for run %s, skipping exact coverage", run_name
        )
        output.touch()
        return

    # Find R1 and R2 FASTQ files
    r1_pattern = f"{sample}_*_R1_*.fastq*"
    r2_pattern = f"{sample}_*_R2_*.fastq*"

    r1_files = list(fastq_folder.glob(r1_pattern))
    r2_files = list(fastq_folder.glob(r2_pattern))

    if not r1_files or not r2_files:
        logger.warning(
            "FASTQ files not found for sample %s in %s, skipping exact coverage",
            sample,
            fastq_folder,
        )
        output.touch()
        return

    if len(r1_files) > 1 or len(r2_files) > 1:
        logger.warning(
            "Multiple FASTQ files found for sample %s in %s, using first match",
            sample,
            fastq_folder,
        )

    fastq1_path = r1_files[0]
    fastq2_path = r2_files[0]

    # Calculate exact coverage
    try:
        with open(conseq_stitched_csv_path, "r") as conseq_file:
            coverage_dict, contigs = calculate_exact_coverage(
                fastq1_path, fastq2_path, conseq_file, overlap_size=70
            )

        # Write coverage to CSV
        with open(output, "w", newline="") as out_csv:
            writer = csv.DictWriter(
                out_csv, ["contig", "position", "exact_coverage"], lineterminator="\n"
            )
            writer.writeheader()

            for contig_name in sorted(coverage_dict.keys()):
                contig_coverage = coverage_dict[contig_name]
                for pos, cov in enumerate(contig_coverage, start=1):
                    writer.writerow(
                        {
                            "contig": contig_name,
                            "position": pos,
                            "exact_coverage": int(cov),
                        }
                    )

        logger.debug("Calculated exact coverage for run %r", info_file.name)

    except NoContigsError:
        logger.debug(
            "No contigs found in conseq_stitched file for run %r", info_file.name
        )
        output.touch()
    except Exception as ex:
        logger.error(
            "Failed to calculate exact coverage for run %r: %s", info_file.name, ex
        )
        output.touch()
