
from pathlib import Path
import json
import csv

from micall.utils.dir_path import DirPath
from micall.utils.new_atomic_file import new_atomic_text_file
from .logger import logger


FIELDNAMES = ("app",
              "overlap_size",
              "overlap_mismatches",
              "overlap_pvalue",
              "run_id",
              "sample",
              )


def combine_runs_overlaps(root: DirPath, runs_txt: Path, target: Path) -> None:
    run_ids = runs_txt.read_text().splitlines()
    with new_atomic_text_file(target) as writer:
        dwriter = csv.DictWriter(writer, fieldnames=FIELDNAMES)
        dwriter.writeheader()
        for run_id in run_ids:
            stats = root / "runs" / str(run_id) / "stats.json"
            if not stats.exists():
                logger.debug("No stats file for run %s.", run_id)
                continue

            with open(stats) as stats_reader:
                stats_object = json.load(stats_reader)

            overlaps = stats_object.get("overlaps", [])
            for overlap in overlaps:
                overlap["app"] = stats_object["app"]
                overlap["run_id"] = stats_object["run_id"]
                overlap["sample"] = stats_object["sample"]
                dwriter.writerow(overlap)
