
from pathlib import Path
import json
import csv

from micall.utils.dir_path import DirPath
from micall.utils.new_atomic_file import new_atomic_text_file
from .logger import logger


FIELDNAMES = ("app",
              "concordance",
              "depth",
              "alignment_score",
              "mlen",
              "total_mlen",
              "overlap_count",
              "number_of_contigs",
              "avg_contigs_size",
              "run_time",
              "run_id",
              "sample",
              )


def combine_runs_stats(root: DirPath, runs_txt: Path, target: Path) -> None:
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

            selected = {key: value for key, value in stats_object.items()
                        if key in FIELDNAMES}
            dwriter.writerow(selected)
