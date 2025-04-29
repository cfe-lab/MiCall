
from pathlib import Path
import json
import csv

from micall.utils.dir_path import DirPath
from micall.utils.new_atomic_file import new_atomic_text_file


FIELDNAMES = ("assembler",
              "concordance",
              "depth",
              "mlen",
              "total_mlen",
              "overlap_count",
              "number_of_contigs",
              "avg_contigs_size",
              "run_time",
              "runid",
              "sample",
              )


def combine_runs_stats(root: DirPath, runs_json: Path, target: Path) -> None:
    with runs_json.open() as reader:
        runs = json.load(reader)
        run_ids = tuple(x["id"] for x in runs)

    with new_atomic_text_file(target) as writer:
        dwriter = csv.DictWriter(writer, fieldnames=FIELDNAMES)
        dwriter.writeheader()
        for run_id in run_ids:
            stats = root / "runs" / str(run_id) / "stats.json"
            with open(stats) as stats_reader:
                stats_object = json.load(stats_reader)

            selected = {key: value for key, value in stats_object.items()
                        if key in FIELDNAMES}
            dwriter.writerow(selected)
