
import json
from typing import Mapping
from pathlib import Path
import subprocess

from micall.utils.dir_path import DirPath
from micall.utils.new_atomic_directory import new_atomic_directory
from .logger import logger


FILEFILTER = '.*((sample_info)|(coverage_score)|(genome_co)|(contigs)|(conseq)).*'


def process_info(root: DirPath, info: Mapping[str, object]) -> None:
    id = info["id"]
    output = root / "runs" / str(id)
    info_path = output / "info.json"

    if info_path.exists():
        logger.debug("Directory for RUN_ID %s already exists.", id)

        with info_path.open() as reader:
            existing_info = json.load(reader)

        if existing_info["state"] == info["state"]:
            logger.debug("Run %s has no updates.", id)
            return
        else:
            logger.info("Run %s has new updates.", id)

    try:
        with new_atomic_directory(output) as output:
            subprocess.check_call(["kivecli",
                                   "download",
                                   "--debug",
                                   "--run_id", str(id),
                                   "--output", str(output),
                                   "--filefilter", FILEFILTER,
                                   ])

            info_path = output / "info.json"
            with info_path.open("w") as writer:
                json.dump(info, writer, indent='\t')

    except BaseException as ex:
        logger.warning("Could not download run %s: %s", id, ex)


def download(root: DirPath, runs_json: Path, runs_txt: Path) -> None:
    run_ids = frozenset(runs_txt.read_text().splitlines())
    with runs_json.open() as reader:
        for info in json.load(reader):
            run_id = str(info["id"])
            if run_id in run_ids:
                process_info(root, info)
