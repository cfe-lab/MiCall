
import json
from typing import Mapping
from pathlib import Path
import subprocess

from micall.utils.dir_path import DirPath
from micall.utils.new_atomic_directory import new_atomic_directory
from .logger import logger


FILEFILTER = '.*((sample_info)|(coverage_score)|(genome_co)|(contigs)).*'


def process_info(root: DirPath, info: Mapping[str, object]) -> None:
    id = info["id"]
    output = root / "runs" / str(id)
    info_path = output / "info.json"

    if info_path.exists():
        logger.debug("Directory for RUN_ID %s already exists. Skipping...", id)
        return

    with new_atomic_directory(output) as output:
        subprocess.check_call(["kivecli",
                               "download",
                               "--debug",
                               "--run_id", str(id),
                               "--output", str(output),
                               "--filefilter", FILEFILTER,
                               ])

        with info_path.open("w") as writer:
            json.dump(info, writer, indent='\t')


def download(root: DirPath, json_file: Path) -> None:
    with json_file.open() as reader:
        for info in json.load(reader):
            process_info(root, info)
