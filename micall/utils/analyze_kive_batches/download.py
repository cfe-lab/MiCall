
import json
from typing import Mapping
from pathlib import Path
import subprocess

from micall.utils.dir_path import DirPath
from .logger import logger


FILEFILTER = '.*((sample_info)|(coverage_score)|(genome_co)|(contigs)).*'


def process_info(root: DirPath, info: Mapping[str, object]) -> None:
    id = info["id"]
    output = root / "runs" / str(id)
    if output.exists():
        logger.debug("Directory for RUN_ID %s already exists. Skipping...", id)
        return

    subprocess.check_call(["kivecli",
                           "download",
                           "--debug",
                           "--run_id", str(id),
                           "--output", str(output),
                           "--filefilter", FILEFILTER,
                           ])

    output.mkdir(exist_ok=True, parents=True)
    with (output / "info.json").open("w") as writer:
        json.dump(info, writer, indent='\t')


def download(root: DirPath, json_file: Path) -> None:
    with json_file.open() as reader:
        for info in json.load(reader):
            process_info(root, info)
