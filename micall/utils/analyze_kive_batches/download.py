
import json
from typing import Mapping, Iterator
from pathlib import Path
import kivecli.download
from kivecli.login import login
import kivecli.dirpath
import kivecli.runfilesfilter
import kivecli.logger
import kivecli.findrun
import logging

from micall.utils.dir_path import DirPath
from micall.utils.new_atomic_directory import new_atomic_directory
from micall.utils.new_atomic_file import new_atomic_text_file
from .logger import logger


FILEFILTER = kivecli.runfilesfilter.RunFilesFilter.parse(
    '.*((sample_info)|(coverage_score)|(genome_co)|(contigs)|(conseq)).*')


kivecli.logger.logger.setLevel(logging.DEBUG)


def process_info(root: DirPath, info: Mapping[str, object]) -> bool:
    run_id = int(str(info["id"]))
    output = root / "runs" / str(run_id)
    info_path = output / "info.json"
    failed_path = output / "failed"

    if failed_path.exists():
        logger.warning("Skipping RUN_ID %s - download failed last time.", run_id)
        return False

    if info["end_time"]:
        if info_path.exists():
            logger.debug("Directory for RUN_ID %s already exists.", run_id)
            return True

    else:
        if info_path.exists():
            with info_path.open() as reader:
                info = json.load(reader)

        if info["end_time"]:
            # Already processed this, no changes possible.
            return True
        else:
            info = kivecli.findrun.find_run(run_id=run_id)

        if not info["end_time"]:
            logger.warning("Run %s is still processing.", run_id)
            return False

    try:
        with new_atomic_directory(output) as output:
            kivecli.download.main_parsed(
                output=kivecli.dirpath.DirPath(output),
                run_id=run_id,
                nowait=False,
                filefilter=FILEFILTER,
            )

            info_path = output / "info.json"
            with info_path.open("w") as writer:
                json.dump(info, writer, indent='\t')
            return True

    except BaseException as ex:
        logger.warning("Could not download run %s: %s", run_id, ex)
        failed_path.touch()
        return False


def download(root: DirPath, runs_json: Path, runs_txt: Path) -> None:
    with runs_json.open() as reader:
        data = json.load(reader)

    def collect_run_ids() -> Iterator[str]:
        with login():
            for info in data:
                if process_info(root, info):
                    run_id = info["id"]
                    state = info["state"]
                    if state != "C":
                        logger.warning("Skipping run %s: state '%s' != 'C'.", run_id, state)
                        continue

                    yield str(run_id)

    new_content = '\n'.join(collect_run_ids())

    if runs_txt.exists():
        previous_content = runs_txt.read_text()
        if previous_content == new_content:
            return

    with new_atomic_text_file(runs_txt) as writer:
        writer.write(new_content)
