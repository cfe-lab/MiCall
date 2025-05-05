
import json
from typing import Iterator, Optional, Iterable
from pathlib import Path
import kivecli.download
from kivecli.login import login
import kivecli.dirpath
import kivecli.runfilesfilter
import kivecli.logger
import kivecli.findrun
from kivecli.kiverun import KiveRun
from kivecli.runstate import RunState
import logging

from micall.utils.dir_path import DirPath
from micall.utils.new_atomic_directory import new_atomic_directory
from micall.utils.new_atomic_file import new_atomic_text_file
from .logger import logger


FILEFILTER = kivecli.runfilesfilter.RunFilesFilter.parse(
    '.*((sample_info)|(coverage_score)|(genome_co)|(contigs)|(conseq)).*')


kivecli.logger.logger.setLevel(logging.DEBUG)


def process_info(root: DirPath, info: KiveRun) -> Optional[KiveRun]:
    run_id = info.id
    output = root / "runs" / str(run_id)
    info_path = output / "info.json"
    failed_path = output / "failed"

    if failed_path.exists():
        logger.warning("Skipping RUN_ID %s - download failed last time.", run_id)
        return None

    if info.is_finished:
        if info_path.exists():
            logger.debug("Directory for RUN_ID %s already exists.", run_id)
            if info.state != RunState.COMPLETE:
                logger.warning("Skipping run %s: state '%s' != 'C'.", run_id, info.state.value)
                return None

            return info

    else:
        if info_path.exists():
            with info_path.open() as reader:
                info = KiveRun.from_json(json.load(reader))

            if info.is_finished:
                logger.debug("Directory for RUN_ID %s already exists.", run_id)
                return info

        info = kivecli.findrun.find_run(run_id=run_id.value)
        if not info.is_finished:
            logger.warning("Run %s is still processing.", run_id)
            return None

    if info.state != RunState.COMPLETE:
        logger.warning("Skipping run %s: state '%s' != 'C'.", run_id, info.state.value)
        return None

    try:
        with new_atomic_directory(output) as output:
            kivecli.download.main_parsed(
                output=kivecli.dirpath.DirPath(output),
                run_id=run_id.value,
                nowait=False,
                filefilter=FILEFILTER,
            )

            info_path = output / "info.json"
            with info_path.open("w") as writer:
                json.dump(info, writer, indent='\t')
            return info

    except BaseException as ex:
        logger.warning("Could not download run %s: %s", run_id, ex)
        with new_atomic_directory(output) as output:
            failed_path = output / "failed"
            failed_path.touch()
        return None


def collect_run_ids(root: DirPath, runs: Iterable[KiveRun]) -> Iterator[KiveRun]:
    with login():
        for info in runs:
            ret = process_info(root, info)
            if ret is not None:
                yield ret


def download(root: DirPath, runs_json: Path, runs_txt: Path) -> None:
    with runs_json.open() as reader:
        runs_raw = json.load(reader)

    runs: Iterable[KiveRun] = tuple(KiveRun.from_json(run) for run in runs_raw)
    new_runs = tuple(collect_run_ids(root, runs))
    new_content = json.dumps(list(run.raw for run in new_runs), indent='\t')

    if runs_txt.exists():
        previous_content = runs_txt.read_text()
        if previous_content == new_content:
            return

    with new_atomic_text_file(runs_txt) as writer:
        writer.write(new_content)
