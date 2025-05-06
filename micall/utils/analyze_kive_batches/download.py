
import json
from typing import Iterator, Iterable, Callable
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



def skip_if_failed(root: DirPath, run: KiveRun) -> Iterator[KiveRun]:
    output = root / "runs" / str(run.id)
    failed_path = output / "failed"

    if failed_path.exists():
        logger.warning("Skipping run %s - download failed last time.", run.id)
        return

    yield run


def try_fetch_info(root: DirPath, run: KiveRun) -> Iterator[KiveRun]:
    output = root / "runs" / str(run.id)
    failed_path = output / "failed"

    if run.is_finished:
        yield run
        return

    try:
        logger.debug("Fetching run info for %s - it has not finished last time.", run.id)
        run = kivecli.findrun.find_run(run_id=run.id.value)

    except BaseException as ex:
        logger.warning("Could not fetch run info %s: %s", run.id, ex)
        with new_atomic_directory(output) as output:
            failed_path = output / "failed"
            failed_path.touch()
            return

    if not run.is_finished:
        logger.warning("Run %s is still processing.", run.id)

    yield run


def skip_incomplete(root: DirPath, run: KiveRun) -> Iterator[KiveRun]:
    if run.state != RunState.COMPLETE:
        logger.warning("Skipping run %s: state '%s' != '%s'.",
                       run.id, run.state.value, RunState.COMPLETE.value)
        return

    yield run


def try_download(root: DirPath, run: KiveRun) -> Iterator[KiveRun]:
    output = root / "runs" / str(run.id)
    info_path = output / "info.json"

    if info_path.exists():
        logger.warning("Skipping run %s - downloaded last time.", run.id)
        return

    try:
        with new_atomic_directory(output) as output:
            kivecli.download.main_parsed(
                output=kivecli.dirpath.DirPath(output),
                run_id=run.id.value,
                nowait=False,
                filefilter=FILEFILTER,
            )

            info_path = output / "info.json"
            with info_path.open("w") as writer:
                json.dump(run, writer, indent='\t')

            yield run

    except BaseException as ex:
        logger.warning("Could not download run %s: %s", run.id, ex)
        with new_atomic_directory(output) as output:
            failed_path = output / "failed"
            failed_path.touch()


def pipeline(root: DirPath,
             runs: Iterable[KiveRun],
             *fns: Callable[[DirPath, KiveRun], Iterator[KiveRun]],
             ) -> Iterator[KiveRun]:

    runs = tuple(runs)
    for fn in fns:
        runs = tuple(result for run in runs for result in fn(root, run))

    yield from runs


def collect_run_ids(root: DirPath, runs: Iterable[KiveRun]) -> Iterator[KiveRun]:
    with login():
        yield from pipeline(root, runs,
                            skip_if_failed,
                            try_fetch_info,
                            skip_incomplete,
                            try_download,
                            )


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
