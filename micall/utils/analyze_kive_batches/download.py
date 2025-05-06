
import json
from typing import Iterator, Iterable, Callable, Optional
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


def mark_run_as_failed(root: DirPath, run: KiveRun) -> None:
    output = root / "runs" / str(run.id)
    with new_atomic_directory(output) as tmpout:
        failed_path = tmpout / "failed"
        failed_path.touch()


def skip_if_failed(root: DirPath, run: KiveRun) -> Optional[KiveRun]:
    output = root / "runs" / str(run.id)
    failed_path = output / "failed"

    if failed_path.exists():
        logger.warning("Skipping run %s - download failed last time.", run.id)
        return None

    return run


def try_fetch_info(root: DirPath, run: KiveRun) -> Optional[KiveRun]:
    if run.is_finished:
        return run

    try:
        logger.debug("Fetching run info for %s - it has not finished last time.", run.id)
        run = kivecli.findrun.find_run(run_id=run.id.value)

    except BaseException as ex:
        logger.warning("Could not fetch run info %s: %s", run.id, ex)
        mark_run_as_failed(root, run)
        return None

    if not run.is_finished:
        logger.warning("Run %s is still processing.", run.id)
        return None

    return run


def skip_incomplete(root: DirPath, run: KiveRun) -> Optional[KiveRun]:
    if run.state != RunState.COMPLETE:
        logger.warning("Skipping run %s: state '%s' != '%s'.",
                       run.id, run.state.value, RunState.COMPLETE.value)
        return None

    return run


def download_run_files(root: DirPath, run: KiveRun) -> None:
    output = root / "runs" / str(run.id)
    with new_atomic_directory(output) as tmpout:
        kivecli.download.main_parsed(
            output=kivecli.dirpath.DirPath(tmpout),
            run_id=run.id.value,
            nowait=False,
            filefilter=FILEFILTER,
        )

        info_path = tmpout / "info.json"
        with info_path.open("w") as writer:
            run.dump(writer)


def try_download(root: DirPath, run: KiveRun) -> Optional[KiveRun]:
    output = root / "runs" / str(run.id)
    info_path = output / "info.json"

    if info_path.exists():
        logger.warning("Skipping run %s - downloaded last time.", run.id)
        return run

    try:
        download_run_files(root, run)
        return run
    except BaseException as ex:
        logger.warning("Could not download run %s: %s", run.id, ex)
        mark_run_as_failed(root, run)
        return None


def pipeline(root: DirPath,
             runs: Iterable[KiveRun],
             *fns: Callable[[DirPath, KiveRun], Optional[KiveRun]],
             ) -> Iterator[KiveRun]:

    runs = tuple(runs)
    for fn in fns:
        maybes = tuple(fn(root, run) for run in runs)
        runs = tuple(mayberun for mayberun in maybes if mayberun is not None)

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
