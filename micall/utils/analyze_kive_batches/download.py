
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


def is_downloaded(root: DirPath, run: KiveRun) -> bool:
    output = root / "runs" / str(run.id)
    download_path = output / "downloaded"
    return download_path.exists()


def skip_if_failed(root: DirPath, run: KiveRun) -> Optional[KiveRun]:
    output = root / "runs" / str(run.id)
    failed_path = output / "failed"

    if failed_path.exists():
        logger.warning("Skipping run %s - download failed last time.", run.id)
        return None

    return run


def save_run_info(root: DirPath, run: KiveRun) -> None:
    info_path = root / "runs" / f"{run.id}.json"
    with new_atomic_text_file(info_path) as writer:
        run.dump(writer)


def try_load_run_info(root: DirPath, run: KiveRun) -> KiveRun:
    info_path = root / "runs" / f"{run.id}.json"
    if info_path.exists():
        with info_path.open() as reader:
            data = json.load(reader)
            return KiveRun.from_json(data)
    else:
        return run


def try_fetch_info(root: DirPath, run: KiveRun) -> Optional[KiveRun]:
    info_path = root / "runs" / f"{run.id}.json"
    existing = info_path.exists()

    if not existing:
        save_run_info(root, run)

    if run.is_finished:
        return run

    if existing:
        run = try_load_run_info(root, run)
        if run.is_finished:
            return run

    try:
        logger.debug("Fetching run info for %s - it has not finished last time.", run.id)
        run = run.update()

    except BaseException as ex:
        logger.warning("Could not fetch run info %s: %s", run.id, ex)
        mark_run_as_failed(root, run)
        return None

    save_run_info(root, run)
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
        (tmpout / "downloaded").touch()


def try_download(root: DirPath, run: KiveRun) -> Optional[KiveRun]:
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

    yield from (run for run in runs if is_downloaded(root, run))
    runs = tuple(run for run in runs if not is_downloaded(root, run))

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


def download(root: DirPath, runs_json: Path, downloaded_runs_json: Path) -> None:
    with runs_json.open() as reader:
        runs_raw = json.load(reader)

    runs: Iterable[KiveRun] = tuple(KiveRun.from_json(run) for run in runs_raw)
    new_runs = tuple(collect_run_ids(root, runs))
    new_content = json.dumps(list(run.raw for run in new_runs), indent='\t')

    if downloaded_runs_json.exists():
        previous_content = downloaded_runs_json.read_text()
        if previous_content == new_content:
            return

    with new_atomic_text_file(downloaded_runs_json) as writer:
        writer.write(new_content)
