
from pathlib import Path
from typing import Iterator
import subprocess

from micall.utils.dir_path import DirPath
from micall.utils.user_error import UserError
from .batch import BatchName
from .generate_setup_stage_ninjafile import generate_setup_stage_ninjafile


def read_batches(batches_list: Path) -> Iterator[BatchName]:
    with batches_list.open() as reader:
        for line in reader:
            yield BatchName(line.strip())


def run_all(batches_list: Path, root: DirPath, properties: Path) -> None:
    root.mkdir(exist_ok=True, parents=True)
    batches = tuple(read_batches(batches_list))

    setup_stage_ninjafile = root / "setup.ninja"
    runs_json = root / "runs.json"
    generate_setup_stage_ninjafile(root,
                                   batches,
                                   target=setup_stage_ninjafile,
                                   runs_json=runs_json,
                                   )

    try:
        subprocess.check_call(["ninja", "-f", str(setup_stage_ninjafile)])
    except BaseException as ex:
        raise UserError("Work failed: %s", str(ex)) from ex
