
from pathlib import Path
from typing import Iterator, Iterable
import subprocess
import json

from micall.utils.dir_path import DirPath
from micall.utils.user_error import UserError
from .batch import BatchName
from .generate_setup_stage_ninjafile import generate_setup_stage_ninjafile


def read_batches(batches_list: Path) -> Iterator[BatchName]:
    with batches_list.open() as reader:
        for line in reader:
            yield BatchName(line.strip())


def combine_run_jsons(outputs: Iterable[Path], target: Path) -> None:
    result = []

    for output in outputs:
        with output.open() as reader:
            batches_json = json.load(reader)
            for batch_json in batches_json:
                runs = batch_json["runs"]
                result.extend(runs)

    with target.open("w") as writer:
        json.dump(result, writer, indent='\t')


def run_all(batches_list: Path, root: DirPath, properties: Path) -> None:
    root.mkdir(exist_ok=True, parents=True)
    batches = tuple(read_batches(batches_list))

    setup_stage_ninjafile = root / "setup.ninja"
    outputs = tuple(generate_setup_stage_ninjafile(root, batches, target=setup_stage_ninjafile))

    try:
        subprocess.check_call(["ninja", "-f", str(setup_stage_ninjafile)])
    except BaseException as ex:
        raise UserError("Work failed: %s", str(ex)) from ex

    runs_json = root / "runs.json"
    combine_run_jsons(outputs, runs_json)
