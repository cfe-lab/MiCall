
from pathlib import Path
import subprocess

from .batch import BatchName


def get_batch_runs(batch: BatchName, target: Path) -> None:
    with target.open("wt") as writer:
        subprocess.check_call(
            ["kivecli", "findbatches", "--debug", "--json", "--name", str(batch)],
            stdout=writer,
        )
