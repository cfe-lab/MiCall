
from pathlib import Path
import subprocess

from micall.utils.user_error import UserError
from .batch import BatchName


def get_batch(batch: BatchName, target: Path) -> None:
    with target.open("wt") as writer:
        try:
            subprocess.check_call(
                ["kivecli", "findbatches", "--debug", "--json", "--name", str(batch)],
                stdout=writer,
            )
        except BaseException as ex:
            raise UserError("Work failed: %s", str(ex)) from ex
