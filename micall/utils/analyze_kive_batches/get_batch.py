
from pathlib import Path
import json
import kivecli.findbatches
import kivecli.login

from micall.utils.user_error import UserError
from .batch import BatchName


def get_batch(batch: BatchName, target: Path) -> None:
    with target.open("wt") as writer:
        try:
            with kivecli.login.login():
                batches = list(kivecli.findbatches.findbatches(name=str(batch)))
                json.dump(batches, writer, indent='\t')
        except BaseException as ex:
            raise UserError("Work failed: %s", str(ex)) from ex
