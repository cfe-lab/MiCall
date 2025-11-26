
from pathlib import Path
import json
import kivecli.findbatches
import kivecli.login

from micall.utils.user_error import UserError
from micall.utils.new_atomic_file import new_atomic_text_file
from .batch import BatchName


def get_batch(batch: BatchName, target: Path) -> None:
    with new_atomic_text_file(target) as writer:
        try:
            with kivecli.login.login():
                batches = list(kivecli.findbatches.findbatches(name=str(batch)))
                json.dump([b.raw for b in batches], writer, indent='\t')
        except BaseException as ex:
            raise UserError("Work failed: %s", str(ex)) from ex
