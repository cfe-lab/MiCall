
import json
from pathlib import Path

from micall.utils.new_atomic_file import new_atomic_text_file
from .logger import logger


def extract_run_ids(input: Path, output: Path) -> None:
    with input.open() as reader:
        runs = json.load(reader)
        with new_atomic_text_file(output) as writer:
            for run in runs:
                run_id = run["id"]
                state = run["state"]
                if state != "C":
                    logger.warning("Skipping run %s: state %s != 'C'.", run_id, state)
                    continue

                writer.write(str(run_id) + "\n")
