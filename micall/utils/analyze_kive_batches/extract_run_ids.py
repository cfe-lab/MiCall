
import json
from pathlib import Path

from micall.utils.new_atomic_file import new_atomic_text_file


def extract_run_ids(input: Path, output: Path) -> None:
    with input.open() as reader:
        runs = json.load(reader)
        with new_atomic_text_file(output) as writer:
            for run in runs:
                run_id = run["id"]
                writer.write(str(run_id) + "\n")
