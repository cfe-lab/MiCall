
import json
from typing import Iterator
from pathlib import Path
from kivecli.kiverun import KiveRun

from micall.utils.new_atomic_file import new_atomic_text_file


def extract_run_ids(input: Path, output: Path) -> None:
    def collect_ids() -> Iterator[str]:
        with input.open() as reader:
            runs = json.load(reader)
            for run in runs:
                run = KiveRun.from_json(run)
                yield str(run.id)

    new_content = '\n'.join(collect_ids())

    if output.exists():
        previous_content = output.read_text()
        if previous_content == new_content:
            return

    with new_atomic_text_file(output) as writer:
        writer.write(new_content)
