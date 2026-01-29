
from typing import Iterable
from pathlib import Path
import json


def combine_batches_runs(batches: Iterable[Path], target: Path) -> None:
    result = []

    for output in batches:
        with output.open() as reader:
            batches_json = json.load(reader)
            for batch_json in batches_json:
                runs = batch_json["runs"]
                result.extend(runs)

    with target.open("w") as writer:
        json.dump(result, writer, indent='\t')
