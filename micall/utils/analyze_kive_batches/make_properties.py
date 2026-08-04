
import csv
import tomllib
from pathlib import Path
from typing import Iterator

from micall.utils.new_atomic_file import new_atomic_text_file


def make_properties(input: Path, output: Path) -> None:
    with input.open("rb") as reader:
        obj = tomllib.load(reader)

    def collect_properties() -> Iterator[str]:
        seen = set()
        yield "app"  # The index key.
        for app, props in obj.items():
            for key in props.keys():
                if key not in seen:
                    seen.add(key)
                    yield key

    fieldnames = tuple(collect_properties())
    with new_atomic_text_file(output) as writer:
        dwriter = csv.DictWriter(writer, fieldnames=fieldnames)
        dwriter.writeheader()

        for app, props in obj.items():
            fields = {
                key: props.get(key)
                for key in fieldnames
                if key != "app"
            }
            fields["app"] = app
            dwriter.writerow(fields)
