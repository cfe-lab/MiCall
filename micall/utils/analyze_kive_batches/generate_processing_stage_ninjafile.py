
from typing import Iterator
from pathlib import Path
import json

from micall.utils.dir_path import DirPath
from .ninjamaker import Build, Deref, Command, Statement, Rule, Description, Recipe


def generate_builds(root: DirPath,
                    runs_json: Path,
                    ) -> Iterator[Build]:

    with runs_json.open() as reader:
        runs = json.load(reader)
        run_ids = tuple(run["id"] for run in runs)

    for run_id in run_ids:
        dir = root / "runs" / str(run_id)
        output = dir / "stats.json"
        input = dir / "info.json"
        yield Build(outputs=[output],
                    rule="analyze",
                    inputs=[input],
                    )

def generate_statements(root: DirPath,
                        runs_json: Path,
                        ) -> Iterator[Statement]:

    yield Rule(name="analyze",
               command=Command.make(
                   "micall",
                   "analyze_kive_batches",
                   "make-stats-1",
                   "--input", Deref("in"),
                   "--output", Deref("out"),
               ),
               description=Description.make("analyze {}", Deref("in")),
               )

    yield Rule(name="combine",
               command=Command.make(
                   "micall",
                   "analyze_kive_batches",
                   "combine-runs-stats",
                   "--root", root,
                   "--runs-json", runs_json,
                   "--target", Deref("out"),
               ),
               description=Description.make("combine"),
               )

    builds = tuple(generate_builds(root, runs_json))

    inputs = [input for build in builds for input in build.outputs]
    output = root / "stats.csv"
    yield Build(rule="combine",
                outputs=[output],
                inputs=inputs,
                )

    yield from builds


def generate_processing_stage_ninjafile(
        root: DirPath,
        target: Path,
        runs_json: Path,
) -> None:
    statements = tuple(generate_statements(root, runs_json))
    builds = [s for s in statements if isinstance(s, Build)]
    outputs = [o for build in builds for o in build.outputs]
    recipe = Recipe(statements, default=outputs)
    target.write_text(recipe.compile())
