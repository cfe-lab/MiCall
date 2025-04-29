
from typing import Iterator
from pathlib import Path

from micall.utils.dir_path import DirPath
from .ninjamaker import Build, Deref, Command, Statement, Rule, Description, Recipe


def generate_builds(root: DirPath,
                    runs_txt: Path,
                    ) -> Iterator[Build]:

    run_ids = runs_txt.read_text().splitlines()
    for run_id in run_ids:
        dir = root / "runs" / str(run_id)
        output = dir / "stats.json"
        input = dir / "info.json"
        yield Build(outputs=[output],
                    rule="stats",
                    inputs=[input],
                    )

def generate_statements(root: DirPath,
                        runs_txt: Path,
                        ) -> Iterator[Statement]:

    yield Rule(name="stats",
               command=Command.make(
                   "micall",
                   "analyze_kive_batches",
                   "make-stats-1",
                   "--input", Deref("in"),
                   "--output", Deref("out"),
               ),
               description=Description.make("analyze {}", Deref("in")),
               )

    yield Rule(name="combine_stats",
               command=Command.make(
                   "micall",
                   "analyze_kive_batches",
                   "combine-runs-stats",
                   "--root", root,
                   "--runs-txt", runs_txt,
                   "--target", Deref("out"),
               ),
               description=Description.make("combine stats"),
               )

    yield Rule(name="aggregate_stats",
               command=Command.make(
                   "micall",
                   "analyze_kive_batches",
                   "aggregate-runs-stats",
                   "--input", Deref("in"),
                   "--output", Deref("out"),
               ),
               description=Description.make("aggregate stats"),
               )

    builds = tuple(generate_builds(root, runs_txt))

    inputs = [input for build in builds for input in build.outputs]
    stats = root / "stats.csv"
    aggregated_stats = root / "agg-stats.csv"

    yield Build(rule="combine_stats",
                outputs=[stats],
                inputs=inputs,
                )

    yield Build(rule="aggregate_stats",
                outputs=[aggregated_stats],
                inputs=[stats],
                )

    yield from builds


def generate_processing_stage_ninjafile(
        root: DirPath,
        target: Path,
        runs_txt: Path,
) -> None:
    statements = tuple(generate_statements(root, runs_txt))
    builds = [s for s in statements if isinstance(s, Build)]
    outputs = [o for build in builds for o in build.outputs]
    recipe = Recipe(statements, default=outputs)
    target.write_text(recipe.compile())
