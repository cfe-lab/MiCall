
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
        stats_output = dir / "stats.json"
        stitcher_output = dir / "stitched"
        input = dir / "info.json"
        yield Build(outputs=[stitcher_output],
                    rule="stitch",
                    inputs=[input],
                    )
        yield Build(outputs=[stats_output],
                    rule="stats",
                    inputs=[input],
                    implicit=[stitcher_output],
                    )

def generate_statements(root: DirPath,
                        runs_txt: Path,
                        ) -> Iterator[Statement]:

    yield Rule(name="stats",
               command=Command.make(
                   "micall",
                   "analyze_kive_batches",
                   "make-stats",
                   "--input", Deref("in"),
                   "--output", Deref("out"),
               ),
               description=Description.make("make stats {}", Deref("in")),
               )

    yield Rule(name="stitch",
               command=Command.make(
                   "micall",
                   "analyze_kive_batches",
                   "stitch-contigs",
                   "--info-file", Deref("in"),
                   "--output", Deref("out"),
               ),
               description=Description.make("stitch contigs {}", Deref("in")),
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

    yield Rule(name="combine_overlaps",
               command=Command.make(
                   "micall",
                   "analyze_kive_batches",
                   "combine-runs-overlaps",
                   "--root", root,
                   "--runs-txt", runs_txt,
                   "--target", Deref("out"),
               ),
               description=Description.make("combine overlaps"),
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

    stats = root / "stats.csv"
    inputs = [input for build in builds for input in build.outputs if build.rule == "stats"]
    overlaps = root / "overlaps.csv"
    aggregated_stats = root / "agg-stats.csv"

    yield Build(rule="combine_stats",
                outputs=[stats],
                inputs=inputs,
                )

    yield Build(rule="combine_overlaps",
                outputs=[overlaps],
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
