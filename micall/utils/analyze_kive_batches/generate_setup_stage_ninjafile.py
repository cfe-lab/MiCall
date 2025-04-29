
from typing import Iterable, Iterator
from pathlib import Path

from micall.utils.dir_path import DirPath
from .batch import BatchName
from .ninjamaker import Build, Deref, Command, Statement, Default, Rule, Description, Recipe


def generate_builds(root: DirPath, batches: Iterable[BatchName]) -> Iterator[Build]:
    for batch in batches:
        filename = str(batch) + ".json"
        output = root / filename
        yield Build(outputs=[output],
                    rule="get",
                    inputs=[],
                    bindings=[("batch", str(batch))],
                    )


def generate_statements(root: DirPath, batches: Iterable[BatchName]) -> Iterator[Statement]:
    builds = tuple(generate_builds(root, batches))
    outputs = tuple(x for build in builds for x in build.outputs)

    yield Rule(name="get",
               command=Command(head="micall", arguments=[
                   "analyze_kive_batches",
                   "get_batch_runs",
                   "--batch", Deref("batch"),
                   "--target", Deref("out"),
               ]),
               description=Description("get {}", [Deref("batch")]),
               )
    yield from builds
    yield Default(outputs)


def generate_setup_stage_ninjafile(
        root: DirPath,
        batches: Iterable[BatchName],
        target: Path,
) -> None:
    recipe = Recipe(tuple(generate_statements(root, batches)))
    target.write_text(recipe.compile())
