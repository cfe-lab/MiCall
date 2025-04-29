
from typing import Iterable, Iterator, Tuple
from pathlib import Path

from micall.utils.dir_path import DirPath
from .batch import BatchName
from .ninjamaker import Build, Deref, Command, Statement, Rule, Description, Recipe


def get_outputs(root: DirPath, batches: Iterable[BatchName]) -> Iterator[Tuple[BatchName, Path]]:
    for batch in batches:
        filename = str(batch) + ".json"
        output = root / filename
        yield (batch, output)


def generate_builds(root: DirPath, batches: Iterable[Tuple[BatchName, Path]]) -> Iterator[Build]:
    for batch, output in batches:
        yield Build(outputs=[output],
                    rule="get",
                    inputs=[],
                    bindings=[("batch", str(batch))],
                    )


def generate_statements(root: DirPath, pairs: Iterable[Tuple[BatchName, Path]]) -> Iterator[Statement]:
    yield Rule(name="get",
               command=Command(head="micall", arguments=[
                   "analyze_kive_batches",
                   "get-batch-runs",
                   "--batch", Deref("batch"),
                   "--target", Deref("out"),
               ]),
               description=Description("get {}", [Deref("batch")]),
               )
    yield from generate_builds(root, pairs)


def generate_setup_stage_ninjafile(
        root: DirPath,
        batches: Iterable[BatchName],
        target: Path,
) -> Iterator[Path]:
    outputs_batches_pairs = tuple(get_outputs(root, batches))
    outputs = tuple(output for batch, output in outputs_batches_pairs)
    recipe = Recipe(tuple(generate_statements(root, outputs_batches_pairs)),
                    default=outputs)
    target.write_text(recipe.compile())
    yield from outputs
