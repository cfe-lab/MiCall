
from typing import Iterable, Iterator, Tuple
from pathlib import Path

from micall.utils.dir_path import DirPath
from .batch import BatchName
from .ninjamaker import Build, Deref, Command, Statement, Rule, Description, Recipe


def get_outputs(root: DirPath, batches: Iterable[BatchName]) -> Iterator[Tuple[BatchName, Path]]:
    for batch in batches:
        filename = str(batch) + ".json"
        output = root / "batches" / filename
        yield (batch, output)


def generate_builds(root: DirPath,
                    runs_json: Path,
                    pairs: Iterable[Tuple[BatchName, Path]],
                    ) -> Iterator[Build]:

    outputs = tuple(output for batch, output in pairs)
    yield Build(rule="combine",
                outputs=[runs_json],
                inputs=outputs,
                )

    for batch, output in pairs:
        yield Build(outputs=[output],
                    rule="get",
                    inputs=[],
                    bindings=[("batch", str(batch))],
                    )


def generate_statements(root: DirPath,
                        runs_json: Path,
                        pairs: Iterable[Tuple[BatchName, Path]],
                        ) -> Iterator[Statement]:

    yield Rule(name="get",
               command=Command.make(
                   "micall",
                   "analyze_kive_batches",
                   "get-batch",
                   "--batch", Deref("batch"),
                   "--target", Deref("out"),
               ),
               description=Description("get {}", [Deref("batch")]),
               )

    yield Rule(name="combine",
               command=Command.make(
                   "micall",
                   "analyze_kive_batches",
                   "combine-batches-runs",
                   "--batches", Deref("in"),
                   "--target", Deref("out"),
               ),
               description=Description("combine {}", [Deref("in")]),
               )

    yield from generate_builds(root, runs_json, pairs)


def generate_setup_stage_ninjafile(
        root: DirPath,
        batches: Iterable[BatchName],
        target: Path,
        runs_json: Path,
) -> None:
    outputs_batches_pairs = tuple(get_outputs(root, batches))
    recipe = Recipe(tuple(generate_statements(root, runs_json, outputs_batches_pairs)),
                    default=[runs_json])
    target.write_text(recipe.compile())
