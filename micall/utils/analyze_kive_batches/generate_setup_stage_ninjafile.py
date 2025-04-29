
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
                    runs_txt: Path,
                    pairs: Iterable[Tuple[BatchName, Path]],
                    ) -> Iterator[Build]:

    outputs = tuple(output for batch, output in pairs)
    yield Build(rule="combine",
                outputs=[runs_json],
                inputs=outputs,
                )

    yield Build(rule="extract",
                outputs=[runs_txt],
                inputs=[runs_json],
                )

    for batch, output in pairs:
        yield Build(outputs=[output],
                    rule="get",
                    inputs=[],
                    bindings=[("batch", str(batch))],
                    )


def generate_statements(root: DirPath,
                        runs_json: Path,
                        runs_txt: Path,
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

    yield Rule(name="extract",
               command=Command.make(
                   "micall",
                   "analyze_kive_batches",
                   "extract-run-ids",
                   "--input", Deref("in"),
                   "--output", Deref("out"),
               ),
               description=Description("extract {}", [Deref("in")]),
               )

    yield from generate_builds(root, runs_json, runs_txt, pairs)


def generate_setup_stage_ninjafile(
        root: DirPath,
        batches: Iterable[BatchName],
        target: Path,
        runs_json: Path,
        runs_txt: Path,
) -> None:
    outputs_batches_pairs = tuple(get_outputs(root, batches))
    recipe = Recipe(tuple(generate_statements(root, runs_json, runs_txt, outputs_batches_pairs)),
                    default=[runs_txt])
    target.write_text(recipe.compile())
