
import pytest

from typing import Callable, Literal, Tuple, Iterator, AbstractSet
from collections import defaultdict
import random
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
from micall.utils.referenceless_contig_stitcher_events import EventType
from pathlib import Path
from micall.utils.contig_stitcher_context import ReferencelessStitcherContext
from micall.utils.referenceless_contig_with_aligner import ContigWithAligner, map_overlap
from micall.utils.referenceless_score import Score
import micall.utils.registry as registry
from micall.utils.referenceless_contig_stitcher import \
    referenceless_contig_stitcher_with_ctx, read_contigs
from micall.core.project_config import ProjectConfig


@pytest.fixture(name='projects', scope="session")
def load_projects():
    yield ProjectConfig.loadDefault()



@pytest.fixture(autouse=True)
def disable_kmer_filter(monkeypatch):
    """Set KMER_SIZE=1 to disable kmer filtering for all tests."""
    monkeypatch.setattr("micall.utils.referenceless_contig_stitcher.KMER_SIZE", 1)


@pytest.fixture
def disable_acceptable_prob_check(monkeypatch):
    monkeypatch.setattr("micall.utils.referenceless_contig_stitcher.MIN_MATCHES", 0)


original_map_overlap = map_overlap


def failing_map_overlap(
    self: ContigWithAligner,
    minimum_score: Score,
    relation: Literal["left", "right", "cover"],
    initial_overlap: str,
) -> Iterator[Tuple[int, int]]:
    if "Z" in getattr(self, "seq", ""):
        # Return an empty iterable (no alignments found)
        yield from ()

    # Defer to the real implementation otherwise.
    yield from original_map_overlap(self, minimum_score, relation, initial_overlap)


@pytest.fixture
def force_failing_map_overlap(monkeypatch):
    monkeypatch.setattr("micall.utils.referenceless_contig_with_aligner.map_overlap", failing_map_overlap)
    monkeypatch.setattr("micall.utils.referenceless_contig_stitcher.map_overlap", failing_map_overlap)


@pytest.fixture
def random_fasta_file(tmp_path: Path, projects) -> Callable[[int, int], Tuple[Path, AbstractSet[str]]]:
    ref_name = "HIV1-B-ZA-KP109515-seed"
    ref_full_seq = projects.getReference(ref_name)

    def ret(n_reads: int, random_seed: int) -> Tuple[Path, AbstractSet[str]]:
        root = tmp_path / str(random_seed)
        root.mkdir(parents=True, exist_ok=True)
        converted_fasta_file = root / "converted.fasta"

        rng = random.Random(random_seed)
        ref_seq = ref_full_seq[2000:3001]

        min_length = round((len(ref_seq) / n_reads)**(1/3) * (100 / (len(ref_seq) / 50)**(1/3)))
        max_length = min(len(ref_seq), min_length * 3)

        assert max_length <= len(ref_seq)

        def generate_indexes() -> Iterator[Tuple[int, int]]:
            for i_read in range(n_reads):
                while True:
                    start = rng.randint(0, len(ref_seq) - 1)
                    read_length = rng.randint(min_length, max_length)
                    end = start + read_length
                    if end >= len(ref_seq):
                        continue

                    yield (start, end)
                    break

        indexes = tuple(generate_indexes())
        coverage_map: dict[int, int] = defaultdict(int)
        for (start, end) in indexes:
            for i in range(start, end + 1):
                coverage_map[i] += 1

        sequences = [SeqRecord(Seq.Seq(ref_seq[start:end+1]),
                               description='',
                               id=f'r{start}-{end}.{i}',
                               name=f'r{start}-{end}.{i}')
                     for i, (start, end) in enumerate(indexes)]
        SeqIO.write(sequences, converted_fasta_file, "fasta")

        def generate_islands() -> Iterator[str]:
            old_key: int = min(coverage_map) - 1
            current = ""
            for key in sorted(coverage_map):
                if key != old_key + 1:
                    yield current
                    current = ""
                else:
                    current += ref_seq[key]
                old_key = key
            if current:
                yield current

        ref_seqs = frozenset(generate_islands())

        return (converted_fasta_file, ref_seqs)

    return ret


@pytest.fixture
def log_check(request, tmp_path: Path):
    """
    This fixture verifies that the exact behaviour of the stitcher has not changed.
    It is very brittle by design.

    It overwrites the expected test results after the first run.
    This way, you will only have to commit them in case that
    the stitcher behaviour is actually supposed to change.

    If this test fails due to unexpected log changes, you may revert the logs with:

        git checkout HEAD micall/tests/data
    """

    ReferencelessStitcherContext.set(ReferencelessStitcherContext())
    registry.set(registry.Registry())

    test_name = request.node.name
    log_name = test_name + ".txt"
    pwd = Path(__file__).parent
    logs_dir = pwd / "data" / "referenceless_stitcher_logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    path_to_produced = logs_dir / log_name
    is_rerun = path_to_produced.exists()

    def pretty_print_event(pair: Tuple[int, EventType]) -> str:
        (i, ev) = pair
        message = str(ev)
        name = str(i)
        words = message.split()
        body = '\n  '.join(words)
        header = f"BEGIN {name}"
        footer = f"END {name}"
        return header + '\n  ' + body + '\n' + footer

    def check():
        logs = ReferencelessStitcherContext.get().events
        produced_logs = '\n'.join(map(pretty_print_event, enumerate(logs)))

        if is_rerun:
            with path_to_produced.open() as reader:
                expected_logs: str = reader.read()

        with path_to_produced.open("w") as writer:
            writer.write(produced_logs)

        if is_rerun:
            are_equal = produced_logs == expected_logs
            assert are_equal

    return check


def run_full_pipeline(log_check, tmp_path: Path, converted_fasta_file: Path, ref_seqs: AbstractSet[str]):
    output_fasta_file = tmp_path / "out.fasta"

    # Read the converted FASTA contigs and run the contig stitcher.
    with converted_fasta_file.open("r") as input_handle, \
         output_fasta_file.open("w") as output_handle:
        referenceless_contig_stitcher_with_ctx(input_handle, output_handle)

    log_check()

    with output_fasta_file.open("r") as output_handle:
        stitched_contigs = tuple(read_contigs(output_handle))

    # We check that one of the resulting contigs are in our original sequence.
    for contig in stitched_contigs:
        assert any(tuple(contig.seq in ref for ref in ref_seqs)), "Stitcher produced nonexisting sequences."

    # For a well-sampled read set, we expect the stitiching algorithm to
    # reconstruct the original sequence.
    # We check that one of the resulting contigs matches our original sequence.
    for ref in ref_seqs:
        reconstructed = None

        for contig in stitched_contigs:
            if contig.seq == ref:
                reconstructed = contig.seq
                break

        assert reconstructed is not None, (
            "The contig stitching did not reconstruct the original sequence."
        )

    # Check that only one contig remains.
    count = len(stitched_contigs)
    assert count == len(ref_seqs), (
        f"Expected {len(ref_seqs)} stitched contigs, got {count}."
    )


