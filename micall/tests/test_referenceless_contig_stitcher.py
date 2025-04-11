import pytest
from pathlib import Path
from typing import Callable, Tuple, Iterator, AbstractSet
from collections import defaultdict
import random
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
from micall.utils.referenceless_contig_stitcher_events import EventType
from micall.utils.referenceless_contig_stitcher import \
    stitch_consensus, ContigWithAligner, \
    referenceless_contig_stitcher_with_ctx, read_contigs
from micall.utils.contig_stitcher_context import ReferencelessStitcherContext
import micall.utils.registry as registry


@pytest.fixture
def disable_acceptable_prob_check(monkeypatch):
    monkeypatch.setattr("micall.utils.referenceless_contig_stitcher.ACCEPTABLE_STITCHING_SCORE", 0)


TTT = 40 * 'T'


@pytest.mark.parametrize(
    "seqs, expected",
    [
        #
        # Singletons
        #

        (('AAAAA' + TTT, TTT + 'GGGGG'), ('AAAAA' + TTT + 'GGGGG',)),
        (('AAAAA' + TTT, TTT + 'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG'), ('AAAAA' + TTT + 'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG',)),
        (('AAAAA' + TTT + 'CCCCCCCCCCCCCCCCCCCC', TTT + 'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG'), ('AAAAA' + TTT + 'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG',)),
        (('AAAAA' + TTT, 'GGGGG' + TTT), ("AAAAA" + TTT,)),
        (('GGGGG' + TTT, 'CCCCCAAAAA' + TTT), ('CCCCCAAAAA' + TTT,)),
        (('AAAAA' + TTT, 'GGGGG' + TTT), ('AAAAA' + TTT,)),
        ((TTT + 'AAAAA', TTT + 'GGGGG'), (TTT + 'AAAAA',)),
        (('AAAAAAAAAAAAAA' + TTT, 'GGGGG' + TTT), ('AAAAAAAAAAAAAA' + TTT,)),
        (('GGGGGGGGGGGGGG' + TTT, 'AAAAA' + TTT), ('GGGGGGGGGGGGGG' + TTT,)),
        (('AAAAA' + TTT, 'GGGGGGGGGGGGGG' + TTT), ('GGGGGGGGGGGGGG' + TTT,)),
        (('GGGGG' + TTT, 'AAAAAAAAAAAAAA' + TTT), ('AAAAAAAAAAAAAA' + TTT,)),
        (('GGGGG' + TTT + 'AAAAAAAAAAAAA', 'CC' + TTT + 'CC'), ('GGGGG' + TTT + 'AAAAAAAAAAAAA',)),
        (('GGGGG' + TTT + 'AAAAAAAAAAAAA', TTT + 'CC'), ('GGGGG' + TTT + 'AAAAAAAAAAAAA',)),
        (('GGGGG' + TTT + 'AAAAAAAAAAAAA', 'CC' + TTT), ('GGGGG' + TTT + 'AAAAAAAAAAAAA',)),
        (('CC' + TTT + 'CC', 'GGGGG' + TTT + 'AAAAAAAAAAAAA'), ('GGGGG' + TTT + 'AAAAAAAAAAAAA',)),
        ((TTT + 'CC', 'GGGGG' + TTT + 'AAAAAAAAAAAAA'), ('GGGGG' + TTT + 'AAAAAAAAAAAAA',)),
        (('CC' + TTT, 'GGGGG' + TTT + 'AAAAAAAAAAAAA'), ('GGGGG' + TTT + 'AAAAAAAAAAAAA',)),
        ((TTT, 'GGGGG' + TTT + 'AAAAAAAAAAAAA'), ('GGGGG' + TTT + 'AAAAAAAAAAAAA',)),
        (('GGGGG' + TTT + 'AAAAAAAAAAAAA', TTT), ('GGGGG' + TTT + 'AAAAAAAAAAAAA',)),
        (('AAAA', 'GGGG'), ('AAAA', 'GGGG',)),
        (('AAAAT', 'TGGGG'), ('AAAAT', 'TGGGG')),

        #
        # Multiple.
        #
        (('AAA' + 'T' * 40,  'T' * 40 + 'GGG' + 'Y' * 40, 'Y' * 40 + 'CCC' + 'M' * 40, 'M' * 40 + 'GGG'),
         ('AAA' + 40 * 'T' + 'GGG' + 'Y' * 40, 'Y' * 40 + 'CCC' + 'M' * 40, 'M' * 40 + 'GGG')),
        (('AAA' + 'T' * 40,
          'T' * 40 + 'GGG' + 'A' * 40,
          'A' * 40 + 'CCC' + 'T' * 40,
          'T' * 40 + 'GGG'),
         ('A' * 40 + 'CCC' + 'T' * 40 + 'GGG' + 'A' * 40,)),
        (('AAA',), ('AAA',)),
        ((), ()),

    ],
)
def test_stitch_simple_cases(seqs, expected):
    contigs = [ContigWithAligner(None, seq) for seq in seqs]
    with ReferencelessStitcherContext.fresh():
        consenses = tuple(sorted(contig.seq for contig in stitch_consensus(contigs)))
    assert consenses == tuple(sorted(expected))


@pytest.fixture
def random_fasta_file(tmp_path: Path) -> Callable[[int, int], Tuple[Path, AbstractSet[str]]]:
    def ret(n_reads: int, random_seed: int) -> Tuple[Path, AbstractSet[str]]:
        root = tmp_path / str(random_seed)
        root.mkdir(parents=True, exist_ok=True)
        converted_fasta_file = root / "converted.fasta"

        rng = random.Random(random_seed)
        ref_seq = ''.join(rng.choices("ACTG", k=1000))

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


@pytest.mark.parametrize("random_seed", [1, 2, 3, 5, 7, 11, 13, 17, 1337])
def test_full_pipeline(log_check, tmp_path: Path, random_fasta_file, random_seed: int):
    assert not ReferencelessStitcherContext.get().is_debug2
    converted_fasta_file, ref_seqs = random_fasta_file(50, random_seed)
    run_full_pipeline(log_check, tmp_path, converted_fasta_file, ref_seqs)


# TODO: ensure that every random seed can be stitched.
@pytest.mark.parametrize("random_seed", sorted(set(range(50)).difference([0, 2, 3, 6, 7, 8, 13, 14, 17, 27, 33, 35])))
def test_full_pipeline_small_values(log_check, tmp_path: Path, random_fasta_file, random_seed: int, monkeypatch):
    monkeypatch.setattr("micall.utils.referenceless_contig_stitcher.MAX_ALTERNATIVES", 1)
    assert not ReferencelessStitcherContext.get().is_debug2
    ReferencelessStitcherContext.get().is_debug2 = True
    converted_fasta_file, ref_seqs = random_fasta_file(6, random_seed)
    run_full_pipeline(log_check, tmp_path, converted_fasta_file, ref_seqs)


# TODO: ensure that every random seed can be stitched.
@pytest.mark.parametrize("random_seed", sorted(set(range(999)).difference([2, 8, 14, 15, 17, 27, 29, 33, 36, 49, 51, 52, 56, 60, 62, 63, 68, 69, 71, 76, 81, 82, 84, 92, 95, 104, 112, 117, 124, 134, 141, 145, 158, 159, 199, 202, 232, 235, 236, 239, 240, 253, 256, 257, 267, 271, 272, 283, 285, 294, 310, 312, 314, 318, 320, 334, 337, 338, 350, 365, 375, 377, 378, 383, 386, 389, 392, 399, 404, 426, 427, 437, 444, 445, 451, 453, 458, 459, 461, 463, 465, 471, 474, 477, 487, 499, 501, 507, 526, 530, 533, 534, 536, 541, 544, 546, 552, 554, 561, 571, 573, 581, 582, 589, 592, 593, 600, 601, 625, 630, 631, 633, 634, 635, 640, 648, 651, 655, 660, 663, 668, 669, 671, 673, 685, 693, 700, 706, 719, 722, 725, 742, 746, 750, 754, 756, 763, 765, 769, 773, 780, 781, 782, 785, 789, 791, 797, 803, 804, 805, 830, 839, 843, 850, 851, 865, 881, 884, 886, 902, 909, 913, 915, 916, 917, 918, 923, 928, 937, 941, 949, 957, 967, 972, 973, 979, 981, 983, 997])))
def test_full_pipeline_tiny_values(log_check, tmp_path: Path, random_fasta_file, random_seed: int, monkeypatch):
    monkeypatch.setattr("micall.utils.referenceless_contig_stitcher.MAX_ALTERNATIVES", 1)
    assert not ReferencelessStitcherContext.get().is_debug2
    ReferencelessStitcherContext.get().is_debug2 = True
    converted_fasta_file, ref_seqs = random_fasta_file(3, random_seed)
    run_full_pipeline(log_check, tmp_path, converted_fasta_file, ref_seqs)
