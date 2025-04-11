import pytest
from pathlib import Path
from typing import Callable, Tuple
import random
import math
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
from micall.utils.referenceless_contig_stitcher_events import EventType
from micall.utils.referenceless_contig_stitcher import \
    stitch_consensus, ContigWithAligner, \
    referenceless_contig_stitcher_with_ctx, read_contigs
from micall.utils.contig_stitcher_context import ReferencelessStitcherContext, context
from micall.utils.fasta_to_fastq import generate_fastq
from micall.utils.fastq_to_fasta import main_typed
import micall.utils.registry as registry
from micall.tests.test_remap import load_projects  # activates the "projects" fixture


assert load_projects is not None


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
        (('AAAAA' + TTT, 'GGGGG' + TTT), ("GGGGG" + TTT,)),
        (('GGGGG' + TTT, 'CCCCCAAAAA' + TTT), ('CCCCCAAAAA' + TTT,)),
        (('AAAAA' + TTT, 'GGGGG' + TTT), ('GGGGG' + TTT,)),
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
def random_fasta_file(tmp_path: Path, projects) -> Callable[[int, int], Tuple[Path, str]]:
    hxb2_name = "HIV1-B-FR-K03455-seed"
    ref_seq = projects.getReference(hxb2_name)
    ref_seq = ref_seq[4000:5000]  # Consider a part of genome, for speed.

    def ret(n_reads: int, random_seed: int) -> Tuple[Path, str]:
        root = tmp_path / str(random_seed)
        root.mkdir(parents=True, exist_ok=True)

        fasta_file = root / "hxb2.fasta"
        fastq_file = root / "hxb2.fastq"
        converted_fasta_file = root / "hxb2_converted.fasta"

        # Step 1. Write the HXB2 sequence into a FASTA file.
        with fasta_file.open("w") as writer:
            records = [SeqRecord(Seq.Seq(ref_seq))]
            SeqIO.write(records, writer, "fasta")

        # Step 2. Use fasta_to_fastq to generate simulated FASTQ reads from the FASTA.
        # We use a fixed random seed for reproducibility.
        rng = random.Random(random_seed)
        # Choose simulation parameters; these should be set so that
        # the reads overlap sufficiently to allow full reconstruction.
        is_reversed = False
        min_length = round(math.sqrt(len(ref_seq) / n_reads) * (100 / math.sqrt(len(ref_seq) / 50)))
        max_length = min(len(ref_seq), min_length * 3)

        generate_fastq(
            fasta=fasta_file,
            fastq=fastq_file,
            n_reads=n_reads,
            is_reversed=is_reversed,
            min_length=min_length,
            max_length=max_length,
            rng=rng,
        )

        # Step 3. Convert the simulated FASTQ file back into a FASTA file.
        main_typed(source_fastq=fastq_file, target_fasta=converted_fasta_file)
        return converted_fasta_file, ref_seq

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

    context.set(ReferencelessStitcherContext.make())
    registry.set(registry.Registry())

    test_name = request.node.name
    log_name = test_name + ".txt"
    pwd = Path(__file__).parent
    logs_dir = pwd / "data" / "referenceless_stitcher_logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    path_to_produced = logs_dir / log_name
    is_rerun = path_to_produced.exists()

    def pretty_print_event(ev: EventType) -> str:
        message = str(ev)
        words = message.split()
        first = str(ev.__class__.__name__)
        rest = words
        return first + '\n  ' + '\n  '.join(rest)

    def check():
        logs = context.get().events
        produced_logs = '\n'.join(map(pretty_print_event, logs))

        if is_rerun:
            with path_to_produced.open() as reader:
                expected_logs: str = reader.read()

        with path_to_produced.open("w") as writer:
            writer.write(produced_logs)

        if is_rerun:
            are_equal = produced_logs == expected_logs
            assert are_equal

    return check


def run_full_pipeline(log_check, tmp_path: Path, converted_fasta_file: Path, ref_seq: str):
    output_fasta_file = tmp_path / "out.fasta"

    # Read the converted FASTA contigs and run the contig stitcher.
    with converted_fasta_file.open("r") as input_handle, \
         output_fasta_file.open("w") as output_handle:
        referenceless_contig_stitcher_with_ctx(input_handle, output_handle)

    log_check()

    with output_fasta_file.open("r") as output_handle:
        stitched_contigs = tuple(read_contigs(output_handle))

    # We check that one of the resulting contigs are in our original HXB2 sequence.
    for contig in stitched_contigs:
        assert contig.seq in ref_seq, "Stitcher produced nonexisting sequences."

    # For a well-sampled read set, we expect the stitiching algorithm to
    # reconstruct the original sequence.
    # We check that one of the resulting contigs matches our original HXB2 sequence.
    reconstructed = None
    for contig in stitched_contigs:
        if contig.seq == ref_seq:
            reconstructed = contig.seq
            break

    assert reconstructed is not None, (
        "The contig stitching did not reconstruct the original HXB2 sequence."
    )

    # Check that only one contig remains.
    count = len(stitched_contigs)
    assert count == 1, (
        "Expected one stitched contig; got "
        f"{count} contigs."
    )


@pytest.mark.parametrize("random_seed", [1, 2, 3, 5, 7, 11, 13, 17, 42, 1337])
def test_full_pipeline(log_check, tmp_path: Path, random_fasta_file, random_seed: int):
    converted_fasta_file, ref_seq = random_fasta_file(50, random_seed)
    run_full_pipeline(log_check, tmp_path, converted_fasta_file, ref_seq)


@pytest.mark.parametrize("random_seed", set(range(200)).difference([2, 20, 22, 24, 26, 31, 62, 79, 82, 86, 97, 101, 113, 124, 143, 146, 179, 183, 195]))
def test_full_pipeline_small_values(log_check, tmp_path: Path, random_fasta_file, random_seed: int, monkeypatch):
    monkeypatch.setattr("micall.utils.referenceless_contig_stitcher.MAX_ALTERNATIVES", 1)
    converted_fasta_file, ref_seq = random_fasta_file(6, random_seed)
    run_full_pipeline(log_check, tmp_path, converted_fasta_file, ref_seq)
