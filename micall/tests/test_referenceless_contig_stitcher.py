import pytest
from pathlib import Path
import random
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
from micall.utils.referenceless_contig_stitcher import \
    stitch_consensus, ContigWithAligner, \
    referenceless_contig_stitcher, read_contigs
from micall.utils.contig_stitcher_context import ReferencelessStitcherContext
from micall.utils.fasta_to_fastq import generate_fastq
from micall.utils.fastq_to_fasta import main_typed
from micall.tests.test_remap import load_projects  # activates the "projects" fixture


assert load_projects is not None


@pytest.fixture(autouse=True)
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


def test_full_pipeline(tmp_path: Path, projects):
    fasta_file = tmp_path / "hxb2.fasta"
    fastq_file = tmp_path / "hxb2.fastq"
    converted_fasta_file = tmp_path / "hxb2_converted.fasta"
    output_fasta_file = tmp_path / "hxb2_reconstructed.fasta"
    hxb2_name = "HIV1-B-FR-K03455-seed"
    ref_seq = projects.getReference(hxb2_name)
    ref_seq = ref_seq[4000:5000]  # Consider a part of genome, for speed.

    # Step 1. Write the HXB2 sequence into a FASTA file.
    with fasta_file.open("w") as writer:
        records = [SeqRecord(Seq.Seq(ref_seq),
                             description='',
                             id='HXB2',
                             name='HXB2')]
        SeqIO.write(records, writer, "fasta")

    # Step 2. Use fasta_to_fastq to generate simulated FASTQ reads from the FASTA.
    # We use a fixed random seed for reproducibility.
    rng = random.Random(42)
    # Choose simulation parameters; these should be set so that
    # the reads overlap sufficiently to allow full reconstruction.
    n_reads = 50      # total number of reads to generate
    is_reversed = False
    min_length = 100
    max_length = 900

    generate_fastq(
        fasta=fasta_file,
        fastq=fastq_file,
        n_reads=n_reads,
        is_reversed=is_reversed,
        min_length=min_length,
        max_length=max_length,
        rng=rng
    )

    # Step 3. Convert the simulated FASTQ file back into a FASTA file.
    main_typed(source_fastq=fastq_file, target_fasta=converted_fasta_file)

    # Step 4. Read the converted FASTA contigs and run the contig stitcher.
    with converted_fasta_file.open("r") as input_handle, \
         output_fasta_file.open("w") as output_handle:
        referenceless_contig_stitcher(input_handle, output_handle)

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
