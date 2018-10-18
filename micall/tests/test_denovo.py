from io import StringIO
from pathlib import Path

from micall.core.denovo import write_genotypes, denovo


def test_write_genotypes_two_sequences(tmpdir):
    contigs_fasta = Path(tmpdir) / "contigs.fasta"
    contigs_fasta.write_text("""\
>foo
TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA
>bar
CAGGGCTCCAGGACTGCACCATGCTCGTGTGTGGCGACGAC
""")
    contigs_csv = StringIO()
    expected_contigs_csv = """\
genotype,match,contig
HCV-1a,1.0,TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA
HCV-1a,1.0,CAGGGCTCCAGGACTGCACCATGCTCGTGTGTGGCGACGAC
"""

    write_genotypes(str(contigs_fasta), contigs_csv)

    assert expected_contigs_csv == contigs_csv.getvalue()


def test_write_genotypes_not_found(tmpdir):
    contigs_fasta = Path(tmpdir) / "contigs.fasta"
    contigs_fasta.write_text("""\
>foo
CATCACATAGGAGA
""")
    contigs_csv = StringIO()
    expected_contigs_csv = """\
genotype,match,contig
unknown,0,CATCACATAGGAGA
"""

    write_genotypes(str(contigs_fasta), contigs_csv)

    assert expected_contigs_csv == contigs_csv.getvalue()


def test_write_genotypes_partial_match(tmpdir):
    contigs_fasta = Path(tmpdir) / "contigs.fasta"
    contigs_fasta.write_text("""\
>foo
TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA
>bar
CATCACATAGGAGACAGGGCTCCAGGACTGCACCATGCTCGTGTGTGGCGACGAC
""")
    contigs_csv = StringIO()
    expected_contigs_csv = """\
genotype,match,contig
HCV-1a,1.0,TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA
HCV-1a,0.75,CATCACATAGGAGACAGGGCTCCAGGACTGCACCATGCTCGTGTGTGGCGACGAC
"""

    write_genotypes(str(contigs_fasta), contigs_csv)

    assert expected_contigs_csv == contigs_csv.getvalue()


def test_write_genotypes_none(tmpdir):
    contigs_fasta = None
    contigs_csv = StringIO()
    expected_contigs_csv = """\
genotype,match,contig
"""

    write_genotypes(contigs_fasta, contigs_csv)

    assert expected_contigs_csv == contigs_csv.getvalue()


def test_denovo(tmpdir):
    microtest_path = Path(__file__).parent / 'microtest'
    contigs_csv = StringIO()
    expected_contigs_csv = """\
genotype,match,contig
HCV-1a,1.0,TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAAGTCCAAGAAGACCCCGATGGGGTTCT\
CGTATGATACCCGCTGTTTTGACTCCACAGTCACTGAGAGCGACATCCGTACGGAGGAGGCAATTTACCAATGTTGTG\
ACCTGGACCCCCAAGCCCGCGTGGCCATCAAGTCCCTCACTGAGAGGCTTTATGTTGGGGGCCCTCTTACCAATTCAA\
GGGGGGAAAACTGCGGCTACCGCAGGTGCCGCGCGAGCGGCGTACTGACAACTAGCTGTGGTAACACCCTCACTTGCT\
ACATCAAGGCCCGGGCAGCCTGTCGAGCCGCAGGGCTCCAGGACTGCACCATGCTCGTGTGTGGCGACGACTTAGTCG\
TTATCTGTGAAAGTGCGGGGGTCCAGGAGGACG
"""

    denovo(microtest_path / '2140A-HCV_S17_L001_R1_001.fastq',
           microtest_path / '2140A-HCV_S17_L001_R2_001.fastq',
           contigs_csv,
           tmpdir)

    assert expected_contigs_csv == contigs_csv.getvalue()
