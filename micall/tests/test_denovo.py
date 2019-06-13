from io import StringIO
from pathlib import Path

from pytest import fixture

from micall.core.denovo import write_genotypes, denovo, DEFAULT_DATABASE, Assembler
from micall.blast_db.make_blast_db import make_blast_db, DEFAULT_PROJECTS


@fixture(scope='session', name='hcv_db')
def check_hcv_db():
    db_path = Path(DEFAULT_DATABASE)
    index_path = db_path.parent / "refs.fasta.nin"
    if not index_path.exists():
        with open(DEFAULT_PROJECTS) as projects_json, \
                open(DEFAULT_DATABASE, 'w') as refs_fasta:
            make_blast_db(projects_json, refs_fasta)
    assert index_path.exists()
    return db_path


def test_write_genotypes_two_sequences(tmpdir, hcv_db):
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


def test_write_genotypes_not_found(tmpdir, hcv_db):
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


def test_write_genotypes_partial_match(tmpdir, hcv_db):
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


def test_write_genotypes_none(tmpdir, hcv_db):
    contigs_fasta = Path(tmpdir) / 'contigs.fasta'
    assert not contigs_fasta.exists()

    contigs_csv = StringIO()
    expected_contigs_csv = """\
genotype,match,contig
"""

    write_genotypes(str(contigs_fasta), contigs_csv)

    assert expected_contigs_csv == contigs_csv.getvalue()


def test_merged_contig(tmpdir, hcv_db):
    contigs_fasta = Path(tmpdir) / 'contigs.fasta'
    assert not contigs_fasta.exists()

    merged_contigs_path = Path(tmpdir) / 'merged_contigs.csv'
    merged_contigs_path.write_text("""\
contig
TGCACAAGACCCAACAACAATACAAGAAAAAGTATAAGGATAGGACCAGGA
""")

    contigs_csv = StringIO()
    expected_contigs_csv = """\
genotype,match,contig
HIV1-C-BR-JX140663-seed,1.0,TGCACAAGACCCAACAACAATACAAGAAAAAGTATAAGGATAGGACCAGGA
"""

    with merged_contigs_path.open() as merged_contigs_csv:
        write_genotypes(str(contigs_fasta),
                        contigs_csv,
                        merged_contigs_csv=merged_contigs_csv)

    assert expected_contigs_csv == contigs_csv.getvalue()


def test_denovo_iva(tmpdir, hcv_db):
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

    denovo(microtest_path / '2150A-HCV_S18_L001_R1_001.fastq',
           microtest_path / '2150A-HCV_S18_L001_R2_001.fastq',
           contigs_csv,
           tmpdir,
           assembler=Assembler.IVA)

    assert expected_contigs_csv == contigs_csv.getvalue()


def test_denovo_savage(tmpdir, hcv_db):
    microtest_path = Path(__file__).parent / 'microtest'
    contigs_csv = StringIO()
    expected_contigs_csv = """\
genotype,match,contig
HCV-1a,1.0,CGACGTGGTTAGCAAGCTCCCCCTGGCCGTGATGGGAAGCTCCTACGGATTCCAATACTCACCAGGA\
CAGCGGGTTGAATTCCTCGTGCAAGCGTGGAAGTCCAAGAAGACCCCGATGGGGTTCTCGTATGATACCCGCTGTTTT\
GACTCCACAGTCACTGAGAGCGACATCCGTACGGAGGAGGCAATTTACCAATGTTGTGACCTGGACCCCCAAGCCCGC\
GTGGCCATCAAGTCCCTCACTGAGAGGCTTTATGTTGGGG
"""

    denovo(microtest_path / '2150A-HCV_S18_L001_R1_001.fastq',
           microtest_path / '2150A-HCV_S18_L001_R2_001.fastq',
           contigs_csv,
           tmpdir,
           assembler=Assembler.SAVAGE)

    assert expected_contigs_csv == contigs_csv.getvalue()
