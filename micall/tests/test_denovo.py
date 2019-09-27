from io import StringIO
from pathlib import Path

from pytest import fixture

from micall.core.denovo import write_contig_refs, denovo, DEFAULT_DATABASE, Assembler
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


def test_write_contig_refs_two_sequences(tmpdir, hcv_db):
    contigs_fasta = Path(tmpdir) / "contigs.fasta"
    contigs_fasta.write_text("""\
>foo
TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA
>bar
CAGGGCTCCAGGACTGCACCATGCTCGTGTGTGGCGACGAC
""")
    contigs_csv = StringIO()
    expected_contigs_csv = """\
ref,match,group_ref,contig
HCV-1a,1.0,HCV-1a,TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA
HCV-1a,1.0,HCV-1a,CAGGGCTCCAGGACTGCACCATGCTCGTGTGTGGCGACGAC
"""

    write_contig_refs(str(contigs_fasta), contigs_csv)

    assert expected_contigs_csv == contigs_csv.getvalue()


def test_write_contig_refs_two_groups(tmpdir, hcv_db):
    contigs_fasta = Path(tmpdir) / "contigs.fasta"
    contigs_fasta.write_text("""\
>foo
ACCCGCCCCTAATAGGGGCGACACTCCGCCATGAATC
>bar
ACCATGGATCACTCCCCTGTGAGGAACTACTGTCTT
>baz
TGCAATGACAGCTTACAGACGGGTTTCCTCGCTTCCTTGTTTTACACCCA
""")
    contigs_csv = StringIO()
    expected_contigs_csv = """\
ref,match,group_ref,contig
HCV-2a,1.0,HCV-2b,ACCCGCCCCTAATAGGGGCGACACTCCGCCATGAATC
HCV-1g,1.0,HCV-1g,ACCATGGATCACTCCCCTGTGAGGAACTACTGTCTT
HCV-2b,1.0,HCV-2b,TGCAATGACAGCTTACAGACGGGTTTCCTCGCTTCCTTGTTTTACACCCA
"""

    write_contig_refs(str(contigs_fasta), contigs_csv)

    assert expected_contigs_csv == contigs_csv.getvalue()


def test_write_contig_refs_not_found(tmpdir, hcv_db):
    contigs_fasta = Path(tmpdir) / "contigs.fasta"
    contigs_fasta.write_text("""\
>foo
CATCACATAGGAGA
""")
    contigs_csv = StringIO()
    expected_contigs_csv = """\
ref,match,group_ref,contig
unknown,0,,CATCACATAGGAGA
"""

    write_contig_refs(str(contigs_fasta), contigs_csv)

    assert expected_contigs_csv == contigs_csv.getvalue()


def test_write_contig_refs_partial_match(tmpdir, hcv_db):
    contigs_fasta = Path(tmpdir) / "contigs.fasta"
    contigs_fasta.write_text("""\
>foo
TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA
>bar
CATCACATAGGAGACAGGGCTCCAGGACTGCACCATGCTCGTGTGTGGCGACGAC
""")
    contigs_csv = StringIO()
    expected_contigs_csv = """\
ref,match,group_ref,contig
HCV-1a,1.0,HCV-1a,TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA
HCV-1a,0.75,HCV-1a,CATCACATAGGAGACAGGGCTCCAGGACTGCACCATGCTCGTGTGTGGCGACGAC
"""

    write_contig_refs(str(contigs_fasta), contigs_csv)

    assert expected_contigs_csv == contigs_csv.getvalue()


def test_write_blast(tmpdir, hcv_db):
    contigs_fasta = Path(tmpdir) / "contigs.fasta"
    contigs_fasta.write_text("""\
>foo
TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA
>bar
CATCACATAGGAGACAGGGCTCCAGGACTGCACCATGCTCGTGTGTGGCGACGAC
""")
    contigs_csv = StringIO()
    expected_contigs_csv = """\
ref,match,group_ref,contig
HCV-1a,1.0,HCV-1a,TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA
HCV-1a,0.75,HCV-1a,CATCACATAGGAGACAGGGCTCCAGGACTGCACCATGCTCGTGTGTGGCGACGAC
"""
    blast_csv = StringIO()
    expected_blast_csv = """\
contig_name,ref_name,score,match,pident,start,end,ref_start,ref_end
bar,HCV-1g,37,0.67,100,19,55,8506,8542
bar,HCV-1a,41,0.75,100,15,55,8518,8558
foo,HCV-1a,41,1.0,100,1,41,8187,8227
"""

    write_contig_refs(str(contigs_fasta), contigs_csv, blast_csv=blast_csv)

    assert expected_contigs_csv == contigs_csv.getvalue()
    assert expected_blast_csv == blast_csv.getvalue()


def test_write_contig_refs_none(tmpdir, hcv_db):
    contigs_fasta = Path(tmpdir) / 'contigs.fasta'
    assert not contigs_fasta.exists()

    contigs_csv = StringIO()
    expected_contigs_csv = """\
ref,match,group_ref,contig
"""

    write_contig_refs(str(contigs_fasta), contigs_csv)

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
ref,match,group_ref,contig
HIV1-C-BR-JX140663-seed,1.0,HIV1-C-BR-JX140663-seed,TGCACAAGACCCAACAACAATACAAGAAAAAGTATAAGGATAGGACCAGGA
"""

    with merged_contigs_path.open() as merged_contigs_csv:
        write_contig_refs(str(contigs_fasta),
                          contigs_csv,
                          merged_contigs_csv=merged_contigs_csv)

    assert expected_contigs_csv == contigs_csv.getvalue()


def test_denovo_iva(tmpdir, hcv_db):
    microtest_path = Path(__file__).parent / 'microtest'
    contigs_csv = StringIO()
    expected_contigs_csv = """\
ref,match,group_ref,contig
HCV-2a,1.0,HCV-2a,TCGACGCTCATTACGACTCAGTCTTAAAGGACATCAAGCTAGCGGCCTCCAAGGTCTCCG\
CGAGGCTCCTCACCCTGGAGGAGGCATGCCAGCTAACTCCACCCCATTCTGCAAGATCCAAATATGGGTTTGGGGCTA\
AGGAGGTGCGCAGCTTGTCCGGGAGGGCCGTTAACCACATCAAGTCCGTGTGGAAGGACCTCCTGGAAGACTCACAAA\
CACCAATTCCCACAACCATCATGGCCAAAAATGAAGTGTTCTGCGTGGACCCCACCAAGGGGGGTAAGAAAGCAGCTC\
GCCTCATCGTTTACCCTGACCTCGG
"""

    denovo(microtest_path / '2160A-HCV_S19_L001_R1_001.fastq',
           microtest_path / '2160A-HCV_S19_L001_R2_001.fastq',
           contigs_csv,
           tmpdir,
           assembler=Assembler.IVA)

    assert expected_contigs_csv == contigs_csv.getvalue()


def test_denovo_savage(tmpdir, hcv_db):
    microtest_path = Path(__file__).parent / 'microtest'
    contigs_csv = StringIO()
    expected_contigs_csv = """\
ref,match,group_ref,contig
HCV-1a,1.0,HCV-1a,CGACGTGGTTAGCAAGCTCCCCCTGGCCGTGATGGGAAGCTCCTACGGATTCCAATACTC\
ACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAAGTCCAAGAAGACCCCGATGGGGTTCTCGTATGATACCCG\
CTGTTTTGACTCCACAGTCACTGAGAGCGACATCCGTACGGAGGAGGCAATTTACCAATGTTGTGACCTGGACCCCCA\
AGCCCGCGTGGCCATCAAGTCCCTCACTGAGAGGCTTTATGTTGGGG
"""

    denovo(microtest_path / '2150A-HCV_S18_L001_R1_001.fastq',
           microtest_path / '2150A-HCV_S18_L001_R2_001.fastq',
           contigs_csv,
           tmpdir,
           assembler=Assembler.SAVAGE)

    assert expected_contigs_csv == contigs_csv.getvalue()
