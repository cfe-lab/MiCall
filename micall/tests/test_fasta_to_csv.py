from io import StringIO
from pathlib import Path

from Bio import SeqIO
import pytest

from micall.utils.fasta_to_csv import default_database, genotype, fasta_to_csv
from micall.blast_db.make_blast_db import make_blast_db
from micall.utils.externals import ProjectsFile


with ProjectsFile().path() as projects_file_path:
    DEFAULT_PROJECTS = str(projects_file_path)


@pytest.fixture(scope='session')
def DEFAULT_DATABASE():
    with default_database() as ret:
        yield ret


@pytest.fixture(scope='session', name='hcv_db')
def check_hcv_db(DEFAULT_DATABASE):
    db_path = DEFAULT_DATABASE
    index_path = db_path.parent / "refs.fasta.nin"
    build_needed = not index_path.exists()
    if not build_needed:
        projects_date = Path(DEFAULT_PROJECTS).stat().st_mtime
        index_date = index_path.stat().st_mtime
        build_needed = index_date < projects_date
    if build_needed:
        with open(DEFAULT_PROJECTS) as projects_json, \
                open(DEFAULT_DATABASE, 'w') as refs_fasta:
            make_blast_db(projects_json, refs_fasta)
    assert index_path.exists()
    return db_path


def test_make_blast_db_excludes_hivgha(hcv_db, DEFAULT_DATABASE):
    fasta_path = DEFAULT_DATABASE
    with fasta_path.open() as f:
        for reference in SeqIO.parse(f, 'fasta'):
            # Exclude the Ghana project, because they're recombinant.
            assert reference.name != 'HIV1-CRF02_AG-GH-AB286855-seed'


def test_genotype(tmpdir, hcv_db):
    contigs_fasta = Path(tmpdir) / "contigs.fasta"
    contigs_fasta.write_text("""\
>foo
TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA
>bar
CATCACATAGGAGACAGGGCTCCAGGACTGCACCATGCTCGTGTGTGGCGACGAC
""")
    blast_csv = StringIO()
    expected_blast_csv = """\
contig_num,ref_name,score,match,pident,start,end,ref_start,ref_end
2,HCV-1g,37,0.67,100,19,55,8506,8542
2,HCV-1a,41,0.75,100,15,55,8518,8558
1,HCV-1a,41,1.0,100,1,41,8187,8227
"""

    genotype(contigs_fasta, blast_csv=blast_csv)

    assert expected_blast_csv == blast_csv.getvalue()


def test_fasta_to_csv_two_sequences(tmpdir, hcv_db):
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

    fasta_to_csv(contigs_fasta, contigs_csv)

    assert expected_contigs_csv == contigs_csv.getvalue()


def test_fasta_to_csv_two_groups(tmpdir, hcv_db):
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

    fasta_to_csv(contigs_fasta, contigs_csv)

    assert expected_contigs_csv == contigs_csv.getvalue()


def test_fasta_to_csv_not_found(tmpdir, hcv_db):
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

    fasta_to_csv(contigs_fasta, contigs_csv)

    assert expected_contigs_csv == contigs_csv.getvalue()


def test_fasta_to_csv_partial_match(tmpdir, hcv_db):
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

    fasta_to_csv(contigs_fasta, contigs_csv)

    assert expected_contigs_csv == contigs_csv.getvalue()


def test_fasta_to_csv_reversed_match(tmpdir, hcv_db):
    """ If BLAST match is reversed, then reverse the contig before reporting. """
    contigs_fasta = Path(tmpdir) / "contigs.fasta"
    contigs_fasta.write_text("""\
>foo
TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA
>bar
GTCGTCGCCACACACGAGCATGGTGCAGTCCTGGAGCCCTGTCTCCTATGTGATG
""")
    contigs_csv = StringIO()
    expected_contigs_csv = """\
ref,match,group_ref,contig
HCV-1a,1.0,HCV-1a,TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA
HCV-1a,0.75,HCV-1a,CATCACATAGGAGACAGGGCTCCAGGACTGCACCATGCTCGTGTGTGGCGACGAC
"""

    fasta_to_csv(contigs_fasta, contigs_csv)

    assert expected_contigs_csv == contigs_csv.getvalue()


def test_fasta_to_csv(tmpdir, hcv_db):
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
contig_num,ref_name,score,match,pident,start,end,ref_start,ref_end
2,HCV-1g,37,0.67,100,19,55,8506,8542
2,HCV-1a,41,0.75,100,15,55,8518,8558
1,HCV-1a,41,1.0,100,1,41,8187,8227
"""

    fasta_to_csv(contigs_fasta, contigs_csv, blast_csv=blast_csv)

    assert expected_contigs_csv == contigs_csv.getvalue()
    assert expected_blast_csv == blast_csv.getvalue()


def test_fasta_to_csv_none(tmpdir, hcv_db):
    contigs_fasta = Path(tmpdir) / 'contigs.fasta'
    assert not contigs_fasta.exists()

    contigs_csv = StringIO()
    expected_contigs_csv = """\
ref,match,group_ref,contig
"""

    fasta_to_csv(contigs_fasta, contigs_csv)

    assert expected_contigs_csv == contigs_csv.getvalue()
