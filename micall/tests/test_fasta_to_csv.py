from io import StringIO
from pathlib import Path

from Bio import SeqIO
import pytest

from micall.utils.fasta_to_csv import DEFAULT_DATABASE, genotype
from micall.blast_db.make_blast_db import make_blast_db, DEFAULT_PROJECTS


@pytest.fixture(scope='session', name='hcv_db')
def check_hcv_db():
    db_path = Path(DEFAULT_DATABASE)
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


def test_make_blast_db_excludes_hivgha(hcv_db):
    fasta_path = Path(DEFAULT_DATABASE)
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

    genotype(str(contigs_fasta), blast_csv=blast_csv)

    assert expected_blast_csv == blast_csv.getvalue()
