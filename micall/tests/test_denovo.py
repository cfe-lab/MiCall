from io import StringIO
from pathlib import Path

from Bio import SeqIO
from pytest import fixture, mark

from micall.utils.fasta_to_csv import DEFAULT_DATABASE, genotype
from micall.core.denovo import write_contig_refs, denovo
from micall.blast_db.make_blast_db import make_blast_db, DEFAULT_PROJECTS


@fixture(scope='session', name='hcv_db')
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

    contigs_stitched_csv = StringIO()
    write_contig_refs(str(contigs_fasta), contigs_csv, contigs_stitched_csv)

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

    contigs_stitched_csv = StringIO()
    write_contig_refs(str(contigs_fasta), contigs_csv, contigs_stitched_csv)

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

    contigs_stitched_csv = StringIO()
    write_contig_refs(str(contigs_fasta), contigs_csv, contigs_stitched_csv)

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

    contigs_stitched_csv = StringIO()
    write_contig_refs(str(contigs_fasta), contigs_csv, contigs_stitched_csv)

    assert expected_contigs_csv == contigs_csv.getvalue()


def test_write_contig_refs_reversed_match(tmpdir, hcv_db):
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

    contigs_stitched_csv = StringIO()
    write_contig_refs(str(contigs_fasta), contigs_csv, contigs_stitched_csv)

    assert expected_contigs_csv == contigs_csv.getvalue()


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


def test_write_contig_refs(tmpdir, hcv_db):
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

    contigs_stitched_csv = StringIO()
    write_contig_refs(str(contigs_fasta), contigs_csv, contigs_stitched_csv, blast_csv=blast_csv)

    assert expected_contigs_csv == contigs_csv.getvalue()
    assert expected_blast_csv == blast_csv.getvalue()


def test_write_contig_refs_none(tmpdir, hcv_db):
    contigs_fasta = Path(tmpdir) / 'contigs.fasta'
    assert not contigs_fasta.exists()

    contigs_csv = StringIO()
    expected_contigs_csv = """\
ref,match,group_ref,contig
"""

    contigs_stitched_csv = StringIO()
    write_contig_refs(str(contigs_fasta), contigs_csv, contigs_stitched_csv)

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
        contigs_stitched_csv = StringIO()
        write_contig_refs(str(contigs_fasta),
                          contigs_csv, contigs_stitched_csv,
                          merged_contigs_csv=merged_contigs_csv)

    assert expected_contigs_csv == contigs_csv.getvalue()


@mark.iva()  # skip with -k-iva
def test_denovo_iva(tmpdir, hcv_db):
    microtest_path = Path(__file__).parent / 'microtest'
    contigs_csv = StringIO()
    expected_contigs_csv = """\
ref,match,group_ref,contig
HCV-2a,1.0,HCV-2a,TGAGGGCCAAAAAGGTAACTTTTGATAGGATGCAAGTGC\
TCGACGCTCATTACGACTCAGTCTTAAAGGACATCAAGCTAGCGGCCTCCAAGGTCTCCG\
CGAGGCTCCTCACCCTGGAGGAGGCATGCCAGCTAACTCCACCCCATTCTGCAAGATCCAAATATGGGTTTGGGGCTA\
AGGAGGTGCGCAGCTTGTCCGGGAGGGCCGTTAACCACATCAAGTCCGTGTGGAAGGACCTCCTGGAAGACTCACAAA\
CACCAATTCCCACAACCATCATGGCCAAAAATGAAGTGTTCTGCGTGGACCCCACCAAGGGGGGTAAGAAAGCAGCTC\
GCCTCATCGTTTACCCTGACCTCGGCGTCAGGGTCTGCGAGAAGATGGCCCTTTATGATGTCACACAAAAGCTTCCTC\
AGGCGGTGATGGGGGCTTCTTATGGATTCCAGTACTCCC
"""

    denovo(str(microtest_path / '2160A-HCV_S19_L001_R1_001.fastq'),
           str(microtest_path / '2160A-HCV_S19_L001_R2_001.fastq'),
           contigs_csv,
           None,
           tmpdir)

    assert contigs_csv.getvalue() == expected_contigs_csv
