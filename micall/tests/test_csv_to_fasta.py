import sys
import pytest
from pathlib import Path
from micall.utils.csv_to_fasta import cli, main


def test_normal_input(tmpdir):
    contigs_fasta = Path(tmpdir) / "contigs.fasta"
    contigs_csv = Path(tmpdir) / "contigs.csv"

    contigs_csv.write_text("""\
ref,match,group_ref,contig
HCV-1a,1.0,HCV-1a,TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA
HCV-1a,0.75,HCV-1a,CATCACATAGGAGACAGGGCTCCAGGACTGCACCATGCTCGTGTGTGGCGACGAC
""")

    expected_contigs_fasta = """\
>1
TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA
>2
CATCACATAGGAGACAGGGCTCCAGGACTGCACCATGCTCGTGTGTGGCGACGAC
"""

    assert main([str(contigs_csv), str(contigs_fasta)]) == 0
    assert expected_contigs_fasta == contigs_fasta.read_text()


def test_empty(tmpdir):
    contigs_fasta = Path(tmpdir) / "contigs.fasta"
    contigs_csv = Path(tmpdir) / "contigs.csv"

    contigs_csv.write_text("""\
ref,match,group_ref,contig
""")

    expected_contigs_fasta = ""
    assert main([str(contigs_csv), str(contigs_fasta)]) == 0
    assert expected_contigs_fasta == contigs_fasta.read_text()


def test_main_invocation(tmpdir, monkeypatch):
    contigs_fasta = Path(tmpdir) / "contigs.fasta"
    contigs_csv = Path(tmpdir) / "contigs.csv"

    contigs_csv.write_text("""\
ref,match,group_ref,contig
HCV-1a,1.0,HCV-1a,TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA
HCV-1a,0.75,HCV-1a,CATCACATAGGAGACAGGGCTCCAGGACTGCACCATGCTCGTGTGTGGCGACGAC
""")

    expected_contigs_fasta = """\
>1
TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA
>2
CATCACATAGGAGACAGGGCTCCAGGACTGCACCATGCTCGTGTGTGGCGACGAC
"""

    argv = ['program', str(contigs_csv), str(contigs_fasta)]
    monkeypatch.setattr(sys, 'argv', argv)

    with pytest.raises(SystemExit):
        cli()

    assert expected_contigs_fasta == contigs_fasta.read_text()
