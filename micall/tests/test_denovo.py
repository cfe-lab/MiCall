from pathlib import Path
import re

from pytest import mark

from micall.core.denovo import denovo
from micall.tests.test_fasta_to_csv import check_hcv_db, DEFAULT_DATABASE
from micall.utils.work_dir import WorkDir  # activates the fixture

# make linters not complain about unused imports.
assert check_hcv_db is not None
assert DEFAULT_DATABASE is not None


def normalize_fasta(content: str) -> str:
    result = re.sub(r'^>.*$', '>',
                    content,
                    flags=re.MULTILINE)
    result = ''.join(result.split('\n'))
    return result


@mark.iva()  # skip with -k-iva
def test_denovo_iva(tmpdir, hcv_db):
    tmpdir = Path(tmpdir)
    microtest_path = Path(__file__).parent / 'microtest'
    contigs_fasta = tmpdir / 'contigs.fasta'
    expected_contigs_fasta = """\
>contig.00001
TGAGGGCCAAAAAGGTAACTTTTGATAGGATGCAAGTGC\
TCGACGCTCATTACGACTCAGTCTTAAAGGACATCAAGCTAGCGGCCTCCAAGGTCTCCG\
CGAGGCTCCTCACCCTGGAGGAGGCATGCCAGCTAACTCCACCCCATTCTGCAAGATCCAAATATGGGTTTGGGGCTA\
AGGAGGTGCGCAGCTTGTCCGGGAGGGCCGTTAACCACATCAAGTCCGTGTGGAAGGACCTCCTGGAAGACTCACAAA\
CACCAATTCCCACAACCATCATGGCCAAAAATGAAGTGTTCTGCGTGGACCCCACCAAGGGGGGTAAGAAAGCAGCTC\
GCCTCATCGTTTACCCTGACCTCGGCGTCAGGGTCTGCGAGAAGATGGCCCTTTATGATGTCACACAAAAGCTTCCTC\
AGGCGGTGATGGGGGCTTCTTATGGATTCCAGTACTCCC
"""

    with WorkDir.using(tmpdir):
        denovo(microtest_path / '2160A-HCV_S19_L001_R1_001.fastq',
               microtest_path / '2160A-HCV_S19_L001_R2_001.fastq',
               contigs_fasta)

    result = contigs_fasta.read_text()
    expected = expected_contigs_fasta
    assert normalize_fasta(result) == normalize_fasta(expected)
