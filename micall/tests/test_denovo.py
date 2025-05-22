from io import StringIO
from pathlib import Path
import re

from pytest import mark

from micall.core.denovo import denovo
from micall.tests.test_fasta_to_csv import check_hcv_db, DEFAULT_DATABASE  # activates the fixture

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
    microtest_path = Path(__file__).parent / 'microtest'
    contigs_fasta = StringIO()
    expected_contigs_fasta = """\
>contig.00001
GTCATACTCCTGGACAGGGGCTCTGATAACTCCTTGTAGCCCCGAAGAGGAAAAGTTGCCGGTCAACCCCTTGAGCAACTCGTTGTTGCGATACCACAACAAGGTGTACAGTACCACATCAAAGACCGCCTCAC\
TGAGGGCCAAAAAGGTAACTTTTGATAGGATGCAAGTGC\
TCGACGCTCATTACGACTCAGTCTTAAAGGACATCAAGCTAGCGGCCTCCAAGGTCTCCG\
CGAGGCTCCTCACCCTGGAGGAGGCATGCCAGCTAACTCCACCCCATTCTGCAAGATCCAAATATGGGTTTGGGGCTA\
AGGAGGTGCGCAGCTTGTCCGGGAGGGCCGTTAACCACATCAAGTCCGTGTGGAAGGACCTCCTGGAAGACTCACAAA\
CACCAATTCCCACAACCATCATGGCCAAAAATGAAGTGTTCTGCGTGGACCCCACCAAGGGGGGTAAGAAAGCAGCTC\
GCCTCATCGTTTACCCTGACCTCGGCGTCAGGGTCTGCGAGAAGATGGCCCTTTATGATGTCACACAAAAGCTTCCTC\
AGGCGGTGATGGGGGCTTCTTATGGATTCCAGTACTCCC\
CCGCTCAGCGGGTGGAGTTTCTCTTGAAGGCATGGGCGGATAAGAGAGACCCTATGGGTTTTTCGTATGATACCCGATGCTTTGACTCAACCGTCAC
"""

    denovo(str(microtest_path / '2160A-HCV_S19_L001_R1_001.fastq'),
           str(microtest_path / '2160A-HCV_S19_L001_R2_001.fastq'),
           contigs_fasta,
           tmpdir)

    result = contigs_fasta.getvalue()
    expected = expected_contigs_fasta
    assert normalize_fasta(result) == normalize_fasta(expected)
