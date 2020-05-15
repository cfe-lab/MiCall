from collections import Counter
from io import StringIO
from pathlib import Path

from micall.core.amplicon_finder import merge_reads, count_kmers, calculate_entropy_from_counts, calculate_entropy, \
    merge_for_entropy

microtest_path = Path(__file__).parent / 'microtest'


def test_merge_reads(tmpdir):
    merge_lengths_csv = StringIO()
    expected_merge_lengths = """\
merge_length,count
114,1
123,1
127,1
134,1
153,1
160,1
166,1
167,1
170,1
183,1
190,1
199,1
215,1
258,1
291,1
298,1
304,1
311,1
318,1
322,1
331,1
342,1
369,1
380,1
394,1
398,1
401,2
412,1
428,1
"""
    merge_reads(microtest_path / '2160A-HCV_S19_L001_R1_001.fastq',
                microtest_path / '2160A-HCV_S19_L001_R2_001.fastq',
                tmpdir/'2160A-joined.fastq',
                merge_lengths_csv)

    assert expected_merge_lengths == merge_lengths_csv.getvalue()


def test_count_kmers():
    sequence = "CCTGGAAGACTCAC"
    counts = Counter()
    expected_counts = Counter({'CCTGGAAGACTC': 1,
                               'CTGGAAGACTCA': 1,
                               'TGGAAGACTCAC': 1})

    count_kmers(sequence, counts)

    assert expected_counts == counts


def test_count_kmers_multiple_sequences():
    sequence1 = "CCTGGAAGACTCAC"
    sequence2 = "TGGAAGACTCACA"
    counts = Counter()
    expected_counts = Counter({'CCTGGAAGACTC': 1,
                               'CTGGAAGACTCA': 1,
                               'TGGAAGACTCAC': 2,
                               'GGAAGACTCACA': 1})

    count_kmers(sequence1, counts)
    count_kmers(sequence2, counts)

    assert expected_counts == counts


def test_entropy_on_single_sequence():
    counts = Counter({'CCTGGAAGACTC': 4000})
    expected_entropy = 0.0

    entropy = calculate_entropy_from_counts(counts)

    assert expected_entropy == entropy


def test_entropy_on_multiple_sequences():
    counts = Counter({'CCTGGAAGACTC': 1000,
                      'CTGGAAGACTCA': 1000,
                      'TGGAAGACTCAC': 2000})
    expected_entropy = 1.5

    entropy = calculate_entropy_from_counts(counts)

    assert expected_entropy == entropy


def test_entropy_on_no_sequences():
    counts = Counter()
    expected_entropy = 0.0

    entropy = calculate_entropy_from_counts(counts)

    assert expected_entropy == entropy


def test_entropy_from_fastq():
    merge_lengths = {12: 3, 13: 2}
    merged_fastq = StringIO("""\
@M01234:01:000000000-AAAAA:1:1101:2160:0082 1:N:0:1
CTGGAAGACTCAC
+
AAAAAAAAAAAAA
@M01234:01:000000000-AAAAA:1:1101:2160:0087 1:N:0:1
TGGAAGACTCACA
+
AAAAAAAAAAAAA
@M01234:01:000000000-AAAAA:1:1101:2160:0089 1:N:0:1
AATTCCCACAAC
+
AAAAAAAAAAAA
@M01234:01:000000000-AAAAA:1:1101:2160:0090 1:N:0:1
AATTCCCACAAC
+
AAAAAAAAAAAA
@M01234:01:000000000-AAAAA:1:1101:2160:0093 1:N:0:1
AATTCCCACAAC
+
AAAAAAAAAAAA
""")
    # Two length 13 reads with shared overlap of 12 bases, should cause entropy
    # of 1.5 from 12-mer prevalences of 0.25, 0.5, and 0.25.
    # Three length 12 reads are identical, so entropy should be 0.
    expected_entropy = {13: 1.5, 12: 0.0}

    entropy = calculate_entropy(merged_fastq, merge_lengths)

    assert expected_entropy == entropy


def test_merge_for_entropy(tmpdir):
    read_entropy_csv = StringIO()
    expected_read_entropy = """\
merge_length,count,entropy
114,1,6.686500527183218
123,1,6.807354922057605
127,1,6.857980995127569
134,1,6.942514505339238
153,1,7.149747119504682
160,1,7.21916852046216
166,1,7.276124405274238
167,1,7.28540221886225
170,1,7.312882955284356
183,1,7.426264754702099
190,1,7.483815777264256
199,1,7.554588851677638
215,1,7.672425341971494
258,1,7.948367231584677
291,1,8.129283016944965
298,1,8.16490692667569
304,1,8.194756854422247
311,1,8.228818690495881
318,1,8.26209484537018
322,1,8.280770770130603
331,1,8.321928094887362
342,1,8.37068740680722
369,1,8.483815777264256
380,1,8.527477006060394
394,1,8.581200581924955
398,1,8.596189756144414
401,2,8.80733031374961
412,1,8.647458426454921
428,1,8.703903573444663
"""
    merge_for_entropy(microtest_path / '2160A-HCV_S19_L001_R1_001.fastq',
                      microtest_path / '2160A-HCV_S19_L001_R2_001.fastq',
                      read_entropy_csv,
                      Path(tmpdir))

    assert expected_read_entropy == read_entropy_csv.getvalue()
