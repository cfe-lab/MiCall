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
157,1
168,1
185,1
219,1
227,1
243,1
255,1
264,1
297,2
332,2
349,1
358,2
378,1
379,1
382,1
387,1
392,1
402,1
404,1
421,1
423,1
426,1
430,1
440,1
455,1
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
157,1,7.1898245588800185
168,1,7.294620748891628
185,1,7.44294349584873
219,1,7.700439718141093
227,1,7.7548875021634665
243,1,7.857980995127571
255,1,7.9307373375628885
264,1,7.98299357469431
297,2,8.540990217897273
332,2,8.637955966873083
349,1,8.400879436282183
358,2,8.444555541339067
378,1,8.519636252843211
379,1,8.523561956057016
382,1,8.535275376620802
387,1,8.554588851677634
392,1,8.573647187493318
402,1,8.611024797307355
404,1,8.61838550225861
421,1,8.679480099505447
423,1,8.686500527183219
426,1,8.696967526234284
430,1,8.710806433699352
440,1,8.744833837499545
455,1,8.79441586635011
"""
    merge_for_entropy(microtest_path / '2160A-HCV_S19_L001_R1_001.fastq',
                      microtest_path / '2160A-HCV_S19_L001_R2_001.fastq',
                      read_entropy_csv,
                      Path(tmpdir))

    assert expected_read_entropy == read_entropy_csv.getvalue()
