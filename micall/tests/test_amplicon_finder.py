from collections import Counter
from io import StringIO
from pathlib import Path

from micall.core.amplicon_finder import merge_reads, count_kmers, calculate_entropy_from_counts, calculate_entropy, \
    merge_for_entropy, plot_merge_lengths

microtest_path = Path(__file__).parent / 'microtest'


def test_merge_reads(tmpdir):
    merge_lengths_csv = StringIO()
    expected_merge_lengths = """\
merge_length,count
99,1
117,2
139,1
154,1
175,1
182,2
191,1
194,1
195,1
196,1
213,1
224,1
226,1
232,2
244,1
249,1
256,1
261,1
265,1
359,1
391,1
420,1
431,1
446,1
450,1
465,1
"""
    merge_reads(microtest_path / '2160A-HCV_S19_L001_R1_001.fastq',
                microtest_path / '2160A-HCV_S19_L001_R2_001.fastq',
                tmpdir/'2160A-joined.fastq',
                merge_lengths_csv)

    assert merge_lengths_csv.getvalue() == expected_merge_lengths


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


def test_plot_merge_lengths():
    read_entropy = StringIO("""\
merge_length,count,entropy
100,5,0.2
101,5,0.2
102,5,0.9
103,5,0.2
""")

    figure = plot_merge_lengths(read_entropy)

    assert len(figure.axes) == 2
    assert figure.axes[0].get_title() == 'Count of merged reads at each length'
    assert figure.axes[1].get_title() == ''


def test_merge_for_entropy(tmpdir):
    read_entropy_csv = StringIO()
    expected_read_entropy = """\
merge_length,count,entropy
99,1,6.459431618637298
117,2,7.727920454563198
139,1,7.0
154,1,7.159871336778388
175,1,7.357552004618085
182,2,8.417852514885896
191,1,7.491853096329677
194,1,7.515699838284042
195,1,7.523561956057013
196,1,7.531381460516312
213,1,7.658211482751794
224,1,7.734709620225842
226,1,7.748192849589461
232,2,7.93269893948193
244,1,7.864186144654279
249,1,7.894817763307945
256,1,7.93663793900257
261,1,7.965784284662087
265,1,7.988684686772166
359,1,8.442943495848729
391,1,8.569855608330949
420,1,8.675957032941751
431,1,8.714245517666122
446,1,8.764871590736089
450,1,8.778077129535358
465,1,8.826548487290914
"""
    merge_for_entropy(microtest_path / '2160A-HCV_S19_L001_R1_001.fastq',
                      microtest_path / '2160A-HCV_S19_L001_R2_001.fastq',
                      read_entropy_csv,
                      Path(tmpdir))

    assert read_entropy_csv.getvalue() == expected_read_entropy
