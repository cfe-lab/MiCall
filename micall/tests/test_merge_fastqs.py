import pytest
import os
import sys
from io import StringIO
import csv

from micall.core.merge_fastqs import main

FASTQ_DATA_A = """@M01234:01:000000000-AAAAA:1:1101:2020:0001 2:N:0:1
AAGTCGACGAATG
+
ABCDEFGHIJKLM"""
FASTQ_DATA_A2 = """@M01234:01:000000000-AAAAA:1:1101:2020:0001 2:N:0:1
AACCTGAAGGGTT
+
ABCDEFGHIJKLM"""
FASTQ_DATA_B = """@M01234:02:000000000-AAAAA:1:1101:2020:0001 2:N:0:1
CCCTAAAGGTTTG
+
ABCDEFGHIJKLM"""
FASTQ_DATA_B2 = """@M01234:02:000000000-AAAAA:1:1101:2020:0001 2:N:0:1
ACGTACGTACGTA
+
ABCDEFGHIJKLM"""


@pytest.fixture
def fastq_files(tmp_path):
    # Create forward and reverse fastq files for two samples: A and B.
    fastq1_a = tmp_path / "fastq1_a.fastq"
    fastq2_a = tmp_path / "fastq2_a.fastq"
    fastq1_b = tmp_path / "fastq1_b.fastq"
    fastq2_b = tmp_path / "fastq2_b.fastq"
    fastq1_a.write_text(FASTQ_DATA_A)
    fastq2_a.write_text(FASTQ_DATA_A2)
    fastq1_b.write_text(FASTQ_DATA_B)
    fastq2_b.write_text(FASTQ_DATA_B2)

    # Paths for output files
    fastq1_result = tmp_path / "fastq1_result.fastq"
    fastq2_result = tmp_path / "fastq2_result.fastq"

    return (fastq1_a, fastq2_a, fastq1_b, fastq2_b, fastq1_result, fastq2_result)


def test_basic_main(fastq_files):
    (fastq1_a, fastq2_a, fastq1_b, fastq2_b, fastq1_result, fastq2_result) = fastq_files
    args = [str(fastq1_a), str(fastq2_a), str(fastq1_b), str(fastq2_b), '--unzipped', str(fastq1_result), str(fastq2_result)]

    main(args)

    # Add some assertions
    assert os.path.exists(fastq_files[-2])  # assert fastq1_result exists
    assert os.path.exists(fastq_files[-1])  # assert fastq2_result exists

    # New checks: Check if the merged files contain the correct sequences and qualities.

    with open(fastq1_result, 'r') as f:
        data = f.read()
        # should contain the sequences and qualities of the first reads of both samples
        assert FASTQ_DATA_A in data
        assert FASTQ_DATA_B in data

    with open(fastq2_result, 'r') as f:
        data = f.read()
        # should contain the sequences and qualities of the second reads of both samples
        assert FASTQ_DATA_A2 in data
        assert FASTQ_DATA_B2 in data
