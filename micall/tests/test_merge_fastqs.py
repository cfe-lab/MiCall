import pytest
import os
import sys
from io import StringIO
import csv
import gzip

from micall.core.merge_fastqs import main, merge_fastqs, BadTrimplanFieldNames, AmbiguousBadCycles

FASTQ_DATA_A_R1 = """@M01234:01:000000000-AAAAA:1:1101:2020:0001 2:N:0:1
AAGTCGACGAATG
+
ABCDEFGHIJKLM"""
FASTQ_DATA_A_R2 = """@M01234:01:000000000-AAAAA:1:1101:2020:0001 2:N:0:1
AACCTGAAGGGTT
+
ABCDEFGHIJKLM"""
FASTQ_DATA_B_R1 = """@M01234:02:000000000-AAAAA:1:1101:2020:0001 2:N:0:1
CCCTAAAGGTTTG
+
ABCDEFGHIJKLM"""
FASTQ_DATA_B_R2 = """@M01234:02:000000000-AAAAA:1:1101:2020:0001 2:N:0:1
ACGTACGTACGTA
+
ABCDEFGHIJKLM"""
FASTQ_DATA_C_R1 = """@M01234:03:000000000-AAAAA:1:1101:2020:0001 2:N:0:1
AAATAAAGGTTTG
+
ABMDEFGZIJKLM"""
FASTQ_DATA_C_R2 = """@M01234:03:000000000-AAAAA:1:1101:2020:0001 2:N:0:1
ABGAACGTACGTA
+
ABCAEFGHIUKLM"""


@pytest.fixture
def fastq_files(tmp_path):
    # Create forward and reverse fastq files for two samples: A and B.
    fastq_a_r1 = tmp_path / "fastq_a_r1.fastq"
    fastq_a_r2 = tmp_path / "fastq_a_r2.fastq"
    fastq_b_r1 = tmp_path / "fastq_b_r1.fastq"
    fastq_b_r2 = tmp_path / "fastq_b_r2.fastq"
    fastq_a_r1.write_text(FASTQ_DATA_A_R1)
    fastq_a_r2.write_text(FASTQ_DATA_A_R2)
    fastq_b_r1.write_text(FASTQ_DATA_B_R1)
    fastq_b_r2.write_text(FASTQ_DATA_B_R2)

    # Paths for output files
    fastq_result_r1 = tmp_path / "fastq_result_r1.fastq"
    fastq_result_r2 = tmp_path / "fastq_result_r2.fastq"

    return (fastq_a_r1, fastq_a_r2, fastq_b_r1, fastq_b_r2, fastq_result_r1, fastq_result_r2)


def test_basic_merge_fastqs(fastq_files):
    (fastq_a_r1, fastq_a_r2, fastq_b_r1, fastq_b_r2, fastq_result_r1, fastq_result_r2) \
        = fastq_files

    trimplan = {(fastq_a_r1, fastq_a_r2),
                (fastq_b_r1, fastq_b_r2),
                }
    mergeplan = {fastq_result_r1: [fastq_a_r1, fastq_b_r1],
                 fastq_result_r2: [fastq_a_r2, fastq_b_r2],
                 }

    merge_fastqs(trimplan, mergeplan)

    # Add some assertions
    assert os.path.exists(fastq_result_r1)
    assert os.path.exists(fastq_result_r2)

    # New checks: Check if the merged files contain the correct sequences and qualities.

    with open(fastq_result_r1, 'r') as f:
        data = f.read()
        # should contain the sequences and qualities of the first reads of both samples
        assert FASTQ_DATA_A_R1 in data
        assert FASTQ_DATA_B_R1 in data
        assert len(data.strip()) == len(FASTQ_DATA_A_R1.strip()) + 1 + len(FASTQ_DATA_B_R1.strip())

    with open(fastq_result_r2, 'r') as f:
        data = f.read()
        # should contain the sequences and qualities of the second reads of both samples
        assert FASTQ_DATA_A_R2 in data
        assert FASTQ_DATA_B_R2 in data
        assert len(data.strip()) == len(FASTQ_DATA_A_R2.strip()) + 1 + len(FASTQ_DATA_B_R2.strip())


@pytest.fixture
def csvs_with_fastq_files(fastq_files, tmp_path):
    (fastq_a_r1, fastq_a_r2, fastq_b_r1, fastq_b_r2, fastq_result_r1, fastq_result_r2) \
        = fastq_files

    trimplan = tmp_path / "trimplan.csv"
    mergeplan = tmp_path / "mergeplan.csv"

    trimplan.write_text(f"""
r1,r2
{fastq_a_r1},{fastq_a_r2}
{fastq_b_r1},{fastq_b_r2}
    """.strip())

    mergeplan.write_text(f"""
input,output
{fastq_a_r1},{fastq_result_r1}
{fastq_b_r1},{fastq_result_r1}
{fastq_a_r2},{fastq_result_r2}
{fastq_b_r2},{fastq_result_r2}
    """.strip())

    return (str(trimplan), str(mergeplan))


def test_basic_main(csvs_with_fastq_files, fastq_files):
    (fastq_a_r1, fastq_a_r2, fastq_b_r1, fastq_b_r2, fastq_result_r1, fastq_result_r2) \
        = fastq_files
    (trimplan, mergeplan) = csvs_with_fastq_files

    main([trimplan, mergeplan])

    # Add some assertions
    assert os.path.exists(fastq_result_r1)
    assert os.path.exists(fastq_result_r2)

    # New checks: Check if the merged files contain the correct sequences and qualities.

    with open(fastq_result_r1, 'r') as f:
        data = f.read()
        # should contain the sequences and qualities of the first reads of both samples
        assert FASTQ_DATA_A_R1 in data
        assert FASTQ_DATA_B_R1 in data
        assert len(data.strip()) == len(FASTQ_DATA_A_R1.strip()) + 1 + len(FASTQ_DATA_B_R1.strip())

    with open(fastq_result_r2, 'r') as f:
        data = f.read()
        # should contain the sequences and qualities of the second reads of both samples
        assert FASTQ_DATA_A_R2 in data
        assert FASTQ_DATA_B_R2 in data
        assert len(data.strip()) == len(FASTQ_DATA_A_R2.strip()) + 1 + len(FASTQ_DATA_B_R2.strip())


@pytest.fixture
def bad_cycles_csv(tmp_path):
    data = """tile,cycle,errorrate
1101,1,1.0
2,2,7.5
"""
    bad_cycles_path = tmp_path / 'bad_cycles.csv'
    with open(bad_cycles_path, 'w') as f:
        f.write(data)
    return str(bad_cycles_path)


@pytest.fixture
def fastq_files_with_bad_cycles(bad_cycles_csv, tmp_path):
    fastq_a_r1 = tmp_path / "fastq_a_r1.fastq"
    fastq_a_r2 = tmp_path / "fastq_a_r2.fastq"
    fastq_b_r1 = tmp_path / "fastq_b_r1.fastq"
    fastq_b_r2 = tmp_path / "fastq_b_r2.fastq"
    fastq_a_r1.write_text(FASTQ_DATA_A_R1)
    fastq_a_r2.write_text(FASTQ_DATA_A_R2)
    fastq_b_r1.write_text(FASTQ_DATA_B_R1)
    fastq_b_r2.write_text(FASTQ_DATA_B_R2)

    fastq_result_r1 = tmp_path / "fastq_result_r1.fastq"
    fastq_result_r2 = tmp_path / "fastq_result_r2.fastq"

    return (fastq_a_r1, fastq_a_r2, fastq_b_r1, fastq_b_r2,
            fastq_result_r1, fastq_result_r2, bad_cycles_csv)


@pytest.fixture
def csvs_with_fastq_files_with_bad_cycles(fastq_files_with_bad_cycles, bad_cycles_csv, tmp_path):
    (fastq_a_r1, fastq_a_r2, fastq_b_r1, fastq_b_r2, fastq_result_r1, fastq_result_r2, bad_cycles_csv) \
        = fastq_files_with_bad_cycles

    trimplan = tmp_path / "trimplan.csv"
    mergeplan = tmp_path / "mergeplan.csv"

    trimplan.write_text(f"""
r1,r2,bad_cycles
{fastq_a_r1},{fastq_a_r2},{bad_cycles_csv}
{fastq_b_r1},{fastq_b_r2},{''}
    """.strip())

    mergeplan.write_text(f"""
input,output
{fastq_a_r1},{fastq_result_r1}
{fastq_b_r1},{fastq_result_r1}
{fastq_a_r2},{fastq_result_r2}
{fastq_b_r2},{fastq_result_r2}
    """.strip())

    return (str(trimplan), str(mergeplan))


def test_rejected_reads(csvs_with_fastq_files_with_bad_cycles, fastq_files_with_bad_cycles):
    (fastq_a_r1, fastq_a_r2, fastq_b_r1, fastq_b_r2, fastq_result_r1, fastq_result_r2, bad_cycles_csv) \
        = fastq_files_with_bad_cycles
    (trimplan, mergeplan) = csvs_with_fastq_files_with_bad_cycles

    main([trimplan, mergeplan])

    with open(fastq_result_r1, 'r') as f:
        data = f.read()
        assert FASTQ_DATA_A_R1 not in data # first nucleotide filtered out.
        assert FASTQ_DATA_B_R1 in data
        assert len(data.strip()) == len(FASTQ_DATA_A_R1.strip()) + 1 + len(FASTQ_DATA_B_R1.strip())
    with open(fastq_result_r2, 'r') as f:
        data = f.read()
        assert FASTQ_DATA_A_R2 in data
        assert FASTQ_DATA_B_R2 in data
        assert len(data.strip()) == len(FASTQ_DATA_A_R2.strip()) + 1 + len(FASTQ_DATA_B_R2.strip())


# Helper function to write data to a gzipped file
def write_to_gzip(filepath, data):
    with gzip.open(filepath, 'wt') as f:
        f.write(data)


@pytest.fixture
def gzipped_fastq_files(tmp_path):
    # Create gzipped forward and reverse fastq files for two samples: A and B.
    fastq_a_r1 = tmp_path / "fastq_a_r1.fastq.gz"
    fastq_a_r2 = tmp_path / "fastq_a_r2.fastq.gz"
    fastq_b_r1 = tmp_path / "fastq_b_r1.fastq"
    fastq_b_r2 = tmp_path / "fastq_b_r2.fastq.gz"
    write_to_gzip(fastq_a_r1, FASTQ_DATA_A_R1)
    write_to_gzip(fastq_a_r2, FASTQ_DATA_A_R2)
    fastq_b_r1.write_text(FASTQ_DATA_B_R1)
    write_to_gzip(fastq_b_r2, FASTQ_DATA_B_R2)

    # Paths for output files
    fastq_result_r1 = tmp_path / "fastq_result_r1.fastq.gz"
    fastq_result_r2 = tmp_path / "fastq_result_r2.fastq"

    return (fastq_a_r1, fastq_a_r2, fastq_b_r1, fastq_b_r2, fastq_result_r1, fastq_result_r2)


@pytest.fixture
def csvs_with_gzipped_fastq_files(gzipped_fastq_files, tmp_path):
    (fastq_a_r1, fastq_a_r2, fastq_b_r1, fastq_b_r2, fastq_result_r1, fastq_result_r2) \
        = gzipped_fastq_files

    trimplan = tmp_path / "trimplan.csv"
    mergeplan = tmp_path / "mergeplan.csv"
    zip_file = tmp_path / "zipfile.csv"

    trimplan.write_text(f"""
r1,r2
{fastq_a_r1},{fastq_a_r2}
{fastq_b_r1},{fastq_b_r2}
    """.strip())

    mergeplan.write_text(f"""
input,output
{fastq_a_r1},{fastq_result_r1}
{fastq_b_r1},{fastq_result_r1}
{fastq_a_r2},{fastq_result_r2}
{fastq_b_r2},{fastq_result_r2}
    """.strip())

    zip_file.write_text(f"""
file
{fastq_a_r1}
{fastq_a_r2}
{fastq_b_r2}
{fastq_result_r1}
    """.strip())

    return (str(trimplan), str(mergeplan), str(zip_file))


def test_gzipped_main(csvs_with_gzipped_fastq_files, gzipped_fastq_files):
    (fastq_a_r1, fastq_a_r2, fastq_b_r1, fastq_b_r2, fastq_result_r1, fastq_result_r2) = gzipped_fastq_files
    (trimplan, mergeplan, zip_file) = csvs_with_gzipped_fastq_files

    main([trimplan, mergeplan, "--zipfile", zip_file])

    # Add some assertions
    assert os.path.exists(fastq_result_r1)  # assert gzipped fastq_result_r1 exists
    assert os.path.exists(fastq_result_r2)  # assert gzipped fastq_result_r2 exists

    # Check if the merged gzipped files contain the correct sequences and qualities.
    with gzip.open(fastq_result_r1, 'rt') as f:
        data = f.read()
        assert FASTQ_DATA_A_R1 in data
        assert FASTQ_DATA_B_R1 in data
        assert len(data.strip()) == len(FASTQ_DATA_A_R1.strip()) + len(FASTQ_DATA_B_R1.strip()) + 1

    with open(fastq_result_r2, 'r') as f:
        data = f.read()
        assert FASTQ_DATA_A_R2 in data
        assert FASTQ_DATA_B_R2 in data
        assert len(data.strip()) == len(FASTQ_DATA_A_R2.strip()) + len(FASTQ_DATA_B_R2.strip()) + 1


@pytest.fixture
def many_fastq_files(tmp_path):
    # Create forward and reverse fastq files for two samples: A and B.
    fastq_a_r1 = tmp_path / "fastq_a_r1.fastq"
    fastq_a_r2 = tmp_path / "fastq_a_r2.fastq"
    fastq_b_r1 = tmp_path / "fastq_b_r1.fastq"
    fastq_b_r2 = tmp_path / "fastq_b_r2.fastq"
    fastq_c_r1 = tmp_path / "fastq_c_r1.fastq"
    fastq_c_r2 = tmp_path / "fastq_c_r2.fastq"
    fastq_a_r1.write_text(FASTQ_DATA_A_R1)
    fastq_a_r2.write_text(FASTQ_DATA_A_R2)
    fastq_b_r1.write_text(FASTQ_DATA_B_R1)
    fastq_b_r2.write_text(FASTQ_DATA_B_R2)
    fastq_c_r1.write_text(FASTQ_DATA_C_R1)
    fastq_c_r2.write_text(FASTQ_DATA_C_R2)

    # Paths for output files
    fastq_result_r1 = tmp_path / "fastq_result_r1.fastq"
    fastq_result_r2 = tmp_path / "fastq_result_r2.fastq"

    return (fastq_a_r1, fastq_a_r2, fastq_b_r1, fastq_b_r2, fastq_c_r1, fastq_c_r2, fastq_result_r1, fastq_result_r2)


@pytest.fixture
def csvs_with_many_fastq_files(many_fastq_files, tmp_path):
    (fastq_a_r1, fastq_a_r2, fastq_b_r1, fastq_b_r2, fastq_c_r1, fastq_c_r2, fastq_result_r1, fastq_result_r2) \
        = many_fastq_files

    trimplan = tmp_path / "trimplan.csv"
    mergeplan = tmp_path / "mergeplan.csv"

    trimplan.write_text(f"""
r1,r2
{fastq_a_r1},{fastq_a_r2}
{fastq_b_r1},{fastq_b_r2}
{fastq_c_r1},{fastq_c_r2}
    """.strip())

    mergeplan.write_text(f"""
input,output
{fastq_a_r1},{fastq_result_r1}
{fastq_b_r1},{fastq_result_r1}
{fastq_c_r1},{fastq_result_r1}
{fastq_a_r2},{fastq_result_r2}
{fastq_b_r2},{fastq_result_r2}
{fastq_c_r2},{fastq_result_r2}
    """.strip())

    return (str(trimplan), str(mergeplan))


def test_many_main(csvs_with_many_fastq_files, many_fastq_files):
    (fastq_a_r1, fastq_a_r2, fastq_b_r1, fastq_b_r2, fastq_c_r1, fastq_c_r2, fastq_result_r1, fastq_result_r2) \
        = many_fastq_files
    (trimplan, mergeplan) = csvs_with_many_fastq_files

    main([trimplan, mergeplan])

    # Add some assertions
    assert os.path.exists(fastq_result_r1)
    assert os.path.exists(fastq_result_r2)

    # New checks: Check if the merged files contain the correct sequences and qualities.

    with open(fastq_result_r1, 'r') as f:
        data = f.read()
        # should contain the sequences and qualities of the first reads of both samples
        assert FASTQ_DATA_A_R1 in data
        assert FASTQ_DATA_B_R1 in data
        assert FASTQ_DATA_C_R1 in data
        assert len(data.strip()) == len(FASTQ_DATA_A_R1.strip()) + 1 + len(FASTQ_DATA_B_R1.strip()) + 1 + len(FASTQ_DATA_C_R1.strip())

    with open(fastq_result_r2, 'r') as f:
        data = f.read()
        # should contain the sequences and qualities of the second reads of both samples
        assert FASTQ_DATA_A_R2 in data
        assert FASTQ_DATA_B_R2 in data
        assert FASTQ_DATA_C_R2 in data
        assert len(data.strip()) == len(FASTQ_DATA_A_R2.strip()) + 1 + len(FASTQ_DATA_B_R2.strip()) + 1 + len(FASTQ_DATA_C_R2.strip())


def test_bad_trimplan_field_names(tmp_path):
    trimplan = tmp_path / "trimplan.csv"
    mergeplan = tmp_path / "mergeplan.csv"

    trimplan.write_text(f"""r1,r2,unknown_field""".strip())
    mergeplan.write_text(f"""input,output""")

    with pytest.raises(BadTrimplanFieldNames):
        main([str(trimplan), str(mergeplan)])


def test_ambiguous_bad_cycles():
    with pytest.raises(AmbiguousBadCycles):
        merge_fastqs(trimplan={("file1", "file2")},
                     mergeplan={},
                     zipped=set(),
                     bad_cycles={"file1": "cycles1", "file2": "cycles2"})
