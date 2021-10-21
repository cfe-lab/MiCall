from io import StringIO

import pytest

from micall.utils.contig_counts import ContigCounts, main


def test_contig_positions():
    genome_coverage_csv = StringIO("""\
contig,coordinates,query_nuc_pos,refseq_nuc_pos
1-HIV1-B-FR-K03455-seed,HIV1-B-FR-K03455-seed,1,99
1-HIV1-B-FR-K03455-seed,HIV1-B-FR-K03455-seed,2,100
1-HIV1-B-FR-K03455-seed,HIV1-B-FR-K03455-seed,3,101
1-HIV1-B-FR-K03455-seed,HIV1-B-FR-K03455-seed,4,102
1-HIV1-B-FR-K03455-seed,HIV1-B-FR-K03455-seed,5,103
""")
    counts = ContigCounts(start=100, end=102)
    expected_map = {'1-HIV1-B-FR-K03455-seed': {2: 100, 3: 101, 4: 102}}

    counts.read_genome_coverage(genome_coverage_csv)

    assert counts.position_map == expected_map


def test_contig_positions_insertion():
    genome_coverage_csv = StringIO("""\
contig,coordinates,query_nuc_pos,refseq_nuc_pos
1-HIV1-B-FR-K03455-seed,HIV1-B-FR-K03455-seed,1,99
1-HIV1-B-FR-K03455-seed,HIV1-B-FR-K03455-seed,2,100
1-HIV1-B-FR-K03455-seed,HIV1-B-FR-K03455-seed,3,101
1-HIV1-B-FR-K03455-seed,HIV1-B-FR-K03455-seed,4,
1-HIV1-B-FR-K03455-seed,HIV1-B-FR-K03455-seed,5,102
""")
    counts = ContigCounts(start=100, end=102)
    expected_map = {'1-HIV1-B-FR-K03455-seed': {2: 100, 3: 101, 5: 102}}

    counts.read_genome_coverage(genome_coverage_csv)

    assert counts.position_map == expected_map


def test_contig_positions_deletion():
    genome_coverage_csv = StringIO("""\
contig,coordinates,query_nuc_pos,refseq_nuc_pos
1-HIV1-B-FR-K03455-seed,HIV1-B-FR-K03455-seed,1,99
1-HIV1-B-FR-K03455-seed,HIV1-B-FR-K03455-seed,2,100
1-HIV1-B-FR-K03455-seed,HIV1-B-FR-K03455-seed,3,101
1-HIV1-B-FR-K03455-seed,HIV1-B-FR-K03455-seed,,102
1-HIV1-B-FR-K03455-seed,HIV1-B-FR-K03455-seed,4,103
""")
    counts = ContigCounts(start=100, end=103)
    expected_map = {'1-HIV1-B-FR-K03455-seed': {2: 100, 3: 101, 4: 103}}

    counts.read_genome_coverage(genome_coverage_csv)

    assert counts.position_map == expected_map


def test_reference_prefix_good():
    genome_coverage_csv = StringIO("""\
contig,coordinates,query_nuc_pos,refseq_nuc_pos
1-HIV1-B-FR-K03455-seed,HIV1-B-FR-K03455-seed,1,99
1-HIV1-B-FR-K03455-seed,HIV1-B-FR-K03455-seed,2,100
1-HIV1-B-FR-K03455-seed,HIV1-B-FR-K03455-seed,3,101
1-HIV1-B-FR-K03455-seed,HIV1-B-FR-K03455-seed,4,102
1-HIV1-B-FR-K03455-seed,HIV1-B-FR-K03455-seed,5,103
""")
    counts = ContigCounts(start=100, end=102)
    counts.reference_prefix = 'HIV'
    expected_map = {'1-HIV1-B-FR-K03455-seed': {2: 100, 3: 101, 4: 102}}

    counts.read_genome_coverage(genome_coverage_csv)

    assert counts.position_map == expected_map


def test_reference_prefix_bad():
    genome_coverage_csv = StringIO("""\
contig,coordinates,query_nuc_pos,refseq_nuc_pos
1-HIV1-B-FR-K03455-seed,HIV1-B-FR-K03455-seed,1,99
1-HIV1-B-FR-K03455-seed,HIV1-B-FR-K03455-seed,2,100
1-HIV1-B-FR-K03455-seed,HIV1-B-FR-K03455-seed,3,101
1-HIV1-B-FR-K03455-seed,HIV1-B-FR-K03455-seed,4,102
1-HIV1-B-FR-K03455-seed,HIV1-B-FR-K03455-seed,5,103
""")
    counts = ContigCounts(start=100, end=102)
    counts.reference_prefix = 'HCV'
    expected_map = {}

    counts.read_genome_coverage(genome_coverage_csv)

    assert counts.position_map == expected_map


def test_read_aligned():
    aligned_csv = StringIO("""\
refname,qcut,rank,count,offset,seq
1-HIV1-B-FR-K03455-seed,15,0,1,0,ACGT
""")
    counts = ContigCounts(start=100, end=102)
    counts.reference_prefix = 'HCV'
    counts.position_map = {'1-HIV1-B-FR-K03455-seed': {2: 100, 3: 101, 4: 102}}
    expected_counts = {'1-HIV1-B-FR-K03455-seed': {100: {'C': 1},
                                                   101: {'G': 1},
                                                   102: {'T': 1}}}

    counts.read_aligned(aligned_csv)

    assert counts.counts == expected_counts


def test_read_aligned_skip_middle():
    aligned_csv = StringIO("""\
refname,qcut,rank,count,offset,seq
1-HIV1-B-FR-K03455-seed,15,0,1,0,ACnT
""")
    counts = ContigCounts(start=100, end=102)
    counts.reference_prefix = 'HCV'
    counts.position_map = {'1-HIV1-B-FR-K03455-seed': {2: 100, 3: 101, 4: 102}}
    expected_counts = {'1-HIV1-B-FR-K03455-seed': {100: {'C': 1},
                                                   102: {'T': 1}}}

    counts.read_aligned(aligned_csv)

    assert counts.counts == expected_counts


def test_read_aligned_count():
    aligned_csv = StringIO("""\
refname,qcut,rank,count,offset,seq
1-HIV1-B-FR-K03455-seed,15,0,2,0,ACGT
1-HIV1-B-FR-K03455-seed,15,0,3,0,ACCT
""")
    counts = ContigCounts(start=100, end=102)
    counts.reference_prefix = 'HCV'
    counts.position_map = {'1-HIV1-B-FR-K03455-seed': {2: 100, 3: 101, 4: 102}}
    expected_counts = {'1-HIV1-B-FR-K03455-seed': {100: {'C': 5},
                                                   101: {'C': 3, 'G': 2},
                                                   102: {'T': 5}}}

    counts.read_aligned(aligned_csv)

    assert counts.counts == expected_counts


def test_read_aligned_offset():
    aligned_csv = StringIO("""\
refname,qcut,rank,count,offset,seq
1-HIV1-B-FR-K03455-seed,15,0,1,10,ACGT
""")
    counts = ContigCounts(start=100, end=102)
    counts.reference_prefix = 'HCV'
    counts.position_map = {'1-HIV1-B-FR-K03455-seed': {12: 100, 13: 101, 14: 102}}
    expected_counts = {'1-HIV1-B-FR-K03455-seed': {100: {'C': 1},
                                                   101: {'G': 1},
                                                   102: {'T': 1}}}

    counts.read_aligned(aligned_csv)

    assert counts.counts == expected_counts


def test_read_aligned_past_ends():
    aligned_csv = StringIO("""\
refname,qcut,rank,count,offset,seq
1-HIV1-B-FR-K03455-seed,15,0,1,8,GGAC
1-HIV1-B-FR-K03455-seed,15,0,2,12,GTGG
""")
    counts = ContigCounts(start=100, end=102)
    counts.reference_prefix = 'HCV'
    counts.position_map = {'1-HIV1-B-FR-K03455-seed': {12: 100, 13: 101, 14: 102}}
    expected_counts = {'1-HIV1-B-FR-K03455-seed': {100: {'C': 1},
                                                   101: {'G': 2},
                                                   102: {'T': 2}}}

    counts.read_aligned(aligned_csv)

    assert counts.counts == expected_counts


def test_display():
    counts = ContigCounts(start=100, end=102)
    counts.counts = {'1-HIV1-B-FR-K03455-seed': {100: {'C': 5000},
                                                 101: {'C': 2000, 'G': 3000},
                                                 103: {'T': 5000, 'A': 1}}}
    expected_display = """\
1-HIV1-B-FR-K03455-seed:
100: C (1.00)
101: G (0.60), C (0.40)
103: T (1.00)
"""

    display = counts.display()

    assert display == expected_display


def test_integration(tmp_path, capsys):
    coverage_path = tmp_path / 'genome_coverage.csv'
    coverage_path.write_text("""\
contig,coordinates,query_nuc_pos,refseq_nuc_pos
1-HIV1-B-FR-K03455-seed,HIV1-B-FR-K03455-seed,1,99
1-HIV1-B-FR-K03455-seed,HIV1-B-FR-K03455-seed,2,100
1-HIV1-B-FR-K03455-seed,HIV1-B-FR-K03455-seed,3,101
1-HIV1-B-FR-K03455-seed,HIV1-B-FR-K03455-seed,4,102
1-HIV1-B-FR-K03455-seed,HIV1-B-FR-K03455-seed,5,103
""")
    aligned_path = tmp_path / 'aligned.csv'
    aligned_path.write_text("""\
refname,qcut,rank,count,offset,seq
1-HIV1-B-FR-K03455-seed,15,0,1,0,ACGT
""")
    expected_display = """\
1-HIV1-B-FR-K03455-seed:
100: C (1.00)
101: G (1.00)
102: T (1.00)
"""

    main(['100', '102', str(tmp_path)])

    captured = capsys.readouterr()

    assert captured.out == expected_display


def test_integration_missing_file(tmp_path, capsys):
    coverage_path = tmp_path / 'genome_coverage.csv'
    coverage_path.write_text("junk")

    # Don't write aligned.csv
    aligned_path = tmp_path / 'aligned.csv'
    expected_error = f"error: argument scratch: can't open {str(aligned_path)!r}"

    with pytest.raises(SystemExit):
        main(['100', '102', str(tmp_path)])

    captured = capsys.readouterr()

    assert expected_error in captured.err
