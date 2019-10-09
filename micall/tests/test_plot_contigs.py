from io import StringIO

import pytest
from genetracks import Figure, Track, Multitrack, Label, Coverage

from micall.core.plot_contigs import summarize_figure, \
    build_coverage_figure, SmoothCoverage, add_partial_banner

HCV_HEADER = ('C[342-915], E1[915-1491], E2[1491-2580], P7[2580-2769], '
              'NS2[2769-3420], NS3[3420-5313], NS4A[5313-5475], '
              'NS4B[5475-6258], NS5A[6258-7602], NS5B[7602-9378]')
HIV_HEADER = '''\
5' LTR[0-634], gag[790-2292], vif[5041-5619], tat[8379-8469], nef[8797-9417]
tat[5831-6045], vpu[6062-6310], rev[8379-8653], 3' LTR[9086-9719]
pol[2085-5096], vpr[5559-5850], rev[5970-6045], env[6225-8795]'''


def test_summarize_labels():
    figure = Figure()
    figure.add(Track(1, 200, label="Foo"))
    figure.add(Track(1, 200, label="Bar", color='none'))
    expected_summary = """\
Foo[1-200]
Bar(1-200)
"""

    summary = summarize_figure(figure)

    assert expected_summary == summary


def test_summarize_label_objects():
    figure = Figure()
    figure.add(Track(0, 0, label=Label(25, "Foo:")))
    figure.add(Track(0, 0, label="Bar:"))
    expected_summary = """\
Foo:
Bar:
"""

    summary = summarize_figure(figure)

    assert expected_summary == summary


def test_summarize_multitracks():
    figure = Figure()
    figure.add(Track(0, 0, label="Foo:"))
    figure.add(Multitrack([Track(10, 20, label="Bar"),
                           Track(30, 40, label="Baz")]))
    expected_summary = """\
Foo:
Bar[10-20], Baz[30-40]
"""

    summary = summarize_figure(figure)

    assert expected_summary == summary


def test_summarize_multitracks_with_separate_label():
    figure = Figure()
    figure.add(Track(0, 0, label="Foo:"))
    figure.add(Multitrack([Track(10, 20),
                           Track(30, 40),
                           Track(10, 40, label="Bar", color='none')]))
    expected_summary = """\
Foo:
[10-20], [30-40], Bar(10-40)
"""

    summary = summarize_figure(figure)

    assert expected_summary == summary


def test_summarize_regions():
    figure = Figure()
    figure.add(Track(1, 200, label="Foo", regions=[(50, 100, 'lightgreen'),
                                                   (110, 120, 'red')]))
    expected_summary = """\
Foo[1-200], lightgreen{50-100}, red{110-120}
"""

    summary = summarize_figure(figure)

    assert expected_summary == summary


def test_summarize_coverage():
    figure = Figure()
    figure.add(Coverage(10, 20, [11, 11, 21, 1, 1, 1]), gap=-4)
    figure.add(Track(12, 22, label="Bar"))
    expected_summary = """\
Coverage 11, 11, 21, 1, 1, 1
Bar[12-22]
"""

    summary = summarize_figure(figure)

    assert expected_summary == summary


def test_summarize_zero_coverage():
    figure = Figure()
    figure.add(Coverage(10, 20, [0, 0, 0]), gap=-4)
    figure.add(Track(10, 20, label="Bar"))

    with pytest.raises(ZeroDivisionError):
        summarize_figure(figure)


def test_summarize_smooth_coverage():
    figure = Figure()
    figure.add(SmoothCoverage(10, 20, [11, 11, 21, 1, 1, 1]), gap=-4)
    figure.add(Track(12, 22, label="Bar"))
    expected_summary = """\
Coverage 11x2, 21, 1x3
Bar[12-22]
"""

    summary = summarize_figure(figure)

    assert expected_summary == summary


def test_summarize_smooth_coverage_ten_percent():
    figure = Figure()
    figure.add(SmoothCoverage(10, 20, [100, 110, 111, 50]), gap=-4)
    figure.add(Track(12, 22, label="Bar"))
    expected_summary = """\
Coverage 100x2, 111, 50
Bar[12-22]
"""

    summary = summarize_figure(figure)

    assert expected_summary == summary


def test_plot_genome_coverage():
    genome_coverage_csv = StringIO("""\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,ins,dels,coverage
1-HCV-1a,HCV-1a,1,1,0,0,5
1-HCV-1a,HCV-1a,2,2,0,0,5
1-HCV-1a,HCV-1a,3,3,0,0,7
1-HCV-1a,HCV-1a,4,4,0,0,5
1-HCV-1a,HCV-1a,5,5,0,0,5
1-HCV-1a,HCV-1a,6,6,0,0,5
""")
    expected_figure = """\
5'[1-341], C[342-914], E1[915-1490], E2[1491-2579], p7[2580-2768], \
NS2[2769-3419], NS3[3420-5312], NS4b[5475-6257], NS4a[5313-5474], \
NS5a[6258-7601], NS5b[7602-9374], 3'[9375-9646]
Coverage 5x2, 7, 5x3
[1-6], 1-HCV-1a - depth 7(1-9646)
"""

    figure = build_coverage_figure(genome_coverage_csv)

    assert expected_figure == summarize_figure(figure)


def test_plot_genome_coverage_offset():
    """ When a contig extends before the reference start, offset everything. """
    genome_coverage_csv = StringIO("""\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,ins,dels,coverage
1-HCV-1a,HCV-1a,1,-2,0,0,5
1-HCV-1a,HCV-1a,2,-1,0,0,5
1-HCV-1a,HCV-1a,3,0,0,0,7
1-HCV-1a,HCV-1a,4,1,0,0,5
1-HCV-1a,HCV-1a,5,2,0,0,5
1-HCV-1a,HCV-1a,6,3,0,0,5
1-HCV-1a,HCV-1a,7,,0,0,5
1-HCV-1a,HCV-1a,8,,0,0,5
1-HCV-1a,HCV-1a,9,,0,0,5
1-HCV-1a,HCV-1a,10,4,0,0,5
1-HCV-1a,HCV-1a,11,5,0,0,5
1-HCV-1a,HCV-1a,12,6,0,0,5
2-unknown-partial,,1,,0,0,6
2-unknown-partial,,2,,0,0,6
2-unknown-partial,,3,,0,0,6
""")
    expected_figure = """\
5'[4-344], C[345-917], E1[918-1493], E2[1494-2582], p7[2583-2771], \
NS2[2772-3422], NS3[3423-5315], NS4b[5478-6260], NS4a[5316-5477], \
NS5a[6261-7604], NS5b[7605-9377], 3'[9378-9649]
Coverage 5x2, 7, 5x6
[1-9], 1-HCV-1a - depth 7(1-9649), lightgreen{7-9}
[4-503], [1004-1503], [2004-2503], [3004-3503], [4004-4503], \
[5004-5503], [6004-6503], [7004-7503], [8004-8503], [9004-9503], \
Partial Blast Results(4-9649)
Coverage 6x3
[4-6], 2-unknown-partial - depth 6(1-9649)
"""

    figure = build_coverage_figure(genome_coverage_csv)

    assert expected_figure == summarize_figure(figure)


def test_plot_genome_coverage_partial():
    genome_coverage_csv = StringIO("""\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,ins,dels,coverage
1-HCV-1a-partial,,1,,0,0,5
1-HCV-1a-partial,,2,,0,0,5
1-HCV-1a-partial,,3,,0,0,7
1-HCV-1a-partial,,4,,0,0,5
1-HCV-1a-partial,,5,,0,0,5
1-HCV-1a-partial,,6,,0,0,5
""")
    expected_figure = """\
[1-500], Partial Blast Results(1-500)
Coverage 5x2, 7, 5x3
[1-6], 1-HCV-1a-partial - depth 7(1-500)
"""

    figure = build_coverage_figure(genome_coverage_csv)

    assert expected_figure == summarize_figure(figure)


def test_plot_genome_coverage_partial_header():
    """ Last dash in the header banner can be less than 500 wide. """
    genome_coverage_csv = StringIO("""\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,ins,dels,coverage
""")
    genome_coverage_csv.seek(0, 2)  # EOF
    for i in range(1010):
        genome_coverage_csv.write(f'1-HCV-1a-partial,,{i+1},,0,0,5\n')
    genome_coverage_csv.seek(0)
    expected_figure = """\
[1-500], [1001-1010], Partial Blast Results(1-1010)
Coverage 5x1010
[1-1010], 1-HCV-1a-partial - depth 5(1-1010)
"""

    figure = build_coverage_figure(genome_coverage_csv)

    assert expected_figure == summarize_figure(figure)


def test_add_partial_banner():
    """ Last dash in the header banner can be less than 500 wide. """
    figure = Figure()
    add_partial_banner(figure, 0, 500)
    add_partial_banner(figure, 0, 700)
    add_partial_banner(figure, 0, 1200)

    expected_figure = """\
[1-500], Partial Blast Results(1-500)
[1-500], Partial Blast Results(1-700)
[1-500], [1001-1200], Partial Blast Results(1-1200)
"""

    assert expected_figure == summarize_figure(figure)


def test_plot_genome_coverage_sorted():
    genome_coverage_csv = StringIO("""\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage
1_3-HCV-1a,HCV-1a,1,1,0,5
1_3-HCV-1a,HCV-1a,2,2,0,5
1_3-HCV-1a,HCV-1a,3,3,0,7
1_3-HCV-1a,HCV-1a,4,4,0,5
1_3-HCV-1a,HCV-1a,5,5,0,5
1_3-HCV-1a,HCV-1a,6,6,0,5
1_3-HCV-1a,HCV-1a,101,101,0,15
1_3-HCV-1a,HCV-1a,102,102,0,15
1_3-HCV-1a,HCV-1a,103,103,0,17
1_3-HCV-1a,HCV-1a,104,104,0,15
1_3-HCV-1a,HCV-1a,105,105,0,15
1_3-HCV-1a,HCV-1a,106,106,0,15
contig-1-HCV-1a,HCV-1a,1,1,,
contig-1-HCV-1a,HCV-1a,2,2,,
contig-1-HCV-1a,HCV-1a,3,3,,
contig-1-HCV-1a,HCV-1a,4,4,,
contig-1-HCV-1a,HCV-1a,5,5,,
contig-1-HCV-1a,HCV-1a,6,6,,
contig-3-HCV-1a,HCV-1a,1,101,,
contig-3-HCV-1a,HCV-1a,2,102,,
contig-3-HCV-1a,HCV-1a,3,103,,
contig-3-HCV-1a,HCV-1a,4,104,,
contig-3-HCV-1a,HCV-1a,5,105,,
contig-3-HCV-1a,HCV-1a,6,106,,
2-HCV-1b-partial,,1,,0,5
2-HCV-1b-partial,,2,,0,5
2-HCV-1b-partial,,3,,0,29
2-HCV-1b-partial,,4,,0,5
2-HCV-1b-partial,,5,,0,5
2-HCV-1b-partial,,6,,0,5
""")
    expected_figure = """\
5'[1-341], C[342-914], E1[915-1490], E2[1491-2579], p7[2580-2768], \
NS2[2769-3419], NS3[3420-5312], NS4b[5475-6257], NS4a[5313-5474], \
NS5a[6258-7601], NS5b[7602-9374], 3'[9375-9646]
Coverage 5x2, 7, 5x3, 0x94, 15x2, 17, 15x3
[1-6], [101-106], 1_3-HCV-1a - depth 17(1-9646)
[1-6], contig-1-HCV-1a(1-9646)
[101-106], contig-3-HCV-1a(1-9646)
[1-500], [1001-1500], [2001-2500], [3001-3500], [4001-4500], \
[5001-5500], [6001-6500], [7001-7500], [8001-8500], [9001-9500], \
Partial Blast Results(1-9646)
Coverage 5x2, 29, 5x3
[1-6], 2-HCV-1b-partial - depth 29(1-9646)
"""

    figure = build_coverage_figure(genome_coverage_csv)

    assert expected_figure == summarize_figure(figure)


def test_plot_genome_coverage_zero():
    genome_coverage_csv = StringIO("""\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage
1-HCV-1b-partial,,1,,0,0
1-HCV-1b-partial,,2,,0,0
1-HCV-1b-partial,,3,,0,0
1-HCV-1b-partial,,4,,0,0
1-HCV-1b-partial,,5,,0,0
1-HCV-1b-partial,,6,,0,0
""")
    expected_figure = """\
[1-500], Partial Blast Results(1-500)
[1-6], 1-HCV-1b-partial(1-500)
"""

    figure = build_coverage_figure(genome_coverage_csv)

    assert expected_figure == summarize_figure(figure)


def test_plot_genome_coverage_empty():
    genome_coverage_csv = StringIO("""\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,ins,dels,coverage
""")
    expected_figure = """\
No contigs found.(1-500)
"""

    figure = build_coverage_figure(genome_coverage_csv)

    assert expected_figure == summarize_figure(figure)


def test_plot_genome_coverage_gap():
    genome_coverage_csv = StringIO("""\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,ins,dels,coverage
1-HCV-1a,HCV-1a,1,1,0,0,5
1-HCV-1a,HCV-1a,2,2,0,0,5
1-HCV-1a,HCV-1a,4,4,0,0,6
1-HCV-1a,HCV-1a,5,5,0,0,6
1-HCV-1a,HCV-1a,6,6,0,0,6
contig-1-HCV-1a,HCV-1a,1,1,0,,
contig-1-HCV-1a,HCV-1a,2,2,0,,
contig-1-HCV-1a,HCV-1a,4,4,0,,
contig-1-HCV-1a,HCV-1a,5,5,0,,
contig-1-HCV-1a,HCV-1a,6,6,0,,
""")
    expected_figure = """\
5'[1-341], C[342-914], E1[915-1490], E2[1491-2579], p7[2580-2768], \
NS2[2769-3419], NS3[3420-5312], NS4b[5475-6257], NS4a[5313-5474], \
NS5a[6258-7601], NS5b[7602-9374], 3'[9375-9646]
Coverage 5x2, 0, 6x3
[1-2], [4-6], 1-HCV-1a - depth 6(1-9646)
[1-2], [4-6], contig-1-HCV-1a(1-9646)
"""

    figure = build_coverage_figure(genome_coverage_csv)

    assert expected_figure == summarize_figure(figure)


def test_plot_genome_coverage_insertion():
    genome_coverage_csv = StringIO("""\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,ins,dels,coverage
1-HCV-1a,HCV-1a,1,1,0,0,5
1-HCV-1a,HCV-1a,2,2,0,0,5
1-HCV-1a,HCV-1a,3,3,0,0,5
1-HCV-1a,HCV-1a,4,,0,0,6
1-HCV-1a,HCV-1a,5,,0,0,6
1-HCV-1a,HCV-1a,6,,0,0,6
1-HCV-1a,HCV-1a,7,4,0,0,7
1-HCV-1a,HCV-1a,8,5,0,0,7
1-HCV-1a,HCV-1a,9,6,0,0,7
1-HCV-1a,HCV-1a,10,7,0,0,8
1-HCV-1a,HCV-1a,11,8,0,0,8
1-HCV-1a,HCV-1a,12,9,0,0,8
""")
    expected_figure = """\
5'[1-341], C[342-914], E1[915-1490], E2[1491-2579], p7[2580-2768], \
NS2[2769-3419], NS3[3420-5312], NS4b[5475-6257], NS4a[5313-5474], \
NS5a[6258-7601], NS5b[7602-9374], 3'[9375-9646]
Coverage 5x3, 7x3, 8x3
[1-9], 1-HCV-1a - depth 8(1-9646), lightgreen{4-6}
"""

    figure = build_coverage_figure(genome_coverage_csv)

    assert expected_figure == summarize_figure(figure)
