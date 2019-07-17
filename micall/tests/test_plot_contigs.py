from csv import DictReader
from io import StringIO

from genetracks import Figure, Track, Multitrack, Label, Coverage

from micall.core.plot_contigs import build_contigs_figure, summarize_figure, build_coverage_figure

HCV_HEADER = ('C[342-915], E1[915-1491], E2[1491-2580], P7[2580-2769], '
              'NS2[2769-3420], NS3[3420-5313], NS4A[5313-5475], NS4B[5475-6258], '
              'NS5A[6258-7602], NS5B[7602-9378]')
HIV_HEADER = '''\
5' LTR[0-634], gag[790-2292], vif[5041-5619], tat[8379-8469], nef[8797-9417]
tat[5831-6045], vpu[6062-6310], rev[8379-8653], 3' LTR[9086-9719]
pol[2085-5096], vpr[5559-5850], rev[5970-6045], env[6225-8795]'''


def test_summarize_labels():
    figure = Figure()
    figure.add(Track(0, 0, label="Foo:"))
    figure.add(Track(0, 0, label="Bar:"))
    expected_summary = """\
Foo:
Bar:
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


def test_summarize_spans():
    figure = Figure()
    figure.add(Track(10, 20, label="Foo:", color=''))
    figure.add(Track(10, 20, label="Bar"))
    expected_summary = """\
Foo:
Bar[10-20]
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


def test_summarize_coverage():
    figure = Figure()
    figure.add(Coverage(10, 20, [10, 20, 1]), gap=-4)
    figure.add(Track(10, 20, label="Bar"))
    expected_summary = """\
Coverage 10, 20, 1
Bar[10-20]
"""

    summary = summarize_figure(figure)

    assert expected_summary == summary


def test_simple():
    blast_rows = DictReader(StringIO("""\
contig_name,ref_name,score,match,start,end,ref_start,ref_end
contig_1,HCV-1a,41,1.0,1,41,8187,8227
"""))
    expected_figure = f"""\
{HCV_HEADER}
HCV-1a:
contig_1[8187-8227]
"""

    figure = build_contigs_figure(blast_rows)

    assert expected_figure == summarize_figure(figure)


def test_multiple_organisms():
    blast_rows = DictReader(StringIO("""\
contig_name,ref_name,score,match,start,end,ref_start,ref_end
contig_1,HCV-1a,41,1.0,1,41,8187,8227
contig_2,HIV1-A1-CD-AM000053-seed,100,1.0,1,100,3000,3100
"""))
    expected_figure = f"""\
{HIV_HEADER}
HIV1-A1-CD-AM000053-seed:
contig_2[3000-3100]
{HCV_HEADER}
HCV-1a:
contig_1[8187-8227]
"""

    figure = build_contigs_figure(blast_rows)

    assert expected_figure == summarize_figure(figure)


def test_multiple_refs():
    blast_rows = DictReader(StringIO("""\
contig_name,ref_name,score,match,start,end,ref_start,ref_end
contig_1,HCV-1a,41,1.0,1,41,8187,8227
contig_2,HCV-2a,50,1.0,51,100,8300,8350
"""))
    expected_figure = f"""\
{HCV_HEADER}
HCV-2a:
contig_2[8300-8350]
HCV-1a:
contig_1[8187-8227]
"""

    figure = build_contigs_figure(blast_rows)

    assert expected_figure == summarize_figure(figure)


def test_organisms_sorted_by_score():
    blast_rows = DictReader(StringIO("""\
contig_name,ref_name,score,match,start,end,ref_start,ref_end
contig_1,HCV-1a,30,1.0,1,30,8001,8030
contig_2,HCV-2c,40,1.0,1,40,8101,8140
contig_3,HIV1-A1-CD-AM000053-seed,50,1.0,1,50,3001,3050
"""))
    expected_figure = f"""\
{HCV_HEADER}
HCV-2c:
contig_2[8101-8140]
HCV-1a:
contig_1[8001-8030]
{HIV_HEADER}
HIV1-A1-CD-AM000053-seed:
contig_3[3001-3050]
"""

    figure = build_contigs_figure(blast_rows)

    assert expected_figure == summarize_figure(figure)


def test_multiple_contigs():
    blast_rows = DictReader(StringIO("""\
contig_name,ref_name,score,match,start,end,ref_start,ref_end
contig_1,HCV-1a,41,1.0,1,41,8187,8227
contig_2,HCV-1a,50,1.0,51,100,8300,8350
"""))
    expected_figure = f"""\
{HCV_HEADER}
HCV-1a:
contig_1[8187-8227]
contig_2[8300-8350]
"""

    figure = build_contigs_figure(blast_rows)

    assert expected_figure == summarize_figure(figure)


def test_contigs_limit():
    blast_rows = DictReader(StringIO("""\
contig_name,ref_name,score,match,start,end,ref_start,ref_end
contig_1,HCV-1a,10,1.0,1,10,8001,8010
contig_1,HCV-1a,20,1.0,101,120,8301,8320
contig_1,HCV-2c,50,1.0,1,50,8001,8050
contig_2,HIV1-A1-CD-AM000053-seed,15,1.0,1,15,3001,3015
contig_2,HCV-1a,60,1.0,101,160,8701,8760
"""))
    expected_figure = f"""\
{HCV_HEADER}
HCV-1a:
contig_1[8301-8320]
contig_2[8701-8760]
HCV-2c:
contig_1[8001-8050]
"""

    figure = build_contigs_figure(blast_rows, limit=3)

    assert expected_figure == summarize_figure(figure)


def test_no_contigs():
    blast_rows = DictReader(StringIO("""\
contig_name,ref_name,score,match,start,end,ref_start,ref_end
"""))
    expected_figure = """\
No contigs assembled.
"""

    figure = build_contigs_figure(blast_rows, limit=3)

    assert expected_figure == summarize_figure(figure)


def test_plot_contig_coverage():
    contig_coverage_csv = StringIO("""\
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
Coverage 5, 5, 7, 5, 5, 5
1-HCV-1a - depth 7[1-6]
"""

    figure = build_coverage_figure(contig_coverage_csv)

    assert expected_figure == summarize_figure(figure)


def test_plot_contig_coverage_offset():
    """ When a contig extends before the reference start, offset everything. """
    contig_coverage_csv = StringIO("""\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,ins,dels,coverage
1-HCV-1a,HCV-1a,1,-2,0,0,5
1-HCV-1a,HCV-1a,2,-1,0,0,5
1-HCV-1a,HCV-1a,3,0,0,0,7
1-HCV-1a,HCV-1a,4,1,0,0,5
1-HCV-1a,HCV-1a,5,2,0,0,5
1-HCV-1a,HCV-1a,6,3,0,0,5
""")
    expected_figure = """\
5'[4-344], C[345-917], E1[918-1493], E2[1494-2582], p7[2583-2771], \
NS2[2772-3422], NS3[3423-5315], NS4b[5478-6260], NS4a[5316-5477], \
NS5a[6261-7604], NS5b[7605-9377], 3'[9378-9649]
Coverage 5, 5, 7, 5, 5, 5
1-HCV-1a - depth 7[1-6]
"""

    figure = build_coverage_figure(contig_coverage_csv)

    assert expected_figure == summarize_figure(figure)


def test_plot_contig_coverage_partial():
    contig_coverage_csv = StringIO("""\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,ins,dels,coverage
1-HCV-1a-partial,,1,,0,0,5
1-HCV-1a-partial,,2,,0,0,5
1-HCV-1a-partial,,3,,0,0,7
1-HCV-1a-partial,,4,,0,0,5
1-HCV-1a-partial,,5,,0,0,5
1-HCV-1a-partial,,6,,0,0,5
""")
    expected_figure = """\
Partial Blast Results
Coverage 5, 5, 7, 5, 5, 5
1-HCV-1a-partial - depth 7[1-6]
"""

    figure = build_coverage_figure(contig_coverage_csv)

    assert expected_figure == summarize_figure(figure)


def test_plot_contig_coverage_empty():
    contig_coverage_csv = StringIO("""\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,ins,dels,coverage
""")
    expected_figure = """\
No contigs found.
"""

    figure = build_coverage_figure(contig_coverage_csv)

    assert expected_figure == summarize_figure(figure)
