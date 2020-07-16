from base64 import standard_b64encode
from io import StringIO, BytesIO
from pathlib import Path
from turtle import Turtle

# noinspection PyPep8Naming
import drawSvg as draw
import pytest
from PIL import Image
from drawSvg import Drawing, Line, Lines, Circle, Text, Rectangle
from genetracks import Figure, Track, Multitrack, Label, Coverage

from micall.core.plot_contigs import summarize_figure, \
    build_coverage_figure, SmoothCoverage, add_partial_banner, Arrow, ArrowGroup


class SvgDiffer:
    def __init__(self):
        self.work_dir: Path = Path(__file__).parent / 'svg_diffs'
        self.work_dir.mkdir(exist_ok=True)
        self.mismatch_found = False
        for work_file in self.work_dir.iterdir():
            if work_file.name == '.gitignore':
                continue
            assert work_file.suffix in ('.svg', '.png')
            work_file.unlink()

    def diff_pixel(self, actual_pixel, expected_pixel):
        ar, ag, ab, aa = actual_pixel
        er, eg, eb, ea = expected_pixel
        if actual_pixel != expected_pixel:
            self.mismatch_found = True
            # Colour
            dr = 0xff
            dg = (ag + eg) // 5
            db = (ab + eb) // 5

            # Opacity
            da = 0xff
        else:
            # Colour
            dr, dg, db = ar, ag, ab

            # Opacity
            da = aa // 3
        return dr, dg, db, da

    def assert_equal(self,
                     svg_actual: Drawing,
                     svg_expected: Drawing,
                     name: str):
        png_actual = drawing_to_image(svg_actual)
        png_expected = drawing_to_image(svg_expected)
        w = max(png_actual.width, png_expected.width)
        h = max(png_actual.height, png_expected.height)

        png_actual_padded = Image.new(png_actual.mode, (w, h))
        png_expected_padded = Image.new(png_expected.mode, (w, h))
        png_actual_padded.paste(png_actual)
        png_expected_padded.paste(png_expected)
        png_diff = Image.new(png_actual.mode, (w, h))
        self.mismatch_found = False
        png_diff.putdata([self.diff_pixel(actual_pixel, expected_pixel)
                          for actual_pixel, expected_pixel in zip(
                            png_actual_padded.getdata(),
                            png_expected_padded.getdata())])

        # Display image when in live turtle mode.
        display_image = getattr(Turtle, 'display_image', None)
        if display_image is not None:
            t = Turtle()
            try:
                # noinspection PyUnresolvedReferences
                w = t.screen.cv.cget('width')
                # noinspection PyUnresolvedReferences
                h = t.screen.cv.cget('height')
                ox, oy = w/2, h/2
                text_height = 20
                t.penup()
                t.goto(-ox, oy)
                t.right(90)
                t.forward(text_height)
                t.write(f'Actual')
                display_image(ox+t.xcor(), oy-t.ycor(),
                              image=encode_image(png_actual))
                t.forward(png_actual.height)
                t.forward(text_height)
                t.write(f'Diff')
                display_image(ox+t.xcor(), oy-t.ycor(),
                              image=encode_image(png_diff))
                t.forward(png_diff.height)
                t.forward(text_height)
                t.write('Expected')
                display_image(ox+t.xcor(), oy-t.ycor(),
                              image=encode_image(png_expected))
                t.forward(png_expected.height)
            except Exception as ex:
                t.write(str(ex))

        if not self.mismatch_found:
            return
        text_actual = svg_actual.asSvg()
        (self.work_dir / (name+'_actual.svg')).write_text(text_actual)
        text_expected = svg_expected.asSvg()
        (self.work_dir / (name+'_expected.svg')).write_text(text_expected)
        with (self.work_dir / (name+'_diff.png')) as f:
            png_diff.save(f)
        assert text_actual == text_expected


def drawing_to_image(drawing: Drawing) -> Image:
    png = drawing.rasterize()
    png_bytes = BytesIO(png.pngData)
    image = Image.open(png_bytes)
    return image


def encode_image(image: Image) -> bytes:
    writer = BytesIO()
    image.save(writer, format='PNG')
    encoded = standard_b64encode(writer.getvalue())
    return encoded.decode('UTF-8')


@pytest.fixture(scope='session')
def svg_differ():
    return SvgDiffer()


def test_summarize_labels():
    figure = Figure()
    figure.add(Track(1, 200, label="Foo"))
    figure.add(Track(1, 200, label="Bar", color='none'))
    expected_summary = """\
Foo[1-200]
Bar(1-200)
"""

    summary = summarize_figure(figure)

    assert summary == expected_summary


def test_summarize_label_objects():
    figure = Figure()
    figure.add(Track(0, 0, label=Label(25, "Foo:")))
    figure.add(Track(0, 0, label="Bar:"))
    expected_summary = """\
Foo:
Bar:
"""

    summary = summarize_figure(figure)

    assert summary == expected_summary


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

    assert summary == expected_summary


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

    assert summary == expected_summary


def test_summarize_regions():
    figure = Figure()
    figure.add(Track(1, 200, label="Foo", regions=[(50, 100, 'lightgreen'),
                                                   (110, 120, 'red')]))
    expected_summary = """\
Foo[1-200], lightgreen{50-100}, red{110-120}
"""

    summary = summarize_figure(figure)

    assert summary == expected_summary


def test_summarize_coverage():
    figure = Figure()
    figure.add(Coverage(10, 20, [11, 11, 21, 1, 1, 1]), gap=-4)
    figure.add(Track(12, 22, label="Bar"))
    expected_summary = """\
Coverage 11, 11, 21, 1, 1, 1
Bar[12-22]
"""

    summary = summarize_figure(figure)

    assert summary == expected_summary


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

    assert summary == expected_summary


def test_summarize_smooth_coverage_ten_percent():
    figure = Figure()
    figure.add(SmoothCoverage(10, 20, [100, 110, 111, 50]), gap=-4)
    figure.add(Track(12, 22, label="Bar"))
    expected_summary = """\
Coverage 100x2, 111, 50
Bar[12-22]
"""

    summary = summarize_figure(figure)

    assert summary == expected_summary


def test_summarize_arrow():
    figure = Figure()
    figure.add(Arrow(10, 50, label='Foo'))
    figure.add(Arrow(60, 30, label='Bar'))
    expected_summary = """\
10--Foo->50
30<-Bar--60
"""

    summary = summarize_figure(figure)

    assert summary == expected_summary


def test_summarize_arrow_group():
    figure = Figure()
    figure.add(ArrowGroup([Arrow(10, 50, label='Foo'),
                           Arrow(60, 30, label='Bar')]))
    figure.add(ArrowGroup([Arrow(1, 50, label='Baz'),
                           Arrow(90, 100, label='Boom')]))
    expected_summary = """\
10--Foo->50, 30<-Bar--60
1--Baz->50, 90--Boom->100
"""

    summary = summarize_figure(figure)

    assert summary == expected_summary


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


def test_plot_genome_coverage_offset_blast():
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
    blast_csv = StringIO("""\
contig_num,ref_name,score,match,pident,start,end,ref_start,ref_end
1,HCV-1a,30,0.33,90,10,12,4,6
""")
    expected_figure = """\
5'[4-344], C[345-917], E1[918-1493], E2[1494-2582], p7[2583-2771], \
NS2[2772-3422], NS3[3423-5315], NS4b[5478-6260], NS4a[5316-5477], \
NS5a[6261-7604], NS5b[7605-9377], 3'[9378-9649]
7--1.1->9
7--1.1->9
Coverage 5x2, 7, 5x6
[1-9], 1-HCV-1a - depth 7(1-9649), lightgreen{7-9}
[4-503], [1004-1503], [2004-2503], [3004-3503], [4004-4503], \
[5004-5503], [6004-6503], [7004-7503], [8004-8503], [9004-9503], \
Partial Blast Results(4-9649)
Coverage 6x3
[4-6], 2-unknown-partial - depth 6(1-9649)
"""

    figure = build_coverage_figure(genome_coverage_csv, blast_csv)

    assert summarize_figure(figure) == expected_figure


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


def test_plot_genome_coverage_unmapped():
    genome_coverage_csv = StringIO("""\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,ins,dels,coverage,link
1-HCV-1a,HCV-1a,1,1,0,0,5,M
1-HCV-1a,HCV-1a,2,2,0,0,5,M
1-HCV-1a,HCV-1a,3,3,0,0,5,M
1-HCV-1a,HCV-1a,4,4,0,0,6,M
1-HCV-1a,HCV-1a,5,5,0,0,6,U
1-HCV-1a,HCV-1a,6,6,0,0,6,U
1-HCV-1a,HCV-1a,7,7,0,0,7,U
1-HCV-1a,HCV-1a,8,8,0,0,7,U
1-HCV-1a,HCV-1a,9,9,0,0,7,U
1-HCV-1a,HCV-1a,10,10,0,0,8,M
1-HCV-1a,HCV-1a,11,11,0,0,8,M
1-HCV-1a,HCV-1a,12,12,0,0,8,M
""")
    expected_figure = """\
5'[1-341], C[342-914], E1[915-1490], E2[1491-2579], p7[2580-2768], \
NS2[2769-3419], NS3[3420-5312], NS4b[5475-6257], NS4a[5313-5474], \
NS5a[6258-7601], NS5b[7602-9374], 3'[9375-9646]
Coverage 5x3, 6x3, 7x3, 8x3
[1-12], 1-HCV-1a - depth 8(1-9646), yellow{5-9}
"""

    figure = build_coverage_figure(genome_coverage_csv)

    assert summarize_figure(figure) == expected_figure


def test_plot_genome_coverage_two_contigs():
    genome_coverage_csv = StringIO("""\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage
1-HCV-1a,HCV-1a,1,1,0,5
1-HCV-1a,HCV-1a,2,2,0,5
1-HCV-1a,HCV-1a,3,3,0,7
1-HCV-1a,HCV-1a,4,4,0,5
1-HCV-1a,HCV-1a,5,5,0,5
1-HCV-1a,HCV-1a,6,6,0,5
contig-1-HCV-1a,HCV-1a,1,1,,
contig-1-HCV-1a,HCV-1a,2,2,,
contig-1-HCV-1a,HCV-1a,3,3,,
contig-1-HCV-1a,HCV-1a,4,4,,
contig-1-HCV-1a,HCV-1a,5,5,,
contig-1-HCV-1a,HCV-1a,6,6,,
2-HCV-1b-partial,,1,,0,5
2-HCV-1b-partial,,2,,0,5
2-HCV-1b-partial,,3,,0,29
2-HCV-1b-partial,,4,,0,5
2-HCV-1b-partial,,5,,0,5
2-HCV-1b-partial,,6,,0,5
3-HCV-1a,HCV-1a,101,101,0,15
3-HCV-1a,HCV-1a,102,102,0,15
3-HCV-1a,HCV-1a,103,103,0,17
3-HCV-1a,HCV-1a,104,104,0,15
3-HCV-1a,HCV-1a,105,105,0,15
3-HCV-1a,HCV-1a,106,106,0,15
contig-3-HCV-1a,HCV-1a,1,101,,
contig-3-HCV-1a,HCV-1a,2,102,,
contig-3-HCV-1a,HCV-1a,3,103,,
contig-3-HCV-1a,HCV-1a,4,104,,
contig-3-HCV-1a,HCV-1a,5,105,,
contig-3-HCV-1a,HCV-1a,6,106,,
""")
    expected_figure = """\
5'[1-341], C[342-914], E1[915-1490], E2[1491-2579], p7[2580-2768], \
NS2[2769-3419], NS3[3420-5312], NS4b[5475-6257], NS4a[5313-5474], \
NS5a[6258-7601], NS5b[7602-9374], 3'[9375-9646]
Coverage 15x2, 17, 15x3
[101-106], 3-HCV-1a - depth 17(1-9646)
[101-106], contig-3-HCV-1a(1-9646)
Coverage 5x2, 7, 5x3
[1-6], 1-HCV-1a - depth 7(1-9646)
[1-6], contig-1-HCV-1a(1-9646)
[1-500], [1001-1500], [2001-2500], [3001-3500], [4001-4500], \
[5001-5500], [6001-6500], [7001-7500], [8001-8500], [9001-9500], \
Partial Blast Results(1-9646)
Coverage 5x2, 29, 5x3
[1-6], 2-HCV-1b-partial - depth 29(1-9646)
"""

    figure = build_coverage_figure(genome_coverage_csv)

    assert summarize_figure(figure) == expected_figure


def test_plot_genome_coverage_blast():
    genome_coverage_csv = StringIO("""\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,ins,dels,coverage
1-HCV-1a,HCV-1a,1,8001,0,0,5
1-HCV-1a,HCV-1a,2,8002,0,0,5
1-HCV-1a,HCV-1a,3,8003,0,0,7
1-HCV-1a,HCV-1a,4,8004,0,0,5
1-HCV-1a,HCV-1a,5,8005,0,0,5
1-HCV-1a,HCV-1a,6,8006,0,0,5
""")
    blast_csv = StringIO("""\
contig_num,ref_name,score,match,pident,start,end,ref_start,ref_end
1,HCV-1g,30,0.33,90,1,2,5001,5002
1,HCV-1a,40,0.33,100,5,6,7006,7005
1,HCV-1a,50,0.5,100,1,3,8001,8003
""")
    expected_figure = """\
5'[1-341], C[342-914], E1[915-1490], E2[1491-2579], p7[2580-2768], \
NS2[2769-3419], NS3[3420-5312], NS4b[5475-6257], NS4a[5313-5474], \
NS5a[6258-7601], NS5b[7602-9374], 3'[9375-9646]
7005<-1.2--7006, 8001--1.1->8003
8001--1.1->8003, 8005--1.2->8006
Coverage 5x2, 7, 5x3
[8001-8006], 1-HCV-1a - depth 7(1-9646)
"""

    figure = build_coverage_figure(genome_coverage_csv, blast_csv)

    assert summarize_figure(figure) == expected_figure


def test_plot_genome_coverage_blast_past_end():
    genome_coverage_csv = StringIO("""\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,ins,dels,coverage
1-HCV-1a,HCV-1a,1,8001,0,0,5
1-HCV-1a,HCV-1a,2,8002,0,0,5
1-HCV-1a,HCV-1a,3,8003,0,0,7
1-HCV-1a,HCV-1a,4,8004,0,0,5
1-HCV-1a,HCV-1a,5,8005,0,0,5
1-HCV-1a,HCV-1a,6,8006,0,0,5
""")
    blast_csv = StringIO("""\
contig_num,ref_name,score,match,pident,start,end,ref_start,ref_end
1,HCV-1g,30,0.33,90,1,2,5001,5002
1,HCV-1a,40,0.33,100,5,10,7010,7005
1,HCV-1a,50,0.5,100,1,3,8001,8003
""")
    expected_figure = """\
5'[1-341], C[342-914], E1[915-1490], E2[1491-2579], p7[2580-2768], \
NS2[2769-3419], NS3[3420-5312], NS4b[5475-6257], NS4a[5313-5474], \
NS5a[6258-7601], NS5b[7602-9374], 3'[9375-9646]
7005<-1.2--7010, 8001--1.1->8003
8001--1.1->8003, 8005--1.2->8006
Coverage 5x2, 7, 5x3
[8001-8006], 1-HCV-1a - depth 7(1-9646)
"""

    figure = build_coverage_figure(genome_coverage_csv, blast_csv)

    assert summarize_figure(figure) == expected_figure


def test_plot_genome_coverage_blast_past_start():
    genome_coverage_csv = StringIO("""\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,ins,dels,coverage
1-HCV-1a,HCV-1a,3,8003,0,0,7
1-HCV-1a,HCV-1a,4,8004,0,0,5
1-HCV-1a,HCV-1a,5,8005,0,0,5
1-HCV-1a,HCV-1a,6,8006,0,0,5
""")
    blast_csv = StringIO("""\
contig_num,ref_name,score,match,pident,start,end,ref_start,ref_end
1,HCV-1g,30,0.33,90,1,2,5001,5002
1,HCV-1a,40,0.33,100,5,6,7006,7005
1,HCV-1a,50,0.5,100,1,3,8001,8003
""")
    expected_figure = """\
5'[1-341], C[342-914], E1[915-1490], E2[1491-2579], p7[2580-2768], \
NS2[2769-3419], NS3[3420-5312], NS4b[5475-6257], NS4a[5313-5474], \
NS5a[6258-7601], NS5b[7602-9374], 3'[9375-9646]
7005<-1.2--7006, 8001--1.1->8003
8003--1.1->8003, 8005--1.2->8006
Coverage 7, 5x3
[8003-8006], 1-HCV-1a - depth 7(1-9646)
"""

    figure = build_coverage_figure(genome_coverage_csv, blast_csv)

    assert summarize_figure(figure) == expected_figure


def test_plot_genome_coverage_blast_aligns_refs():
    genome_coverage_csv = StringIO("""\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,ins,dels,coverage
1-HIV1-G-CM-KP718923-seed,HIV1-B-FR-K03455-seed,1,2261,0,0,5
1-HIV1-G-CM-KP718923-seed,HIV1-B-FR-K03455-seed,2,2262,0,0,5
1-HIV1-G-CM-KP718923-seed,HIV1-B-FR-K03455-seed,3,2263,0,0,5
1-HIV1-G-CM-KP718923-seed,HIV1-B-FR-K03455-seed,4,2264,0,0,5
1-HIV1-G-CM-KP718923-seed,HIV1-B-FR-K03455-seed,5,2265,0,0,5
1-HIV1-G-CM-KP718923-seed,HIV1-B-FR-K03455-seed,6,2266,0,0,5
1-HIV1-G-CM-KP718923-seed,HIV1-B-FR-K03455-seed,7,2267,0,0,5
1-HIV1-G-CM-KP718923-seed,HIV1-B-FR-K03455-seed,8,2268,0,0,5
1-HIV1-G-CM-KP718923-seed,HIV1-B-FR-K03455-seed,9,2269,0,0,5
1-HIV1-G-CM-KP718923-seed,HIV1-B-FR-K03455-seed,10,2270,0,0,5
""")
    blast_csv = StringIO("""\
contig_num,ref_name,score,match,pident,start,end,ref_start,ref_end
1,HIV1-G-CM-KP718923-seed,300,1,90,1,10,1653,1662
""")
    expected_figure = """\
5' LTR[1-634], gag[789-2289], vif[5040-5616], tat[8379-8469], nef[8797-9417]
tat[5831-6045], vpu[6061-6307], rev[8379-8653], 3' LTR[9086-9719]
pol[2085-5096], vpr[5558-5847], rev[5970-6045], env[6225-8795]
PR[2252-2549], RT[2549-3869], INT[4229-5093], V3[7109-7217], GP41[7757-8792]
2261--1.1->2270
2261--1.1->2270
Coverage 5x10
[2261-2270], 1-HIV1-G-CM-KP718923-seed - depth 5(1-9719)
"""

    figure = build_coverage_figure(genome_coverage_csv, blast_csv)

    assert summarize_figure(figure) == expected_figure


def test_plot_genome_coverage_blast_collision():
    """ Two blast results end at the same position. """
    genome_coverage_csv = StringIO("""\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,ins,dels,coverage
1-HCV-1a,HCV-1a,1,8001,0,0,5
1-HCV-1a,HCV-1a,2,8002,0,0,5
1-HCV-1a,HCV-1a,3,8003,0,0,7
1-HCV-1a,HCV-1a,4,8004,0,0,5
1-HCV-1a,HCV-1a,5,8005,0,0,5
1-HCV-1a,HCV-1a,6,8006,0,0,5
""")
    blast_csv = StringIO("""\
contig_num,ref_name,score,match,pident,start,end,ref_start,ref_end
1,HCV-1g,30,0.33,90,1,2,5001,5002
1,HCV-1a,40,0.33,100,3,6,7003,7006
1,HCV-1a,50,0.5,100,1,6,8001,8006
""")
    expected_figure = """\
5'[1-341], C[342-914], E1[915-1490], E2[1491-2579], p7[2580-2768], \
NS2[2769-3419], NS3[3420-5312], NS4b[5475-6257], NS4a[5313-5474], \
NS5a[6258-7601], NS5b[7602-9374], 3'[9375-9646]
7003--1.2->7006, 8001--1.1->8006
8001--1.1->8006, 8003--1.2->8006
Coverage 5x2, 7, 5x3
[8001-8006], 1-HCV-1a - depth 7(1-9646)
"""

    figure = build_coverage_figure(genome_coverage_csv, blast_csv)

    assert summarize_figure(figure) == expected_figure


# noinspection DuplicatedCode
def test_arrow(svg_differ):
    f, expected_svg = start_drawing(200, 55)
    expected_svg.append(Line(0, 20, 168, 20, stroke='black'))
    expected_svg.append(Circle(175/2, 20, 10, stroke='black', fill='ivory'))
    expected_svg.append(Text('1.2',
                             11,
                             175/2, 20,
                             text_anchor='middle',
                             dy="0.35em"))
    expected_svg.append(Lines(175, 20,
                              168, 23.5,
                              168, 16.5,
                              175, 20,
                              fill='black'))
    f.add(Arrow(0, 175, h=20, label='1.2'))
    svg = f.show()

    svg_differ.assert_equal(svg, expected_svg, 'test_arrow')


# noinspection DuplicatedCode
def test_arrow_bottom(svg_differ):
    f, expected_svg = start_drawing(200, 55)
    expected_svg.append(Line(0, 10, 168, 10, stroke='black'))
    expected_svg.append(Circle(175/2, 20, 10, stroke='black', fill='ivory'))
    expected_svg.append(Text('1.2',
                             11,
                             175/2, 20,
                             text_anchor='middle',
                             dy="0.35em"))
    expected_svg.append(Lines(175, 10,
                              168, 13.5,
                              168, 6.5,
                              175, 10,
                              fill='black'))
    f.add(Arrow(0, 175, h=20, elevation=-1, label='1.2'))
    svg = f.show()

    svg_differ.assert_equal(svg, expected_svg, 'test_arrow')


# noinspection DuplicatedCode
def test_reverse_arrow(svg_differ):
    f, expected_svg = start_drawing(200, 55)
    expected_svg.append(Line(7, 10, 175, 10, stroke='black'))
    expected_svg.append(Circle(175/2, 20, 10, stroke='black', fill='ivory'))
    expected_svg.append(Text('X',
                             11,
                             175/2, 20,
                             text_anchor='middle',
                             dy="0.35em"))
    expected_svg.append(Lines(0, 10,
                              7, 13.5,
                              7, 6.5,
                              0, 10,
                              fill='black'))
    f.add(Arrow(175, 0, h=20, elevation=-1, label='X'))
    svg = f.show()

    svg_differ.assert_equal(svg, expected_svg, 'test_arrow')


# noinspection DuplicatedCode
def test_scaled_arrow(svg_differ):
    expected_svg = Drawing(100, 35, origin=(0, 0))
    expected_svg.append(Line(0, 10, 93, 10, stroke='black'))
    expected_svg.append(Circle(50, 20, 10, stroke='black', fill='ivory'))
    expected_svg.append(Text('2.3',
                             11,
                             50, 20,
                             text_anchor='middle',
                             dy="0.35em"))
    expected_svg.append(Lines(100, 10,
                              93, 13.5,
                              93, 6.5,
                              100, 10,
                              fill='black'))

    f = Figure()
    f.add(Arrow(0, 200, h=20, elevation=-1, label='2.3'))
    svg = f.show(w=100)

    svg_differ.assert_equal(svg, expected_svg, 'test_arrow')


# noinspection DuplicatedCode
def test_small_arrow(svg_differ):
    f, expected_svg = start_drawing(200, 55)
    expected_svg.append(Line(100, 10, 125, 10, stroke='black'))
    expected_svg.append(Circle(116, 20, 10, stroke='black', fill='ivory'))
    expected_svg.append(Text('2.3',
                             11,
                             116, 20,
                             text_anchor='middle',
                             dy="0.35em"))
    expected_svg.append(Lines(132, 10,
                              125, 13.5,
                              125, 6.5,
                              132, 10,
                              fill='black'))

    f.add(Arrow(100, 132, h=20, elevation=-1, label='2.3'))
    svg = f.show()

    svg_differ.assert_equal(svg, expected_svg, 'test_arrow')


# noinspection DuplicatedCode
def test_tiny_arrow(svg_differ):
    f, expected_svg = start_drawing(200.0, 55.0)
    expected_svg.append(Circle(102, 20, 10, stroke='black', fill='ivory'))
    expected_svg.append(Text('2.3',
                             11,
                             102, 20,
                             text_anchor='middle',
                             dy="0.35em"))
    expected_svg.append(Lines(104, 10,
                              100, 13.5,
                              100, 6.5,
                              104, 10,
                              fill='black'))

    f.add(Arrow(100, 104, h=20, elevation=-1, label='2.3'))
    svg = f.show()

    svg_differ.assert_equal(svg, expected_svg, 'test_arrow')


# noinspection DuplicatedCode
def test_tiny_arrow_at_edge(svg_differ):
    expected_svg = Drawing(210, 35, origin=(0, 0))
    expected_svg.append(Circle(197.5, 20, 10, stroke='black', fill='ivory'))
    expected_svg.append(Text('2.3',
                             11,
                             197.5, 20,
                             text_anchor='middle',
                             dy="0.35em"))
    expected_svg.append(Lines(200, 10,
                              195, 13.5,
                              195, 6.5,
                              200, 10,
                              fill='black'))

    f = Figure()
    f.add(ArrowGroup([Arrow(195, 200, h=20, elevation=-1, label='2.3')]))
    svg = f.show()

    svg_differ.assert_equal(svg, expected_svg, 'test_arrow')


def start_drawing(width, height):
    expected_svg = Drawing(width, height, origin=(0, 0))
    expected_svg.append(Rectangle(0, height-15,
                                  200, 10,
                                  stroke='lightgrey',
                                  fill='lightgrey'))
    expected_svg.append(Text('Header',
                             10,
                             width/2, height-15,
                             font_family='monospace',
                             text_anchor='middle'))
    f = Figure()
    f.add(Track(0, width, label='Header'))
    return f, expected_svg


def test_arrow_repr(svg_differ):
    arrow = Arrow(175, 0, label='X')

    assert repr(arrow) == "Arrow(175, 0, label='X')"


# noinspection DuplicatedCode
def test_arrow_group(svg_differ):
    expected_figure = Figure()
    expected_figure.add(Track(1, 500, label='Header'))
    h = 30
    expected_figure.add(Arrow(1, 200, label='X', h=h), gap=-h)
    expected_figure.add(Arrow(300, 500, label='Y', h=h))
    expected_svg = expected_figure.show()

    f = Figure()
    f.add(Track(1, 500, label='Header'))
    f.add(ArrowGroup([Arrow(1, 200, label='X', h=h),
                      Arrow(300, 500, label='Y', h=h)]))
    svg = f.show()

    svg_differ.assert_equal(svg, expected_svg, 'test_arrow_group')


# noinspection DuplicatedCode
def test_arrow_group_overlap(svg_differ):
    expected_figure = Figure()
    expected_figure.add(Track(1, 500, label='Header'))
    h = 20
    expected_figure.add(Arrow(1, 300, label='X', h=h), gap=3)
    expected_figure.add(Arrow(1, 300, label='Y', h=h))
    expected_svg = expected_figure.show()

    f = Figure()
    f.add(Track(1, 500, label='Header'))
    f.add(ArrowGroup([Arrow(1, 300, label='X', h=h),
                      Arrow(1, 300, label='Y', h=h)]))
    svg = f.show()

    svg_differ.assert_equal(svg, expected_svg, 'test_arrow_group')


# noinspection DuplicatedCode
def test_arrow_group_reverse_overlap(svg_differ):
    expected_figure = Figure()
    expected_figure.add(Track(1, 500, label='Header'))
    h = 20
    expected_figure.add(Arrow(1, 300, label='X', h=h), gap=3)
    expected_figure.add(Arrow(400, 250, label='Y', h=h))
    expected_svg = expected_figure.show()

    f = Figure()
    f.add(Track(1, 500, label='Header'))
    f.add(ArrowGroup([Arrow(1, 300, label='X', h=h),
                      Arrow(400, 250, label='Y', h=h)]))
    svg = f.show()

    svg_differ.assert_equal(svg, expected_svg, 'test_arrow_group')


# noinspection DuplicatedCode
def test_arrow_group_small_neighbour(svg_differ):
    expected_figure = Figure()
    expected_figure.add(Track(1, 500, label='Header'))
    h = 20
    expected_figure.add(Arrow(301, 315, elevation=-1, label='1.2', h=h), gap=-h)
    expected_figure.add(Arrow(1, 300, elevation=-1, label='1.1', h=h))
    expected_svg = expected_figure.show()

    f = Figure()
    f.add(Track(1, 500, label='Header'))
    f.add(ArrowGroup([Arrow(1, 300, elevation=-1, label='1.1', h=h),
                      Arrow(301, 315, elevation=-1, label='1.2', h=h)]))
    svg = f.show()

    svg_differ.assert_equal(svg, expected_svg, 'test_arrow_group')


def test_draw_coverage(svg_differ):
    expected_figure = Figure()
    expected_figure.add(Track(0, 1, color='', h=-4))  # Just a spacer.
    expected_figure.add(Track(100, 200, label='Bar'))
    expected_svg = expected_figure.show()
    expected_svg.insert(0, draw.Rectangle(100, 20, 25, 5, fill='blue'))
    expected_svg.insert(1, draw.Rectangle(125, 20, 25, 10, fill='blue'))
    expected_svg.insert(2, draw.Rectangle(175, 20, 25, 1, fill='blue'))

    figure = Figure()
    coverage_depths = 25 * [5] + 25 * [10] + 25 * [0] + 25 * [1]
    figure.add(SmoothCoverage(100, 200, coverage_depths), gap=-4)
    figure.add(Track(100, 200, label="Bar"))

    svg = figure.show()

    svg_differ.assert_equal(svg, expected_svg, 'test_draw_coverage')
