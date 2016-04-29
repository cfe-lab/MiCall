from cStringIO import StringIO
from unittest import TestCase
from micall.core.filter_quality import report_bad_cycles


class FilterQualityTest(TestCase):
    def setUp(self):
        self.addTypeEqualityFunc(str, self.assertMultiLineEqual)

    def test_good(self):
        quality_csv = StringIO("""\
tile,cycle,errorrate
2,1,1.0
2,2,7.49
""")
        expected_bad_cycles_csv = """\
tile,cycle,errorrate
"""
        bad_cycles_csv = StringIO()

        report_bad_cycles(quality_csv, bad_cycles_csv)

        self.assertEqual(expected_bad_cycles_csv, bad_cycles_csv.getvalue())

    def test_bad(self):
        quality_csv = StringIO("""\
tile,cycle,errorrate
2,1,1.0
2,2,7.5
""")
        expected_bad_cycles_csv = """\
tile,cycle,errorrate
2,2,7.5
"""
        bad_cycles_csv = StringIO()

        report_bad_cycles(quality_csv, bad_cycles_csv)

        self.assertEqual(expected_bad_cycles_csv, bad_cycles_csv.getvalue())

    def test_missing(self):
        quality_csv = StringIO("""\
tile,cycle,errorrate
2,1,1.0
2,2
""")
        expected_bad_cycles_csv = """\
tile,cycle,errorrate
2,2,
"""
        bad_cycles_csv = StringIO()

        report_bad_cycles(quality_csv, bad_cycles_csv)

        self.assertEqual(expected_bad_cycles_csv, bad_cycles_csv.getvalue())

    def test_blank(self):
        quality_csv = StringIO("""\
tile,cycle,errorrate
2,1,1.0
2,2,
""")
        expected_bad_cycles_csv = """\
tile,cycle,errorrate
2,2,
"""
        bad_cycles_csv = StringIO()

        report_bad_cycles(quality_csv, bad_cycles_csv)

        self.assertEqual(expected_bad_cycles_csv, bad_cycles_csv.getvalue())

    def test_filter_following(self):
        quality_csv = StringIO("""\
tile,cycle,errorrate
2,1,7.5
2,2,1.0
""")
        expected_bad_cycles_csv = """\
tile,cycle,errorrate
2,1,7.5
2,2,1.0
"""
        bad_cycles_csv = StringIO()

        report_bad_cycles(quality_csv, bad_cycles_csv)

        self.assertEqual(expected_bad_cycles_csv, bad_cycles_csv.getvalue())

    def test_tile_count(self):
        quality_csv = StringIO("""\
tile,cycle,errorrate
1,1,7.5
1,-1,7.5
2,1,1.0
2,-1,7.5
""")
        expected_bad_tiles_csv = """\
tile,bad_cycles
1,2
2,1
"""
        bad_cycles_csv = StringIO()
        bad_tiles_csv = StringIO()

        report_bad_cycles(quality_csv, bad_cycles_csv, bad_tiles_csv)

        self.assertEqual(expected_bad_tiles_csv, bad_tiles_csv.getvalue())
