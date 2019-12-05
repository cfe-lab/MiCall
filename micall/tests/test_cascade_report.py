from io import StringIO
from unittest import TestCase

from micall.core.cascade_report import CascadeReport


class CascadeReportTest(TestCase):
    def test_no_inputs(self):
        report = CascadeReport(StringIO())
        expected_report = """\
demultiplexed,v3loop,g2p,prelim_map,remap,aligned
0,0,0,0,0,0
"""

        report.generate()

        self.assertEqual(expected_report, report.cascade_csv.getvalue())

    def test_g2p_inputs(self):
        report = CascadeReport(StringIO())
        report.g2p_summary_csv = StringIO("""\
mapped,valid,X4calls,X4pct,final,validpct
100,90,0,0.00,,90.00
""")
        expected_report = """\
demultiplexed,v3loop,g2p,prelim_map,remap,aligned
100,100,90,0,0,0
"""

        report.generate()

        self.assertEqual(expected_report, report.cascade_csv.getvalue())

    def test_remap_inputs(self):
        report = CascadeReport(StringIO())
        report.remap_counts_csv = StringIO("""\
type,count,ignored
raw,300,x
prelim R1-seed,200,99
prelim *,100,
remap-1 R1-seed,220,
remap-final R1-seed,240,
unmapped,60,
""")
        expected_report = """\
demultiplexed,v3loop,g2p,prelim_map,remap,aligned
150,0,0,100,120,0
"""

        report.generate()

        self.assertEqual(expected_report, report.cascade_csv.getvalue())

    def test_g2p_inputs_remapped(self):
        report = CascadeReport(StringIO(), is_g2p_remapped=True)
        report.g2p_summary_csv = StringIO("""\
mapped,valid,X4calls,X4pct,final,validpct
100,90,0,0.00,,90.00
""")
        report.remap_counts_csv = StringIO("""\
type,count,ignored
raw,300,x
prelim R1-seed,200,99
prelim *,100,
remap-1 R1-seed,220,
remap-final R1-seed,240,
unmapped,60,
""")
        expected_report = """\
demultiplexed,v3loop,g2p,prelim_map,remap,aligned
150,100,90,100,120,0
"""

        report.generate()

        self.assertEqual(expected_report, report.cascade_csv.getvalue())

    def test_remap_inputs_two_seeds(self):
        report = CascadeReport(StringIO())
        report.remap_counts_csv = StringIO("""\
type,count,ignored
raw,300,x
prelim R1-seed,200,99
prelim R2-seed,20,99
prelim *,80,
remap-1 R1-seed,220,
remap-1 R2-seed,22,
remap-final R1-seed,240,
remap-final R2-seed,24,
unmapped,36,
""")
        expected_report = """\
demultiplexed,v3loop,g2p,prelim_map,remap,aligned
150,0,0,110,132,0
"""

        report.generate()

        self.assertEqual(expected_report, report.cascade_csv.getvalue())

    def test_aligned_inputs(self):
        report = CascadeReport(StringIO())
        report.aligned_csv = StringIO("""\
ignored,count
x,50
9,40
""")
        expected_report = """\
demultiplexed,v3loop,g2p,prelim_map,remap,aligned
0,0,0,0,0,90
"""

        report.generate()

        self.assertEqual(expected_report, report.cascade_csv.getvalue())

    def test_all_inputs(self):
        report = CascadeReport(StringIO())
        report.g2p_summary_csv = StringIO("""\
mapped,valid,X4calls,X4pct,final,validpct
90,80,0,0.00,,90.00
""")
        report.remap_counts_csv = StringIO("""\
type,count,ignored
raw,300,x
prelim R1-seed,200,99
prelim *,100,
remap-1 R1-seed,220,
remap-final R1-seed,240,
unmapped,60,
""")
        report.aligned_csv = StringIO("""\
ignored,count
x,50
9,40
""")
        expected_report = """\
demultiplexed,v3loop,g2p,prelim_map,remap,aligned
240,90,80,100,120,90
"""

        report.generate()

        self.assertEqual(expected_report, report.cascade_csv.getvalue())
