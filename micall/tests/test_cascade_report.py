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
0,100,90,0,0,0
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
        report = CascadeReport(StringIO())
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
150,90,80,100,120,90
"""

        report.generate()

        self.assertEqual(expected_report, report.cascade_csv.getvalue())


    def test_denovo_multimapping_reports_physical_pair_count(self):
        report = CascadeReport(StringIO())
        report.remap_counts_csv = StringIO("""\
type,count
raw,4
remap-final A,2
remap-final B,2
remap-fragments,1
unmapped,2
""")
        report.generate()
        self.assertEqual(
            'demultiplexed,v3loop,g2p,prelim_map,remap,aligned\n'
            '2,0,0,0,1,0\n',
            report.cascade_csv.getvalue())



def test_denovo_cascade_uses_explicit_unique_fragments_with_mapped_singletons():
    from io import StringIO
    from micall.core.cascade_report import CascadeReport

    report = CascadeReport(StringIO())
    report.remap_counts_csv = StringIO(
        'type,count\n'
        'raw,8\n'
        'remap-final A,5\n'
        'remap-final B,2\n'
        'remap-fragments,3\n'
        'unmapped,4\n')
    report.generate()
    assert report.counts['demultiplexed'] == 4
    assert report.counts['remap'] == 3

def test_unique_aligned_fragment_sidecar_prevents_alt_placement_inflation():
    from io import StringIO
    report = CascadeReport(StringIO())
    report.remap_counts_csv = StringIO(
        'type,count\n'
        'raw,8\n'
        'remap-final A,10\n'
        'remap-fragments,4\n')
    report.aligned_csv = StringIO('count\n1\n4\n7\n')
    report.unique_aligned_csv = StringIO(
        'unique_aligned_fragments\n3\n')
    report.generate()
    assert report.counts['demultiplexed'] == 4
    assert report.counts['remap'] == 4
    assert report.counts['aligned'] == 3

def test_unique_aligned_sidecar_rejects_impossible_clinical_qc_counts():
    from io import StringIO
    import pytest

    report = CascadeReport(StringIO())
    report.remap_counts_csv = StringIO(
        'type,count\nraw,8\nremap-fragments,3\n')
    report.unique_aligned_csv = StringIO(
        'unique_aligned_fragments\n4\n')
    with pytest.raises(ValueError, match='exceed mapped'):
        report.generate()


def test_unique_aligned_sidecar_rejects_duplicate_clinical_report_rows():
    from io import StringIO
    import pytest

    report = CascadeReport(StringIO())
    report.remap_counts_csv = StringIO(
        'type,count\nraw,8\nremap-fragments,3\n')
    report.unique_aligned_csv = StringIO(
        'unique_aligned_fragments\n2\n2\n')
    with pytest.raises(ValueError, match='one data row'):
        report.generate()


def test_denovo_fragment_total_overrides_per_contig_counts_in_any_row_order():
    from io import StringIO
    for rows in (
            ['remap-fragments,3', 'remap-final A,16', 'remap-final B,16'],
            ['remap-final A,16', 'remap-final B,16', 'remap-fragments,3']):
        report = CascadeReport(StringIO())
        report.remap_counts_csv = StringIO(
            'type,count\nraw,8\n' + '\n'.join(rows) + '\n')
        report.generate()
        assert report.counts['remap'] == 3
