from micall.utils.report_amino import ReportNucleotide


def test():
    report_nuc = ReportNucleotide(10)

    assert repr(report_nuc) == 'ReportNucleotide(10, SeedNucleotide({}))'
