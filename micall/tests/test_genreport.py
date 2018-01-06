from unittest import TestCase

from micall.hivdb.genreport import ReportPage


class GenReportTest(TestCase):
    def test_unknown_keys(self):
        with self.assertRaisesRegex(
                ValueError,
                r"Unknown configuration: bogus_key, weird_key."):
            ReportPage(dict(weird_key="Hello",
                            known_regions=['R1', 'R2'],
                            bogus_key=42))

    def test_missing_keys(self):
        with self.assertRaisesRegex(
                ValueError,
                r"Missing configuration: generated_by_text, known_drug_classes, "
                r"known_drugs, report_title, resistance_level_colours."):
            ReportPage(dict(known_regions=['R1', 'R2'],
                            disclaimer_text="Hello."),
                       raise_missing=True)

    def test_register_regions(self):
        page1 = ReportPage(dict(known_regions=['R1', 'R2']))
        page2 = ReportPage(dict(known_regions=['R3']))
        expected_regions = {'R1': page1,
                            'R2': page1,
                            'R3': page2}
        regions = {}

        page1.register_regions(regions)
        page2.register_regions(regions)

        self.assertEqual(expected_regions, regions)

    def test_register_drug_classes(self):
        page1 = ReportPage(dict(known_drug_classes=[('C1', 'Class 1'),
                                                    ('C2', 'Class 2')],
                                known_drugs={'C1': [('D1', 'Drug 1')],
                                             'C2': [('D2', 'Drug 2')]}))
        page2 = ReportPage(dict(known_drug_classes=[('C3', 'Class 3')],
                                known_drugs={'C3': [('D3', 'Drug 3')]}))
        expected_drug_classes = {'C1': page1,
                                 'C2': page1,
                                 'C3': page2}
        drug_classes = {}

        page1.register_drug_classes(drug_classes)
        page2.register_drug_classes(drug_classes)

        self.assertEqual(expected_drug_classes, drug_classes)

    def test_repr(self):
        page = ReportPage({'report_title': 'Example Report',
                           'known_regions': ['R1', 'R2']})
        expected_repr = "ReportPage({'report_title': 'Example Report'})"

        r = repr(page)

        self.assertEqual(expected_repr, r)

    def test_get_reported_drug_classes(self):
        page = ReportPage(dict(known_drug_classes=[('C1', 'Class 1'),
                                                   ('C2', 'Class 2')],
                               known_drugs={'C1': [('D1', 'Drug 1')],
                                            'C2': [('D2', 'Drug 2')]}))
        page.resistance_calls['D1'] = 'Some resistance data'
        expected_drug_classes = {'C1'}

        reported_drug_classes = page.get_reported_drug_classes()

        self.assertEqual(expected_drug_classes, reported_drug_classes)