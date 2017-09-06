from io import StringIO
from unittest import TestCase

from micall.core.project_config import ProjectConfig
from micall.utils.fetch_sequences import extract_sequence, parse_compendium, \
    compare_config


class FetchSequencesTest(TestCase):
    def test(self):
        source = """\
preceding text
ORIGIN
       1 actggattca tcagtttgca
      21 taccg
//
following text
"""
        expected_sequence = 'ACTGGATTCATCAGTTTGCATACCG'

        sequence = extract_sequence(source)

        self.assertEqual(expected_sequence, sequence)


class ParseCompendiumTest(TestCase):
    def test_single(self):
        fasta = StringIO("""\
>A.B.C.R1
ACTGA
""")
        expected_sequences = {'R1': 'ACTGA'}

        sequences = parse_compendium(fasta)

        self.assertEqual(expected_sequences, sequences)

    def test_no_dots(self):
        fasta = StringIO("""\
>R1
ACTGA
""")
        expected_sequences = {'R1': 'ACTGA'}

        sequences = parse_compendium(fasta)

        self.assertEqual(expected_sequences, sequences)

    def test_multiple_sequences(self):
        fasta = StringIO("""\
>A.B.C.R1
ACTGA
>X.Y.R2
GATA
""")
        expected_sequences = {'R1': 'ACTGA', 'R2': 'GATA'}

        sequences = parse_compendium(fasta)

        self.assertEqual(expected_sequences, sequences)

    def test_multiple_lines(self):
        fasta = StringIO("""\
>A.B.C.R1
ACTGA
TTACA
""")
        expected_sequences = {'R1': 'ACTGATTACA'}

        sequences = parse_compendium(fasta)

        self.assertEqual(expected_sequences, sequences)


class CompareConfigTest(TestCase):
    @staticmethod
    def build_config(project, sequences):
        projects = ProjectConfig()
        projects.config = {
            'projects': {
                project: {
                    'regions': [
                        {
                            'seed_region_names': list(sequences.keys())
                        }
                    ]
                }
            },
            'regions': {name: {'reference': [sequence]}
                        for name, sequence in sequences.items()}
        }
        return projects

    def test_match(self):
        sequences = {'R1': 'ACTGATA'}
        project_config = self.build_config('PROJ1',
                                           {'P-A-KX-R1-seed': 'ACTGATA'})
        expected_report = """\
P-A-KX-R1-seed
==

0 differences
"""

        report = ''.join(compare_config('PROJ1', project_config, sequences))

        self.assertEqual(expected_report, report)

    def test_diff(self):
        sequences = {'R1': 'ACTGATT'}
        project_config = self.build_config('PROJ1',
                                           {'P-A-KX-R1-seed': 'ACTGATA'})
        expected_report = """\
P-A-KX-R1-seed
- ACTGATA
?       ^
+ ACTGATT
?       ^

1 differences
"""

        report = ''.join(compare_config('PROJ1', project_config, sequences))

        self.assertEqual(expected_report, report)

    def test_multiple_seeds(self):
        sequences = {'R1': 'ACTGATT', 'R2': 'TTT'}
        project_config = self.build_config('PROJ1',
                                           {'P-A-KX-R1-seed': 'ACTGATT',
                                            'P-A-KY-R2-seed': 'TTT'})
        expected_report = """\
P-A-KX-R1-seed
==

P-A-KY-R2-seed
==

0 differences
"""

        report = ''.join(compare_config('PROJ1', project_config, sequences))

        self.assertEqual(expected_report, report)

    def test_dashes(self):
        sequences = {'R1': 'ACT---GATT'}
        project_config = self.build_config('PROJ1',
                                           {'P-A-KX-R1-seed': 'ACTGATT'})
        expected_report = """\
P-A-KX-R1-seed
==

0 differences
"""

        report = ''.join(compare_config('PROJ1', project_config, sequences))

        self.assertEqual(expected_report, report)

    def test(self):
        sequences = {'R1': 'ACTGATT'}
        project_config = self.build_config('PROJ1',
                                           {'P-A-KX-R1-seed': 'ACTGATT',
                                            'P-A-KY-R2-seed': 'TTT'})
        expected_report = """\
P-A-KX-R1-seed
==

P-A-KY-R2-seed
missing

1 differences
"""

        report = ''.join(compare_config('PROJ1', project_config, sequences))

        self.assertEqual(expected_report, report)
