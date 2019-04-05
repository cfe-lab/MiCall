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
    def build_config(sequences):
        projects = ProjectConfig()
        projects.config = {
            'projects': {
                'PROJ1': {
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
        sequences = {'KX1234': 'ACTGATA'}  # {accession_num: seq}
        project_config = self.build_config({'P-A-KX-KX1234-seed': 'ACTGATA'})
        expected_report = """\
Matching references: P-A-KX-KX1234-seed
"""

        report, error_count = compare_config(['P-A-KX-KX1234-seed'],
                                             project_config,
                                             sequences,
                                             name_part=3)

        assert expected_report == report
        assert 0 == error_count

    def test_diff(self):
        source_sequences = {'R1': 'ACTGATT'}
        project_config = self.build_config({'P-A-KX-R1-seed': 'ACTGATA'})
        expected_report = """\
ERROR: changed references:
P-A-KX-R1-seed
- ACTGATT
?       ^
+ ACTGATA
?       ^

"""

        report, error_count = compare_config(['P-A-KX-R1-seed'],
                                             project_config,
                                             source_sequences,
                                             name_part=3)

        assert expected_report == report
        assert 1 == error_count

    def test_multiple_seeds(self):
        sequences = {'R1': 'ACTGATT', 'R2': 'TTT'}
        project_config = self.build_config({'P-A-KX-R1-seed': 'ACTGATT',
                                            'P-A-KY-R2-seed': 'TTT'})
        expected_report = """\
Matching references: P-A-KX-R1-seed, P-A-KY-R2-seed
"""

        report, error_count = compare_config(['P-A-KX-R1-seed', 'P-A-KY-R2-seed'],
                                             project_config,
                                             sequences,
                                             name_part=3)

        assert expected_report == report
        assert 0 == error_count

    def test_dashes(self):
        sequences = {'R1': 'ACT---GATT'}
        project_config = self.build_config({'P-A-KX-R1-seed': 'ACTGATT'})
        expected_report = """\
Matching references: P-A-KX-R1-seed
"""

        report, error_count = compare_config(['P-A-KX-R1-seed'],
                                             project_config,
                                             sequences,
                                             name_part=3)

        assert expected_report == report
        assert 0 == error_count

    def test_missing_from_source(self):
        sequences = {'R1': 'ACTGATT'}
        project_config = self.build_config({'P-A-KX-R1-seed': 'ACTGATT',
                                            'P-A-KY-R2-seed': 'TTT'})
        expected_report = """\
Matching references: P-A-KX-R1-seed
ERROR: references missing from source: P-A-KY-R2-seed
"""

        report, error_count = compare_config(['P-A-KX-R1-seed', 'P-A-KY-R2-seed'],
                                             project_config,
                                             sequences,
                                             name_part=3)

        assert expected_report == report
        assert 1 == error_count

    def test_many_matches(self):
        sequences = {'R1': 'ACTGATT', 'R2': 'TTT', 'R3': 'A', 'R4': 'A'}
        project_config = self.build_config({'P-A-KX-R1-seed': 'ACTGATT',
                                            'P-A-KY-R2-seed': 'TTT',
                                            'P-A-KY-R3-seed': 'A',
                                            'P-A-KY-R4-seed': 'A'})
        expected_report = """\
Matching references: P-A-KX-R1-seed, P-A-KY-R2-seed, P-A-KY-R3-seed,
P-A-KY-R4-seed
"""

        report, error_count = compare_config(['P-A-KX-R1-seed',
                                              'P-A-KY-R2-seed',
                                              'P-A-KY-R3-seed',
                                              'P-A-KY-R4-seed'],
                                             project_config,
                                             sequences,
                                             name_part=3)

        assert expected_report == report
        assert 0 == error_count

    def test_many_missing(self):
        sequences = {'R1': 'ACTGATT'}
        project_config = self.build_config({'P-A-KX-R1-seed': 'ACTGATT',
                                            'P-A-KY-R2-seed': 'TTT',
                                            'P-A-KY-R3-seed': 'A',
                                            'P-A-KY-R4-seed': 'A'})
        expected_report = """\
Matching references: P-A-KX-R1-seed
ERROR: references missing from source: P-A-KY-R2-seed, P-A-KY-R3-seed,
P-A-KY-R4-seed
"""

        report, error_count = compare_config(['P-A-KX-R1-seed',
                                              'P-A-KY-R2-seed',
                                              'P-A-KY-R3-seed',
                                              'P-A-KY-R4-seed'],
                                             project_config,
                                             sequences,
                                             name_part=3)

        assert expected_report == report
        assert 3 == error_count
