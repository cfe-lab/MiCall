from io import StringIO

import pytest

from micall.data.landmark_reader import LandmarkReader


def test_gene_with_end():
    landmarks_yaml = StringIO("""\
- seed_pattern: R1
  coordinates: R1-seed
  landmarks:
    - {name: gene1, start: 1, end: 634, frame: 0, colour: darkgrey}
    - {name: gene2, start: 789, end: 2289, frame: 0, colour: lightblue}
    - {name: gene3, start: 5040, end: 5616, frame: 0, colour: steelblue}
- seed_pattern: R2
  coordinates: R2-seed
  landmarks:
    - {name: gene1, start: 1001, end: 1634, frame: 0, colour: darkgrey}
    - {name: gene2, start: 1789, end: 3289, frame: 0, colour: lightblue}
""")
    expected_gene = dict(name='gene2',
                         start=789,
                         end=2289,
                         frame=0,
                         colour='lightblue')
    reader = LandmarkReader.load(landmarks_yaml)

    gene = reader.get_gene('R1-seed', 'gene2', drop_stop_codon=False)

    assert gene == expected_gene


def test_gene_without_end():
    landmarks_yaml = StringIO("""\
- seed_pattern: R1
  coordinates: R1-seed
  landmarks:
    - {name: gene1, start: 1, frame: 0, colour: darkgrey}
    - {name: gene2, start: 789, frame: 0, colour: lightblue}
    - {name: gene3, start: 5040, end: 5616, frame: 0, colour: steelblue}
- seed_pattern: R2
  coordinates: R2-seed
  landmarks:
    - {name: gene1, start: 1001, frame: 0, colour: darkgrey}
    - {name: gene2, start: 1789, end: 3289, frame: 0, colour: lightblue}
""")
    expected_gene = dict(name='gene2',
                         start=789,
                         end=5040,
                         frame=0,
                         colour='lightblue')
    reader = LandmarkReader.load(landmarks_yaml)

    gene = reader.get_gene('R1-seed', 'gene2', drop_stop_codon=False)

    assert gene == expected_gene


def test_load_defaults():
    expected_gene = dict(name='PR',
                         start=2252,
                         end=2549,
                         frame=3,
                         colour='orange')
    reader = LandmarkReader.load()

    gene = reader.get_gene('HIV1-B-FR-K03455-seed', 'PR', drop_stop_codon=False)

    assert gene == expected_gene


def test_gene_with_prefix():
    landmarks_yaml = StringIO("""\
- seed_pattern: R1
  coordinates: R1-seed
  prefix: R1-group-
  landmarks:
    - {name: gene1, start: 1, frame: 0, colour: darkgrey}
    - {name: gene2, start: 789, frame: 0, colour: lightblue}
    - {name: gene3, start: 5040, end: 5616, frame: 0, colour: steelblue}
""")
    expected_gene = dict(name='gene2',
                         start=789,
                         end=5040,
                         frame=0,
                         colour='lightblue')
    reader = LandmarkReader.load(landmarks_yaml)

    gene = reader.get_gene('R1-seed', 'R1-group-gene2', drop_stop_codon=False)

    assert gene == expected_gene


def test_gene_drop_stop_codon():
    landmarks_yaml = StringIO("""\
- seed_pattern: R1
  coordinates: R1-seed
  landmarks:
    - {name: gene1, start: 1, frame: 0, colour: darkgrey}
    - {name: gene2, start: 789, frame: 0, colour: lightblue}
    - {name: gene3, start: 5040, end: 5616, frame: 0, colour: steelblue}
""")
    expected_gene = dict(name='gene2',
                         start=789,
                         end=5037,
                         frame=0,
                         colour='lightblue')
    reader = LandmarkReader.load(landmarks_yaml)

    gene = reader.get_gene('R1-seed', 'gene2', drop_stop_codon=True)

    assert gene == expected_gene


def test_gene_missing_prefix():
    landmarks_yaml = StringIO("""\
- seed_pattern: R1
  coordinates: R1-seed
  prefix: R1-group-
  landmarks:
    - {name: gene1, start: 1, frame: 0, colour: darkgrey}
    - {name: gene2, start: 789, frame: 0, colour: lightblue}
    - {name: gene3, start: 5040, end: 5616, frame: 0, colour: steelblue}
""")
    reader = LandmarkReader.load(landmarks_yaml)

    with pytest.raises(ValueError, match=r"Gene name 'gene2' does not start "
                                         r"with prefix 'R1-group-'."):
        reader.get_gene('R1-seed', 'gene2')


def test_gene_with_full_name():
    landmarks_yaml = StringIO("""\
- seed_pattern: R1
  coordinates: R1-seed
  prefix: R1-group-
  landmarks:
    - {name: '1', full_name: gene1, start: 1, frame: 0, colour: darkgrey}
    - {name: '2', full_name: gene2, start: 789, frame: 0, colour: lightblue}
    - {name: '3', full_name: gene3, start: 5040, end: 5616, frame: 0, colour: steelblue}
""")
    expected_gene = dict(name='2',
                         full_name='gene2',
                         start=789,
                         end=5040,
                         frame=0,
                         colour='lightblue')
    reader = LandmarkReader.load(landmarks_yaml)

    gene = reader.get_gene('R1-seed', 'R1-group-gene2', drop_stop_codon=False)

    assert gene == expected_gene


def test_get_coordinates_found():
    landmarks_yaml = StringIO("""\
- seed_pattern: R1.*
  coordinates: R1a
  landmarks:
    - {name: '1', full_name: gene1, start: 1, frame: 0, colour: darkgrey}
""")
    reader = LandmarkReader.load(landmarks_yaml)

    coordinates = reader.get_coordinates('R1b')

    assert coordinates == 'R1a'


def test_get_coordinates_unknown():
    landmarks_yaml = StringIO("""\
- seed_pattern: R1.*
  coordinates: R1a
  landmarks:
    - {name: '1', full_name: gene1, start: 1, frame: 0, colour: darkgrey}
""")
    reader = LandmarkReader.load(landmarks_yaml)

    with pytest.raises(ValueError, match=r"No landmarks match 'R2b'"):
        reader.get_coordinates('R2b')


def test_get_gene_unknown():
    landmarks_yaml = StringIO("""\
- seed_pattern: R1.*
  coordinates: R1a
  landmarks:
    - {name: '1', full_name: gene1, start: 1, frame: 0, colour: darkgrey}
""")
    reader = LandmarkReader.load(landmarks_yaml)

    with pytest.raises(ValueError, match=r"Coordinates unknown: 'R2b'"):
        reader.get_gene('R2b', '2')
