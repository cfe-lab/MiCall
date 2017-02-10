import json
from io import StringIO
import unittest

from micall.monitor.projects_dump import check_key_positions


class CheckKeyPositionsTest(unittest.TestCase):
    def setUp(self):
        self.warningIO = StringIO()

    def testSingleRegion(self):
        projects = json.loads("""\
{
    "R1": {
      "max_variants": 5,
      "regions": [
        {
          "coordinate_region": "R1",
          "key_positions": [
              {
                  "end_pos": null,
                  "start_pos": 42
              }
          ],
          "seed_region_names": [
            "R1-seed"
          ]
        }
      ]
    }
}
""")
        expected_warnings = ""

        check_key_positions(projects, self.warningIO)

        self.assertMultiLineEqual(expected_warnings, self.warningIO.getvalue())

    def testMultipleRegionSingleSetOfPositions(self):
        projects = json.loads("""\
{
    "R1": {
      "max_variants": 5,
      "regions": [
        {
          "coordinate_region": "R1",
          "key_positions": [
              {
                  "end_pos": null,
                  "start_pos": 42
              }
          ],
          "seed_region_names": [
            "R1a-seed"
          ]
        },
        {
          "coordinate_region": "R1",
          "key_positions": [],
          "seed_region_names": [
            "R1b-seed"
          ]
        }
      ]
    }
}
""")
        expected_warnings = ""

        check_key_positions(projects, self.warningIO)

        self.assertMultiLineEqual(expected_warnings, self.warningIO.getvalue())

    def testMultipleRegionTwoSetsOfPositions(self):
        projects = json.loads("""\
{
    "R1": {
      "max_variants": 5,
      "regions": [
        {
          "coordinate_region": "R1",
          "key_positions": [
              {
                  "end_pos": null,
                  "start_pos": 42
              }
          ],
          "seed_region_names": [
            "R1a-seed"
          ]
        },
        {
          "coordinate_region": "R1",
          "key_positions": [
              {
                  "end_pos": 99,
                  "start_pos": 1
              }
          ],
          "seed_region_names": [
            "R1b-seed"
          ]
        }
      ]
    }
}
""")
        expected_warnings = (
            "WARNING: project R1 has multiple sets of key positions for " +
            "coordinate region R1.\n")

        check_key_positions(projects, self.warningIO)

        self.assertMultiLineEqual(expected_warnings, self.warningIO.getvalue())

    def testMultipleRegionThreeSetsOfPositions(self):
        projects = json.loads("""\
{
    "R1": {
      "max_variants": 5,
      "regions": [
        {
          "coordinate_region": "R1",
          "key_positions": [
              {
                  "end_pos": null,
                  "start_pos": 42
              }
          ],
          "seed_region_names": [
            "R1a-seed"
          ]
        },
        {
          "coordinate_region": "R1",
          "key_positions": [
              {
                  "end_pos": 99,
                  "start_pos": 1
              }
          ],
          "seed_region_names": [
            "R1b-seed"
          ]
        },
        {
          "coordinate_region": "R1",
          "key_positions": [
              {
                  "end_pos": 14,
                  "start_pos": 3
              }
          ],
          "seed_region_names": [
            "R1c-seed"
          ]
        }
      ]
    }
}
""")
        expected_warnings = (
            "WARNING: project R1 has multiple sets of key positions for " +
            "coordinate region R1.\n")

        check_key_positions(projects, self.warningIO)

        self.assertMultiLineEqual(expected_warnings, self.warningIO.getvalue())

    def testDuplicateSeedAndCoordinate(self):
        projects = json.loads("""\
{
    "R1": {
      "max_variants": 5,
      "regions": [
        {
          "coordinate_region": "R1",
          "key_positions": [
              {
                  "end_pos": null,
                  "start_pos": 42
              }
          ],
          "seed_region_names": [
            "R1-seed"
          ]
        },
        {
          "coordinate_region": "R1",
          "key_positions": [],
          "seed_region_names": [
            "R1-seed"
          ]
        }
      ]
    }
}
""")
        expected_warnings = (
            "WARNING: project R1 has duplicate seed and coordinate: R1-seed, R1\n")

        check_key_positions(projects, self.warningIO)

        self.assertMultiLineEqual(expected_warnings, self.warningIO.getvalue())
