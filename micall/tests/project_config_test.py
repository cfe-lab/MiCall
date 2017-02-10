import unittest
from io import StringIO
from micall.core.project_config import ProjectConfig


class ProjectConfigurationTest(unittest.TestCase):
    def setUp(self):
        self.defaultJsonIO = StringIO("""\
{
  "projects": {
    "R1": {
      "max_variants": 5,
      "regions": [
        {
          "coordinate_region": "R1",
          "seed_region_names": ["R1-seed"],
          "id": 10042
        }
      ]
    }
  },
  "regions": {
    "R1-seed": {
      "is_nucleotide": true,
      "reference": [
        "ACTGAAA",
        "GGG"
      ],
      "seed_group": "R1-seeds"
    },
    "R1": {
      "is_nucleotide": false,
      "reference": [
        "RWN",
        "NWR"
      ],
      "seed_group": null
    }
  }
}
""")
        self.config = ProjectConfig()

    def testConvert(self):
        expected_fasta = """\
>R1-seed
ACTGAAAGGG
"""
        fasta = StringIO()

        self.config.load(self.defaultJsonIO)
        self.config.writeSeedFasta(fasta)

        self.assertMultiLineEqual(expected_fasta, fasta.getvalue())

    def testSharedRegions(self):
        jsonIO = StringIO("""\
{
  "projects": {
    "R1": {
      "regions": [
        {
          "coordinate_region": null,
          "seed_region_names": ["R1-seed"]
        }
      ]
    },
    "R1 and R2": {
      "regions": [
        {
          "coordinate_region": null,
          "seed_region_names": ["R1-seed"]
        },
        {
          "coordinate_region": null,
          "seed_region_names": ["R2-seed"]
        }
      ]
    }
  },
  "regions": {
    "R1-seed": {
      "is_nucleotide": true,
      "reference": [
        "ACTGAAA",
        "GGG"
      ]
    },
    "R2-seed": {
      "is_nucleotide": true,
      "reference": [
        "TTT"
      ]
    }
  }
}
""")
        expected_fasta = """\
>R1-seed
ACTGAAAGGG
>R2-seed
TTT
"""
        fasta = StringIO()

        self.config.load(jsonIO)
        self.config.writeSeedFasta(fasta)

        self.assertMultiLineEqual(expected_fasta, fasta.getvalue())

    def testUnusedRegion(self):
        jsonIO = StringIO("""\
{
  "projects": {
    "R1": {
      "regions": [
        {
          "coordinate_region": null,
          "seed_region_names": ["R1-seed"]
        }
      ]
    }
  },
  "regions": {
    "R1-seed": {
      "is_nucleotide": true,
      "reference": [
        "ACTGAAA",
        "GGG"
      ]
    },
    "R2-seed": {
      "is_nucleotide": true,
      "reference": [
        "TTT"
      ]
    }
  }
}
""")
        expected_fasta = """\
>R1-seed
ACTGAAAGGG
"""
        fasta = StringIO()

        self.config.load(jsonIO)
        self.config.writeSeedFasta(fasta)

        self.assertMultiLineEqual(expected_fasta, fasta.getvalue())

    def testExcludeSeeds(self):
        jsonIO = StringIO("""\
{
  "projects": {
    "R1": {
      "regions": [
        {
          "coordinate_region": null,
          "seed_region_names": ["R1-seed"]
        }
      ]
    },
    "R2": {
      "regions": [
        {
          "coordinate_region": null,
          "seed_region_names": ["R2-seed"]
        }
      ]
    },
    "R3": {
      "regions": [
        {
          "coordinate_region": null,
          "seed_region_names": ["R3-seed"]
        }
      ]
    }
  },
  "regions": {
    "R1-seed": {
      "is_nucleotide": true,
      "reference": [
        "ACTGAAA",
        "GGG"
      ]
    },
    "R2-seed": {
      "is_nucleotide": true,
      "reference": [
        "TTT"
      ]
    },
    "R3-seed": {
      "is_nucleotide": true,
      "reference": [
        "TAG"
      ]
    }
  }
}
""")
        expected_fasta = """\
>R2-seed
TTT
"""
        fasta = StringIO()

        self.config.load(jsonIO)
        self.config.writeSeedFasta(fasta, excluded_seeds=['R1-seed', 'R3-seed'])

        self.assertMultiLineEqual(expected_fasta, fasta.getvalue())

    def testExcludeUnknownSeed(self):
        expected_fasta = """\
>R1-seed
ACTGAAAGGG
"""
        fasta = StringIO()

        self.config.load(self.defaultJsonIO)
        self.config.writeSeedFasta(fasta, excluded_seeds=['R99-seed'])

        self.assertMultiLineEqual(expected_fasta, fasta.getvalue())

    def testDuplicateReference(self):
        jsonIO = StringIO("""\
{
  "projects": {
    "R1": {
      "regions": [
        {
          "coordinate_region": null,
          "seed_region_names": ["R1a-seed", "R1b-seed"]
        }
      ]
    }
  },
  "regions": {
    "R1a-seed": {
      "is_nucleotide": true,
      "reference": [
        "ACTAAAGGG"
      ]
    },
    "R1b-seed": {
      "is_nucleotide": true,
      "reference": [
        "ACTAAAGGG"
      ]
    }
  }
}
""")
        fasta = StringIO()
        self.config.load(jsonIO)

        self.assertRaisesRegex(RuntimeError,
                               "Duplicate references: R1a-seed and R1b-seed.",
                               self.config.writeSeedFasta,
                               fasta)

    def testGetReference(self):
        self.config.load(self.defaultJsonIO)
        seed_name = 'R1-seed'
        expected_ref = 'ACTGAAAGGG'

        seed_ref = self.config.getReference(seed_name)

        self.assertSequenceEqual(expected_ref, seed_ref)

    def testGetCoordinateReferences(self):
        self.config.load(self.defaultJsonIO)
        seed_name = 'R1-seed'
        expected_refs = {'R1': 'RWNNWR'}

        coordinate_refs = self.config.getCoordinateReferences(seed_name)

        self.assertDictEqual(expected_refs, coordinate_refs)

    def testUnknownReference(self):
        self.config.load(self.defaultJsonIO)
        seed_name = 'R-unknown'

        self.assertRaises(KeyError, self.config.getReference, seed_name)

    def testMaxVariants(self):
        self.config.load(self.defaultJsonIO)
        coordinate_region_name = 'R1'

        self.assertEqual(5, self.config.getMaxVariants(coordinate_region_name))

    def testMaxVariantsUnusedRegion(self):
        jsonIO = StringIO("""\
{
  "projects": {
    "R1": {
      "max_variants": 2,
      "regions": [
        {
          "coordinate_region": "R1",
          "seed_region_names": ["R1-seed"]
        }
      ]
    }
  },
  "regions": {
    "R1-seed": {
      "is_nucleotide": true,
      "reference": [
        "ACTGAAA",
        "GGG"
      ]
    },
    "R1": {
      "is_nucleotide": false,
      "reference": [
        "NSFW"
      ]
    },
    "R2": {
      "is_nucleotide": false,
      "reference": [
        "RSW"
      ]
    }
  }
}
""")
        self.config.load(jsonIO)
        coordinate_region_name = 'R2'

        self.assertEqual(0, self.config.getMaxVariants(coordinate_region_name))

    def testMaxVariantsTwoProjects(self):
        """ If two projects specify a maximum for the same coordinate region,
        use the bigger of the two.
        """
        jsonIO = StringIO("""\
{
  "projects": {
    "R1": {
      "max_variants": 9,
      "regions": [
        {
          "coordinate_region": "R1",
          "seed_region_names": ["R1-seed"]
        }
      ]
    },
    "R1-and-R2": {
      "max_variants": 2,
      "regions": [
        {
          "coordinate_region": "R1",
          "seed_region_names": ["R1-seed"]
        },
        {
          "coordinate_region": "R2",
          "seed_region_names": ["R1-seed"]
        }
      ]
    }
  },
  "regions": {
    "R1-seed": {
      "is_nucleotide": true,
      "reference": [
        "ACTGAAA",
        "GGG"
      ]
    },
    "R1": {
      "is_nucleotide": false,
      "reference": [
        "NSFW"
      ]
    },
    "R2": {
      "is_nucleotide": false,
      "reference": [
        "RSW"
      ]
    }
  }
}
""")
        self.config.load(jsonIO)
        coordinate_region_name = 'R1'

        self.assertEqual(9, self.config.getMaxVariants(coordinate_region_name))

    def testReload(self):
        jsonIO1 = StringIO("""\
{
  "projects": {
    "R1": {
      "regions": [
        {
          "coordinate_region": null,
          "seed_region_names": ["R1-seed"]
        }
      ]
    }
  },
  "regions": {
    "R1-seed": {
      "is_nucleotide": true,
      "reference": [
        "ACTGAAA",
        "GGG"
      ]
    }
  }
}
""")
        jsonIO2 = StringIO("""\
{
  "projects": {
    "R2": {
      "regions": [
        {
          "coordinate_region": null,
          "seed_region_names": ["R2-seed"]
        }
      ]
    }
  },
  "regions": {
    "R2-seed": {
      "is_nucleotide": true,
      "reference": [
        "GACCTA"
      ]
    }
  }
}
""")

        self.config.load(jsonIO1)
        self.config.load(jsonIO2)

        self.assertRaises(KeyError, self.config.getReference, "R1-seed")
        self.assertSequenceEqual("GACCTA", self.config.getReference("R2-seed"))

    def testProjectSeeds(self):
        expected_seeds = set(['R1-seed'])

        self.config.load(self.defaultJsonIO)
        seeds = self.config.getProjectSeeds('R1')

        self.assertSetEqual(expected_seeds, seeds)

    def testSeedGroup(self):
        expected_group = "R1-seeds"

        self.config.load(self.defaultJsonIO)
        group = self.config.getSeedGroup('R1-seed')

        self.assertEqual(expected_group, group)


class ProjectConfigurationProjectRegionsTest(unittest.TestCase):
    def setUp(self):
        self.config = ProjectConfig()
        self.defaultJsonIO = StringIO("""\
{
  "projects": {
    "R1": {
      "max_variants": 0,
      "regions": [
        {
          "coordinate_region": "R1",
          "coordinate_region_length": 3,
          "key_positions": [],
          "min_coverage1": 10,
          "min_coverage2": 50,
          "min_coverage3": 100,
          "seed_region_names": [
            "R1-seed"
          ]
        }
      ]
    },
    "R1 and R2": {
      "max_variants": 0,
      "regions": [
        {
          "coordinate_region": "R1",
          "coordinate_region_length": 3,
          "key_positions": [1, 3],
          "min_coverage1": 10,
          "min_coverage2": 50,
          "min_coverage3": 100,
          "seed_region_names": [
            "R1-seed"
          ]
        },
        {
          "coordinate_region": "R2",
          "coordinate_region_length": 1,
          "key_positions": [],
          "min_coverage1": 10,
          "min_coverage2": 50,
          "min_coverage3": 100,
          "seed_region_names": [
            "R2-seed"
          ]
        }
      ]
    }
  }
}
""")

    def testProjectRegions(self):
        expected_project_regions = [{"project_name": "R1",
                                     "coordinate_region_length": 3,
                                     "key_positions": [],
                                     "min_coverage1": 10,
                                     "min_coverage2": 50,
                                     "min_coverage3": 100},
                                    {"project_name": "R1 and R2",
                                     "coordinate_region_length": 3,
                                     "key_positions": [1, 3],
                                     "min_coverage1": 10,
                                     "min_coverage2": 50,
                                     "min_coverage3": 100}]

        self.config.load(self.defaultJsonIO)
        project_regions = list(self.config.getProjectRegions('R1-seed', 'R1'))

        self.assertEqual(expected_project_regions, project_regions)

    def testProjectExcluded(self):
        excluded_projects = ['R1']
        expected_project_regions = [{"project_name": "R1 and R2",
                                     "coordinate_region_length": 3,
                                     "key_positions": [1, 3],
                                     "min_coverage1": 10,
                                     "min_coverage2": 50,
                                     "min_coverage3": 100}]

        self.config.load(self.defaultJsonIO)
        project_regions = list(self.config.getProjectRegions(
            'R1-seed',
            'R1',
            excluded_projects))

        self.assertEqual(expected_project_regions, project_regions)
