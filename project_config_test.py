import unittest
import StringIO
import project_config

class ProjectConfigurationTest(unittest.TestCase):
    def setUp(self):
        self.defaultJsonIO = StringIO.StringIO("""\
{
  "projects": {
    "R1": {
      "regions": [
        {
          "coordinate_region": "R1",
          "seed_region": "R1-seed"
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
        "RWN",
        "NWR"
      ]
    }
  }
}
""")
        self.config = project_config.ProjectConfig()
        
    def testConvert(self):
        expected_fasta = """\
>R1-seed
ACTGAAAGGG
"""
        fasta = StringIO.StringIO()

        self.config.load(self.defaultJsonIO)
        self.config.writeSeedFasta(fasta)
        
        self.assertMultiLineEqual(expected_fasta, fasta.getvalue())

    def testSharedRegions(self):
        jsonIO = StringIO.StringIO("""\
{
  "projects": {
    "R1": {
      "regions": [
        {
          "coordinate_region": null,
          "seed_region": "R1-seed"
        }
      ]
    },
    "R1 and R2": {
      "regions": [
        {
          "coordinate_region": null,
          "seed_region": "R1-seed"
        },
        {
          "coordinate_region": null,
          "seed_region": "R2-seed"
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
        fasta = StringIO.StringIO()

        self.config.load(jsonIO)
        self.config.writeSeedFasta(fasta)
        
        self.assertMultiLineEqual(expected_fasta, fasta.getvalue())

    def testUnusedRegion(self):
        jsonIO = StringIO.StringIO("""\
{
  "projects": {
    "R1": {
      "regions": [
        {
          "coordinate_region": null,
          "seed_region": "R1-seed"
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
        fasta = StringIO.StringIO()

        self.config.load(jsonIO)
        self.config.writeSeedFasta(fasta)
        
        self.assertMultiLineEqual(expected_fasta, fasta.getvalue())

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
        
    def testReload(self):
        jsonIO1 = StringIO.StringIO("""\
{
  "projects": {
    "R1": {
      "regions": [
        {
          "coordinate_region": null,
          "seed_region": "R1-seed"
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
        jsonIO2 = StringIO.StringIO("""\
{
  "projects": {
    "R2": {
      "regions": [
        {
          "coordinate_region": null,
          "seed_region": "R2-seed"
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
