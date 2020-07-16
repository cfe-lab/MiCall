import csv
import sys
from collections import Counter
from csv import DictReader
from io import StringIO

import pytest

from micall.core import project_config
from micall.core.aln2counts import InsertionWriter, SequenceReport, SeedNucleotide
from micall.core.project_config import ProjectConfig

# noinspection PyUnresolvedReferences
from micall.tests.test_remap import load_projects


def prepare_reads(aligned_reads_text):
    full_text = "refname,qcut,rank,count,offset,seq\n" + aligned_reads_text
    dummy_file = StringIO(full_text)
    return csv.DictReader(dummy_file)


@pytest.fixture
def sequence_report():
    return create_sequence_report()


def create_sequence_report():
    projects = project_config.ProjectConfig()
    # Content of seed regions is irrelevant. For R-NO-COORD, there is
    # no coordinate reference, so we use the seed reference for display, but
    # only the length matters.
    projects.load(StringIO("""\
{
"projects": {
"R1": {
  "max_variants": 10,
  "regions": [
    {
      "coordinate_region": "R1",
      "seed_region_names": ["R1-seed"]
    }
  ]
},
"R2": {
  "max_variants": 10,
  "regions": [
    {
      "coordinate_region": "R2",
      "seed_region_names": ["R2-seed"]
    }
  ]
},
"R3": {
  "max_variants": 10,
  "regions": [
    {
      "coordinate_region": "R3",
      "seed_region_names": ["R3-seed"]
    }
  ]
},
"R4": {
  "max_variants": 10,
  "regions": [
    {
      "coordinate_region": "R4",
      "seed_region_names": ["R4-seed"]
    }
  ]
},
"R5": {
  "max_variants": 10,
  "regions": [
    {
      "coordinate_region": "R5",
      "seed_region_names": ["R5-seed"]
    }
  ]
},
"R6": {
  "max_variants": 10,
  "regions": [
    {
      "coordinate_region": "R6",
      "seed_region_names": ["R6a-seed", "R6b-seed"]
    }
  ]
},
"R7": {
  "max_variants": 10,
  "regions": [
    {
      "coordinate_region": "R7a",
      "seed_region_names": ["R7-seed"]
    },
    {
      "coordinate_region": "R7b",
      "seed_region_names": ["R7-seed"]
    }
  ]
}
},
"regions": {
"R1-seed": {
  "is_nucleotide": true,
  "reference": [
    "AAATTTAGG"
  ]
},
"R1": {
  "is_nucleotide": false,
  "reference": [
    "KFR"
  ]
},
"R2-seed": {
  "is_nucleotide": true,
  "reference": [
    "AAATTTGGCCCGAGA"
  ]
},
"R2": {
  "is_nucleotide": false,
  "reference": [
    "KFGPR"
  ]
},
"R3-seed": {
  "is_nucleotide": true,
  "reference": [
    "AAATTTCAGACCCCACGAGAGCAT"
  ]
},
"R3": {
  "is_nucleotide": false,
  "reference": [
    "KFQTPREH"
  ]
},
"R4-seed": {
  "is_nucleotide": true,
  "reference": [
    "ATGGCAAACTCAATCAAT"
  ]
},
"R4": {
  "is_nucleotide": false,
  "reference": [
    "SIN"
  ]
},
"R5-seed": {
  "comment": "Coord has G that's not in seed.",
  "is_nucleotide": true,
  "reference": [
    "AAATTTCCGAGA"
  ]
},
"R5": {
  "is_nucleotide": false,
  "reference": [
    "KFGPR"
  ]
},
"R6a-seed": {
  "is_nucleotide": true,
  "reference": [
    "AAATTTAGG"
  ],
  "seed_group": "R6-seeds"
},
"R6b-seed": {
  "is_nucleotide": true,
  "reference": [
    "GGGAAATTCAGGACAGGGGGGGGG"
  ],
  "seed_group": "R6-seeds"
},
"R6": {
  "is_nucleotide": false,
  "reference": [
    "KFR"
  ]
},
"R7-seed": {
  "is_nucleotide": true,
  "reference": [
    "AAATTTCAGACCCCACGAGAGCAT"
  ]
},
"R7a": {
  "is_nucleotide": false,
  "reference": [
    "KFQ"
  ]
},
"R7b": {
  "is_nucleotide": false,
  "reference": [
    "REH"
  ]
},
"R-NO-COORD": {
  "is_nucleotide": true,
  "reference": [
    "ACTACTACT"
  ]
}
}
}
"""))
    landmarks_yaml = """\
- seed_pattern: R1-.*
  coordinates: R1-seed
  landmarks:
    - {name: ref, start: 1, end: 9, colour: steelblue}
- seed_pattern: R2-.*
  coordinates: R2-seed
  landmarks:
    - {name: ref, start: 1, end: 15, colour: steelblue}
- seed_pattern: R3-.*
  coordinates: R3-seed
  landmarks:
    - {name: a, start: 1, end: 12, colour: lightblue}
    - {name: z, start: 13, end: 24, colour: steelblue}
"""
    conseq_mixture_cutoffs = [0.1]
    report = StubbedSequenceReport(InsertionWriter(StringIO()),
                                   projects,
                                   conseq_mixture_cutoffs,
                                   landmarks_yaml=landmarks_yaml)
    return report


@pytest.fixture
def sequence_report_overlapping_regions(sequence_report):
    sequence_report.projects.load(StringIO("""\
{
  "projects": {
    "R1": {
      "max_variants": 0,
      "regions": [
        {
          "coordinate_region": "R1",
          "seed_region_names": ["R1-seed"]
        },
        {
          "coordinate_region": "R1-expanded",
          "seed_region_names": ["R1-seed"]
        }
      ]
    }
  },
  "regions": {
    "R1-seed": {
      "is_nucleotide": true,
      "reference": [
        "AAATTTAGG"
      ]
    },
    "R1": {
      "is_nucleotide": false,
      "reference": [
        "KFR"
      ]
    },
    "R1-expanded": {
      "is_nucleotide": false,
      "reference": [
        "GPKFREH"
      ]
    }
  }
}
    """))
    return sequence_report


class StubbedSequenceReport(SequenceReport):
    def __init__(self, *args, **kwargs):
        SequenceReport.__init__(self, *args, **kwargs)
        self.overrides = {}

    def _pair_align(self, reference, query, *args, **kwargs):
        override = self.overrides.get((reference, query))
        return (override
                if override is not None
                else SequenceReport._pair_align(self, reference, query))

    def add_override(self,
                     reference,
                     query,
                     aligned_query,
                     aligned_reference,
                     score=sys.maxsize):
        self.overrides[(reference, query)] = (aligned_reference,
                                              aligned_query,
                                              score)


def test_single_read_amino_report(sequence_report):
    """ In this sample, there is a single read with two codons.
    AAA -> K
    TTT -> F
    The coordinate reference has three codons, so the third position is
    empty.
    """
    # refname,qcut,rank,count,offset,seq
    aligned_reads = prepare_reads("""\
R1-seed,15,0,9,0,AAATTT
""")

    expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
"""
    report_file = StringIO()
    sequence_report.write_amino_header(report_file)
    sequence_report.read(aligned_reads)
    sequence_report.write_amino_counts()

    assert report_file.getvalue() == expected_text


def test_single_read_nucleotide_report(sequence_report):
    """ In this sample, there is a single read with two codons.
    AAA -> K
    TTT -> F
    The coordinate reference has three codons, so the third position is
    empty.
    """
    # refname,qcut,rank,count,offset,seq
    aligned_reads = prepare_reads("""\
R1-seed,15,0,9,0,AAATTT
""")

    expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T,N,del,ins,clip,v3_overlap,coverage
R1-seed,R1,15,1,1,9,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,2,2,9,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,3,3,9,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,4,4,0,0,0,9,0,0,0,0,0,9
R1-seed,R1,15,5,5,0,0,0,9,0,0,0,0,0,9
R1-seed,R1,15,6,6,0,0,0,9,0,0,0,0,0,9
R1-seed,R1,15,,7,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,,8,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,,9,0,0,0,0,0,0,0,0,0,0
"""

    report_file = StringIO()
    sequence_report.write_nuc_header(report_file)
    sequence_report.read(aligned_reads)
    sequence_report.write_nuc_counts()

    assert report_file.getvalue() == expected_text


def test_consensus_from_single_read(sequence_report):
    """ In this sample, there is a single read with two codons.
    AAA -> K
    TTT -> F
    """
    # refname,qcut,rank,count,offset,seq
    aligned_reads = prepare_reads("""\
R1-seed,15,0,9,0,AAATTT
""")
    expected_text = """\
region,q-cutoff,consensus-percent-cutoff,offset,sequence
R1-seed,15,MAX,0,AAATTT
R1-seed,15,0.100,0,AAATTT
"""

    report_file = StringIO()
    sequence_report.write_consensus_header(report_file)
    sequence_report.read(aligned_reads)
    sequence_report.write_consensus()

    assert report_file.getvalue() == expected_text


def test_multiple_prefix_nucleotide_report_overlapping_regions(
        sequence_report_overlapping_regions):
    """ Assemble counts from two contigs when coordinate regions overlap.

    Contig 1-R1 AAATTT -> KF
    Contig 2-R1 TTTAGG -> FR

    Contig 1 and 2 should combine into R1 with KFR, but also be reported in
    R1-expanded.
    """
    # refname,qcut,rank,count,offset,seq
    aligned_reads1 = prepare_reads("1-R1-seed,15,0,5,0,AAATTT")
    aligned_reads2 = prepare_reads("2-R1-seed,15,0,2,0,TTTAGG")

    expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,\
A,C,G,T,N,del,ins,clip,v3_overlap,coverage
R1-seed,R1,15,,1,5,0,0,0,0,0,0,0,0,5
R1-seed,R1,15,,2,5,0,0,0,0,0,0,0,0,5
R1-seed,R1,15,,3,5,0,0,0,0,0,0,0,0,5
R1-seed,R1,15,,4,0,0,0,7,0,0,0,0,0,7
R1-seed,R1,15,,5,0,0,0,7,0,0,0,0,0,7
R1-seed,R1,15,,6,0,0,0,7,0,0,0,0,0,7
R1-seed,R1,15,,7,2,0,0,0,0,0,0,0,0,2
R1-seed,R1,15,,8,0,0,2,0,0,0,0,0,0,2
R1-seed,R1,15,,9,0,0,2,0,0,0,0,0,0,2
R1-seed,R1-expanded,15,,1,0,0,0,0,0,0,0,0,0,0
R1-seed,R1-expanded,15,,2,0,0,0,0,0,0,0,0,0,0
R1-seed,R1-expanded,15,,3,0,0,0,0,0,0,0,0,0,0
R1-seed,R1-expanded,15,,4,0,0,0,0,0,0,0,0,0,0
R1-seed,R1-expanded,15,,5,0,0,0,0,0,0,0,0,0,0
R1-seed,R1-expanded,15,,6,0,0,0,0,0,0,0,0,0,0
R1-seed,R1-expanded,15,,7,5,0,0,0,0,0,0,0,0,5
R1-seed,R1-expanded,15,,8,5,0,0,0,0,0,0,0,0,5
R1-seed,R1-expanded,15,,9,5,0,0,0,0,0,0,0,0,5
R1-seed,R1-expanded,15,,10,0,0,0,7,0,0,0,0,0,7
R1-seed,R1-expanded,15,,11,0,0,0,7,0,0,0,0,0,7
R1-seed,R1-expanded,15,,12,0,0,0,7,0,0,0,0,0,7
R1-seed,R1-expanded,15,,13,2,0,0,0,0,0,0,0,0,2
R1-seed,R1-expanded,15,,14,0,0,2,0,0,0,0,0,0,2
R1-seed,R1-expanded,15,,15,0,0,2,0,0,0,0,0,0,2
R1-seed,R1-expanded,15,,16,0,0,0,0,0,0,0,0,0,0
R1-seed,R1-expanded,15,,17,0,0,0,0,0,0,0,0,0,0
R1-seed,R1-expanded,15,,18,0,0,0,0,0,0,0,0,0,0
R1-seed,R1-expanded,15,,19,0,0,0,0,0,0,0,0,0,0
R1-seed,R1-expanded,15,,20,0,0,0,0,0,0,0,0,0,0
R1-seed,R1-expanded,15,,21,0,0,0,0,0,0,0,0,0,0
"""

    report = sequence_report_overlapping_regions
    report_file = StringIO()
    detail_report_file = StringIO()
    report.write_nuc_header(report_file)
    report.write_nuc_detail_header(detail_report_file)
    report.read(aligned_reads1)
    report.write_nuc_detail_counts()
    report.combine_reports()
    report.read(aligned_reads2)
    report.write_nuc_detail_counts()
    report.combine_reports()
    report.write_nuc_counts()

    assert report_file.getvalue() == expected_text


def test_nucleotide_report_excluded_regions(sequence_report_overlapping_regions):
    """ Although the reads align with two coordinate regions, only report one.

    R1 AAATTTAGG -> KFR

    R1-expanded includes KFR, but will be excluded.
    """
    # refname,qcut,rank,count,offset,seq
    aligned_reads = prepare_reads("R1-seed,15,0,5,0,AAATTTAGG")

    expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,\
A,C,G,T,N,del,ins,clip,v3_overlap,coverage
R1-seed,R1,15,1,1,5,0,0,0,0,0,0,0,0,5
R1-seed,R1,15,2,2,5,0,0,0,0,0,0,0,0,5
R1-seed,R1,15,3,3,5,0,0,0,0,0,0,0,0,5
R1-seed,R1,15,4,4,0,0,0,5,0,0,0,0,0,5
R1-seed,R1,15,5,5,0,0,0,5,0,0,0,0,0,5
R1-seed,R1,15,6,6,0,0,0,5,0,0,0,0,0,5
R1-seed,R1,15,7,7,5,0,0,0,0,0,0,0,0,5
R1-seed,R1,15,8,8,0,0,5,0,0,0,0,0,0,5
R1-seed,R1,15,9,9,0,0,5,0,0,0,0,0,0,5
"""

    report = sequence_report_overlapping_regions
    report_file = StringIO()
    report.write_nuc_header(report_file)
    report.read(aligned_reads, excluded_regions={'R1-expanded'})
    report.write_nuc_counts()

    assert report_file.getvalue() == expected_text


def test_nucleotide_report_included_regions(sequence_report_overlapping_regions):
    """ Although the reads align with two coordinate regions, only report one.

    R1-expanded AAATTTAGG -> KFR

    R1 includes KFR, but is not included.
    """
    # refname,qcut,rank,count,offset,seq
    aligned_reads = prepare_reads("R1-seed,15,0,5,0,AAATTTAGG")

    expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,\
A,C,G,T,N,del,ins,clip,v3_overlap,coverage
R1-seed,R1-expanded,15,,1,0,0,0,0,0,0,0,0,0,0
R1-seed,R1-expanded,15,,2,0,0,0,0,0,0,0,0,0,0
R1-seed,R1-expanded,15,,3,0,0,0,0,0,0,0,0,0,0
R1-seed,R1-expanded,15,,4,0,0,0,0,0,0,0,0,0,0
R1-seed,R1-expanded,15,,5,0,0,0,0,0,0,0,0,0,0
R1-seed,R1-expanded,15,,6,0,0,0,0,0,0,0,0,0,0
R1-seed,R1-expanded,15,1,7,5,0,0,0,0,0,0,0,0,5
R1-seed,R1-expanded,15,2,8,5,0,0,0,0,0,0,0,0,5
R1-seed,R1-expanded,15,3,9,5,0,0,0,0,0,0,0,0,5
R1-seed,R1-expanded,15,4,10,0,0,0,5,0,0,0,0,0,5
R1-seed,R1-expanded,15,5,11,0,0,0,5,0,0,0,0,0,5
R1-seed,R1-expanded,15,6,12,0,0,0,5,0,0,0,0,0,5
R1-seed,R1-expanded,15,7,13,5,0,0,0,0,0,0,0,0,5
R1-seed,R1-expanded,15,8,14,0,0,5,0,0,0,0,0,0,5
R1-seed,R1-expanded,15,9,15,0,0,5,0,0,0,0,0,0,5
R1-seed,R1-expanded,15,,16,0,0,0,0,0,0,0,0,0,0
R1-seed,R1-expanded,15,,17,0,0,0,0,0,0,0,0,0,0
R1-seed,R1-expanded,15,,18,0,0,0,0,0,0,0,0,0,0
R1-seed,R1-expanded,15,,19,0,0,0,0,0,0,0,0,0,0
R1-seed,R1-expanded,15,,20,0,0,0,0,0,0,0,0,0,0
R1-seed,R1-expanded,15,,21,0,0,0,0,0,0,0,0,0,0
"""

    report = sequence_report_overlapping_regions
    report_file = StringIO()
    report.write_nuc_header(report_file)
    report.read(aligned_reads, included_regions={'R1-expanded'})
    report.write_nuc_counts()

    assert report_file.getvalue() == expected_text


# noinspection DuplicatedCode
def test_duplicated_sars_base_amino(sequence_report):
    """ Special case for duplicated base in SARS orf1ab.

    Expect amino sequence AQSFLNRVCG.
    """

    # refname,qcut,rank,count,offset,seq
    aligned_reads = prepare_reads("""\
SARS-CoV-2-seed,15,0,9,0,GCACAATCGTTTTTAAACGGGTTTGCGGTGTAAGTGCAGCCCGTCTTACAC
""")
    # Repeat is here:                     ^

    #                                       A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,...,coverage
    expected_text = """\
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,1,4396,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,4,4397,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,7,4398,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,10,4399,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,13,4400,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,16,4401,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,18,4402,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,21,4403,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,24,4404,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,27,4405,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9"""
    sequence_report.projects = ProjectConfig.loadDefault()
    orf1ab_size = len(sequence_report.projects.getReference('SARS-CoV-2-orf1ab'))
    nsp12_size = len(sequence_report.projects.getReference('SARS-CoV-2-nsp12'))

    report_file = StringIO()
    sequence_report.write_amino_header(report_file)
    sequence_report.read(aligned_reads)
    sequence_report.write_amino_counts()

    report = report_file.getvalue()
    report_lines = report.splitlines()
    expected_size = orf1ab_size + nsp12_size + 1
    if len(report_lines) != expected_size:
        assert (len(report_lines), report) == (expected_size, '')

    key_lines = report_lines[nsp12_size + 4396:nsp12_size + 4406]
    key_report = '\n'.join(key_lines)
    assert key_report == expected_text


# noinspection DuplicatedCode
def test_duplicated_sars_base_amino_offset10(sequence_report):
    """ Special case for duplicated base in SARS orf1ab with offset.

    Expect amino sequence AQSFLNRVCG.
    """

    # refname,qcut,rank,count,offset,seq
    aligned_reads = prepare_reads("""\
SARS-CoV-2-seed,15,0,9,10,GCACAATCGTTTTTAAACGGGTTTGCGGTGTAAGTGCAGCCCGTCTTA
""")

    #                                        A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,...,coverage
    expected_text = """\
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,11,4396,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,14,4397,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,17,4398,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,20,4399,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,23,4400,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,26,4401,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,28,4402,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,31,4403,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,34,4404,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,37,4405,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9"""
    sequence_report.projects = ProjectConfig.loadDefault()
    orf1ab_size = len(sequence_report.projects.getReference('SARS-CoV-2-orf1ab'))
    nsp12_size = len(sequence_report.projects.getReference('SARS-CoV-2-nsp12'))

    report_file = StringIO()
    sequence_report.write_amino_header(report_file)
    sequence_report.read(aligned_reads)
    sequence_report.write_amino_counts()

    report = report_file.getvalue()
    report_lines = report.splitlines()
    expected_size = orf1ab_size + nsp12_size + 1
    if len(report_lines) != expected_size:
        assert (len(report_lines), report) == (expected_size, '')

    key_lines = report_lines[nsp12_size + 4396:nsp12_size + 4406]
    key_report = '\n'.join(key_lines)
    assert key_report == expected_text


# noinspection DuplicatedCode
def test_duplicated_sars_base_amino_offset11(sequence_report):
    """ Special case for duplicated base in SARS orf1ab (reading frame 1).

    Expect amino sequence AQSFLNRVCG.
    """

    # refname,qcut,rank,count,offset,seq
    aligned_reads = prepare_reads("""\
SARS-CoV-2-seed,15,0,9,11,GCACAATCGTTTTTAAACGGGTTTGCGGTGTAAGTGCAGCCCGTCTTA
""")

    #                                        A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,...,coverage
    expected_text = """\
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,12,4396,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,15,4397,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,18,4398,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,21,4399,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,24,4400,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,27,4401,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,29,4402,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,32,4403,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,35,4404,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,38,4405,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9"""
    sequence_report.projects = ProjectConfig.loadDefault()
    orf1ab_size = len(sequence_report.projects.getReference('SARS-CoV-2-orf1ab'))
    nsp12_size = len(sequence_report.projects.getReference('SARS-CoV-2-nsp12'))

    report_file = StringIO()
    sequence_report.write_amino_header(report_file)
    sequence_report.read(aligned_reads)
    sequence_report.write_amino_counts()

    report = report_file.getvalue()
    report_lines = report.splitlines()
    expected_size = orf1ab_size + nsp12_size + 1
    if len(report_lines) != expected_size:
        assert (len(report_lines), report) == (expected_size, '')

    key_lines = report_lines[nsp12_size + 4396:nsp12_size + 4406]
    key_report = '\n'.join(key_lines)
    assert key_report == expected_text


def test_duplicated_sars_base_nuc(sequence_report):
    """ Make sure duplicated base in SARS isn't duplicated in nuc.csv. """

    # refname,qcut,rank,count,offset,seq
    aligned_reads = prepare_reads("""\
SARS-CoV-2-seed,15,0,9,10,ACAATCGTTTTTAAACGGGTTTGCGGTGTAAGTGCAGCCCGTCTTACACCG
""")

    #                  A,C,G,T,N,...,coverage
    expected_text = """\
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,21,13198,0,0,0,9,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,22,13199,0,0,0,9,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,23,13200,9,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,24,13201,9,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,25,13202,9,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,26,13203,0,9,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,27,13205,0,0,9,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,28,13206,0,0,9,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,29,13207,0,0,9,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,30,13208,0,0,0,9,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-orf1ab,15,31,13209,0,0,0,9,0,0,0,0,0,9"""
    sequence_report.projects = ProjectConfig.loadDefault()
    orf1ab_size = len(sequence_report.projects.getReference('SARS-CoV-2-orf1ab'))
    nsp12_size = len(sequence_report.projects.getReference('SARS-CoV-2-nsp12'))

    report_file = StringIO()
    sequence_report.write_nuc_header(report_file)
    sequence_report.read(aligned_reads)
    sequence_report.write_nuc_counts()

    report = report_file.getvalue()
    report_lines = report.splitlines()
    header_size = 1
    skipped_rows = 2
    expected_size = (orf1ab_size + nsp12_size)*3 + header_size - skipped_rows
    if len(report_lines) != expected_size:
        assert (len(report_lines), report) == (expected_size, '')

    key_lines = report_lines[nsp12_size*3 + 13197:nsp12_size*3 + 13208]
    key_report = '\n'.join(key_lines)
    assert key_report == expected_text


# noinspection DuplicatedCode
def test_contig_coverage_report_huge_gap(projects):
    """ A gap so big that Gotoh can't bridge it, but minimap2 can. """
    sequence_report = SequenceReport(InsertionWriter(StringIO()),
                                     projects,
                                     conseq_mixture_cutoffs=[0.1])
    ref = projects.getReference('HIV1-B-FR-K03455-seed')
    seq = ref[100:150] + ref[1000:1050]
    expected_positions = list(range(101, 151)) + list(range(1001, 1051))
    remap_conseq_csv = StringIO(f"""\
region,sequence
HIV1-B-FR-K03455-seed,{seq}
""")
    # refname,qcut,rank,count,offset,seq
    aligned_reads1 = prepare_reads(f"""\
HIV1-B-FR-K03455-seed,15,0,4,0,{seq}
""")

    report_file = StringIO()
    sequence_report.read_remap_conseqs(remap_conseq_csv)
    sequence_report.write_amino_header(StringIO())
    sequence_report.write_genome_coverage_header(report_file)
    sequence_report.read(aligned_reads1)
    sequence_report.write_genome_coverage_counts()
    sequence_report.write_amino_counts()

    report_file.seek(0)
    covered_positions = [int(row['refseq_nuc_pos'])
                         for row in DictReader(report_file)
                         if row['refseq_nuc_pos']]
    assert covered_positions == expected_positions


# noinspection DuplicatedCode
def test_contig_coverage_report_past_reference_end(projects, sequence_report):
    hxb2_name = 'HIV1-B-FR-K03455-seed'
    ref = projects.getReference(hxb2_name)
    landmarks = dict(coordinates=hxb2_name)
    assert len(ref) == 9719
    seq = ref[-100:] + 'CGTAC'
    seed_nucs = [('C', SeedNucleotide(Counter({'C': 1})))] * len(seq)
    expected_tail = """\
1-my-contig,HIV1-B-FR-K03455-seed,99,9718,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,100,9719,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,101,9720,0,1,U
1-my-contig,HIV1-B-FR-K03455-seed,102,9721,0,1,U
1-my-contig,HIV1-B-FR-K03455-seed,103,9722,0,1,U
1-my-contig,HIV1-B-FR-K03455-seed,104,9723,0,1,U
1-my-contig,HIV1-B-FR-K03455-seed,105,9724,0,1,U
"""

    report_file = StringIO()
    sequence_report.write_genome_coverage_header(report_file)
    SequenceReport.write_sequence_coverage_counts(
        projects,
        sequence_report.genome_coverage_writer,
        '1-my-contig',
        landmarks,
        seq,
        seed_nucs=seed_nucs)

    report_text = report_file.getvalue()
    tail = report_text[-len(expected_tail):]
    assert tail == expected_tail


# noinspection DuplicatedCode
def test_contig_coverage_report_past_reference_start(projects, sequence_report):
    hxb2_name = 'HIV1-B-FR-K03455-seed'
    ref = projects.getReference(hxb2_name)
    landmarks = dict(coordinates=hxb2_name)
    seq = 'CGTAC' + ref[:100]
    seed_nucs = [('C', SeedNucleotide(Counter({'C': 1})))] * len(seq)
    # link is (M)apped, (U)nmapped, or (I)nserted
    expected_head = """\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage,link
1-my-contig,HIV1-B-FR-K03455-seed,1,-4,0,1,U
1-my-contig,HIV1-B-FR-K03455-seed,2,-3,0,1,U
1-my-contig,HIV1-B-FR-K03455-seed,3,-2,0,1,U
1-my-contig,HIV1-B-FR-K03455-seed,4,-1,0,1,U
1-my-contig,HIV1-B-FR-K03455-seed,5,0,0,1,U
1-my-contig,HIV1-B-FR-K03455-seed,6,1,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,7,2,0,1,M
"""

    report_file = StringIO()
    sequence_report.write_genome_coverage_header(report_file)
    SequenceReport.write_sequence_coverage_counts(
        projects,
        sequence_report.genome_coverage_writer,
        '1-my-contig',
        landmarks,
        seq,
        seed_nucs=seed_nucs)

    report_text = report_file.getvalue()
    head = report_text[:len(expected_head)]
    assert head == expected_head


# noinspection DuplicatedCode
def test_contig_coverage_report_offset_reads(projects, sequence_report):
    hxb2_name = 'HIV1-B-FR-K03455-seed'
    ref = projects.getReference(hxb2_name)
    landmarks = dict(coordinates=hxb2_name)
    seq = ref[50:150]
    seed_nucs = [('C', SeedNucleotide(Counter({'C': 1})))] * len(seq)
    expected_head = """\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage,link
1-my-contig,HIV1-B-FR-K03455-seed,51,51,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,52,52,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,53,53,0,1,M
"""

    report_file = StringIO()
    sequence_report.write_genome_coverage_header(report_file)
    SequenceReport.write_sequence_coverage_counts(
        projects,
        sequence_report.genome_coverage_writer,
        '1-my-contig',
        landmarks,
        seq,
        consensus_offset=50,
        seed_nucs=seed_nucs)

    report_text = report_file.getvalue()
    head = report_text[:len(expected_head)]
    assert head == expected_head


def test_contig_coverage_report_for_partial_contig(sequence_report):
    """ Contig coverage is reported for partial contigs.

    Reference columns are left blank, though, because they're not aligned.
    See blast.csv for best guess at alignment.
    """
    # refname,qcut,rank,count,offset,seq
    aligned_reads1 = prepare_reads("1-R1-seed-partial,15,0,5,0,CCCCCC")
    contigs_csv = StringIO("""\
ref,match,group_ref,contig
R1-seed,1,R1-seed,CCCCCC
""")

    expected_text = """\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage,link
1-R1-seed-partial,,1,1,0,5,U
1-R1-seed-partial,,2,2,0,5,U
1-R1-seed-partial,,3,3,0,5,U
1-R1-seed-partial,,4,4,0,5,U
1-R1-seed-partial,,5,5,0,5,U
1-R1-seed-partial,,6,6,0,5,U
"""

    report_file = StringIO()
    sequence_report.write_amino_header(StringIO())
    sequence_report.write_amino_detail_header(StringIO())
    sequence_report.read_contigs(contigs_csv)
    sequence_report.write_genome_coverage_header(report_file)
    sequence_report.read(aligned_reads1)
    sequence_report.write_genome_coverage_counts()

    assert report_file.getvalue() == expected_text


def test_contig_coverage_report_for_reversed_contig(sequence_report):
    """ Contig coverage is reported for reversed contigs.

    Reference columns are left blank, though, because they're not aligned.
    See blast.csv for best guess at alignment.
    """
    # refname,qcut,rank,count,offset,seq
    aligned_reads1 = prepare_reads("1-R1-seed-reversed,15,0,5,0,CCCCCC")
    contigs_csv = StringIO("""\
ref,match,group_ref,contig
R1-seed,1,R1-seed,CCCCCC
""")

    expected_text = """\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage,link
1-R1-seed-reversed,,1,1,0,5,U
1-R1-seed-reversed,,2,2,0,5,U
1-R1-seed-reversed,,3,3,0,5,U
1-R1-seed-reversed,,4,4,0,5,U
1-R1-seed-reversed,,5,5,0,5,U
1-R1-seed-reversed,,6,6,0,5,U
"""

    report_file = StringIO()
    sequence_report.read_contigs(contigs_csv)
    sequence_report.write_genome_coverage_header(report_file)
    sequence_report.read(aligned_reads1)
    sequence_report.write_genome_coverage_counts()
    sequence_report.write_amino_detail_counts()

    assert report_file.getvalue() == expected_text


def test_contig_coverage_report_merged_contigs(sequence_report):
    """ Assemble counts from three contigs to two references.

    Reads:
    Contig 1_3-R1 AAATTTAGG -> KFR
    Contig 2-R2 GGCCCG -> GP

    Contigs:
    Contig 1-R1 AAATTT -> KF
    Contig 2-R2 GGCCCG -> GP
    Contig 3-R1 TTTAGG -> FR

    Contig 1 and 3 have been combined into R1 with KFR.
    """
    # refname,qcut,rank,count,offset,seq
    aligned_reads1 = prepare_reads("1_3-R1-seed,15,0,5,0,AAATTT\n"
                                   "1_3-R1-seed,15,0,2,3,TTTAGG")
    aligned_reads2 = prepare_reads("2-R2-seed,15,0,4,0,GGCCCG")
    contigs_csv = StringIO("""\
ref,match,group_ref,contig
R1-seed,1,R1-seed,AAATTT
R2-seed,1,R2-seed,GGCCCG
R1-seed,1,R1-seed,TTTAGG
""")

    expected_text = """\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage,link
1_3-R1-seed,R1-seed,1,1,0,5,U
1_3-R1-seed,R1-seed,2,2,0,5,U
1_3-R1-seed,R1-seed,3,3,0,5,U
1_3-R1-seed,R1-seed,4,4,0,7,U
1_3-R1-seed,R1-seed,5,5,0,7,U
1_3-R1-seed,R1-seed,6,6,0,7,U
1_3-R1-seed,R1-seed,7,7,0,2,U
1_3-R1-seed,R1-seed,8,8,0,2,U
1_3-R1-seed,R1-seed,9,9,0,2,U
contig-1-R1-seed,R1-seed,1,1,,,U
contig-1-R1-seed,R1-seed,2,2,,,U
contig-1-R1-seed,R1-seed,3,3,,,U
contig-1-R1-seed,R1-seed,4,4,,,U
contig-1-R1-seed,R1-seed,5,5,,,U
contig-1-R1-seed,R1-seed,6,6,,,U
contig-3-R1-seed,R1-seed,1,1,,,U
contig-3-R1-seed,R1-seed,2,2,,,U
contig-3-R1-seed,R1-seed,3,3,,,U
contig-3-R1-seed,R1-seed,4,4,,,U
contig-3-R1-seed,R1-seed,5,5,,,U
contig-3-R1-seed,R1-seed,6,6,,,U
2-R2-seed,R2-seed,1,1,0,4,U
2-R2-seed,R2-seed,2,2,0,4,U
2-R2-seed,R2-seed,3,3,0,4,U
2-R2-seed,R2-seed,4,4,0,4,U
2-R2-seed,R2-seed,5,5,0,4,U
2-R2-seed,R2-seed,6,6,0,4,U
contig-2-R2-seed,R2-seed,1,1,,,U
contig-2-R2-seed,R2-seed,2,2,,,U
contig-2-R2-seed,R2-seed,3,3,,,U
contig-2-R2-seed,R2-seed,4,4,,,U
contig-2-R2-seed,R2-seed,5,5,,,U
contig-2-R2-seed,R2-seed,6,6,,,U
"""

    report_file = StringIO()
    sequence_report.read_contigs(contigs_csv)
    sequence_report.write_amino_header(StringIO())
    sequence_report.write_amino_detail_header(StringIO())
    sequence_report.write_genome_coverage_header(report_file)
    sequence_report.read(aligned_reads1)
    sequence_report.write_genome_coverage_counts()
    sequence_report.write_amino_counts()
    sequence_report.read(aligned_reads2)
    sequence_report.write_genome_coverage_counts()
    sequence_report.write_amino_counts()

    assert report_file.getvalue() == expected_text


# noinspection DuplicatedCode
def test_write_sequence_coverage_counts_without_coverage(projects,
                                                         sequence_report):
    hxb2_name = 'HIV1-B-FR-K03455-seed'
    ref = projects.getReference(hxb2_name)
    landmarks = dict(coordinates=hxb2_name)
    seq = ref[100:150] + ref[1000:1050]
    expected_positions = list(range(101, 151)) + list(range(1001, 1051))

    report_file = StringIO()
    sequence_report.write_genome_coverage_header(report_file)
    SequenceReport.write_sequence_coverage_counts(
        projects,
        sequence_report.genome_coverage_writer,
        '1-my-contig',
        landmarks,
        seq)

    report_file.seek(0)
    covered_positions = [int(row['refseq_nuc_pos'])
                         for row in DictReader(report_file)
                         if row['refseq_nuc_pos']]
    assert covered_positions == expected_positions


# noinspection DuplicatedCode
def test_write_sequence_coverage_counts_with_coverage(projects,
                                                      sequence_report):
    hxb2_name = 'HIV1-B-FR-K03455-seed'
    ref = projects.getReference(hxb2_name)
    landmarks = dict(coordinates=hxb2_name)
    seq = ref[100:150] + ref[1000:1050]
    seed_nucs = [('C', SeedNucleotide(Counter({'C': 1})))] * 100
    seed_nucs[2] = ('G', SeedNucleotide(Counter({'G': 4})))
    seed_nucs[98] = ('T', SeedNucleotide(Counter({'T': 5})))
    expected_head = """\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage,link
1-my-contig,HIV1-B-FR-K03455-seed,1,101,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,2,102,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,3,103,0,4,M
1-my-contig,HIV1-B-FR-K03455-seed,4,104,0,1,M
"""
    expected_tail = """\
1-my-contig,HIV1-B-FR-K03455-seed,98,1048,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,99,1049,0,5,M
1-my-contig,HIV1-B-FR-K03455-seed,100,1050,0,1,M
"""

    report_file = StringIO()
    sequence_report.write_genome_coverage_header(report_file)
    SequenceReport.write_sequence_coverage_counts(
        projects,
        sequence_report.genome_coverage_writer,
        '1-my-contig',
        landmarks,
        seq,
        seed_nucs=seed_nucs)

    report_text = report_file.getvalue()
    head = report_text[:len(expected_head)]
    tail = report_text[-len(expected_tail):]
    assert head == expected_head
    assert tail == expected_tail


# noinspection DuplicatedCode
def test_write_sequence_coverage_counts_with_deletion(projects, sequence_report):
    hxb2_name = 'HIV1-B-FR-K03455-seed'
    ref = projects.getReference(hxb2_name)
    landmarks = dict(coordinates=hxb2_name)
    seq = ref[100:110] + ref[115:160]
    seed_nucs = [('C', SeedNucleotide(Counter({'C': 1})))] * len(seq)
    seed_nucs[12] = ('G', SeedNucleotide(Counter({'G': 4})))
    expected_head = """\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage,link
1-my-contig,HIV1-B-FR-K03455-seed,1,101,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,2,102,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,3,103,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,4,104,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,5,105,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,6,106,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,7,107,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,8,108,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,9,109,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,10,110,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,11,116,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,12,117,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,13,118,0,4,M
1-my-contig,HIV1-B-FR-K03455-seed,14,119,0,1,M
"""

    report_file = StringIO()
    sequence_report.write_genome_coverage_header(report_file)
    SequenceReport.write_sequence_coverage_counts(
        projects,
        sequence_report.genome_coverage_writer,
        '1-my-contig',
        landmarks,
        seq,
        seed_nucs=seed_nucs)

    report_text = report_file.getvalue()
    head = report_text[:len(expected_head)]
    assert head == expected_head


# noinspection DuplicatedCode
def test_write_sequence_coverage_counts_with_some_deletions(projects,
                                                            sequence_report):
    """ Some reads had deletions at a position. """
    hxb2_name = 'HIV1-B-FR-K03455-seed'
    ref = projects.getReference(hxb2_name)
    landmarks = dict(coordinates=hxb2_name)
    seq = ref[100:150]
    seed_nucs = [('C', SeedNucleotide(Counter({'C': 1})))] * len(seq)
    seed_nucs[5] = ('G', SeedNucleotide(Counter({'G': 4, '-': 2})))
    expected_head = """\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage,link
1-my-contig,HIV1-B-FR-K03455-seed,1,101,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,2,102,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,3,103,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,4,104,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,5,105,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,6,106,2,6,M
1-my-contig,HIV1-B-FR-K03455-seed,7,107,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,8,108,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,9,109,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,10,110,0,1,M
"""

    report_file = StringIO()
    sequence_report.write_genome_coverage_header(report_file)
    SequenceReport.write_sequence_coverage_counts(
        projects,
        sequence_report.genome_coverage_writer,
        '1-my-contig',
        landmarks,
        seq,
        seed_nucs=seed_nucs)

    report_text = report_file.getvalue()
    head = report_text[:len(expected_head)]
    assert head == expected_head


# noinspection DuplicatedCode
def test_write_sequence_coverage_counts_with_insert(projects, sequence_report):
    hxb2_name = 'HIV1-B-FR-K03455-seed'
    ref = projects.getReference(hxb2_name)
    landmarks = dict(coordinates=hxb2_name)
    seq = ref[100:110] + 'ACTGA' + ref[110:160]
    seed_nucs = [('C', SeedNucleotide(Counter({'C': 1})))] * len(seq)
    seed_nucs[12] = ('T', SeedNucleotide(Counter({'T': 4})))
    expected_head = """\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage,link
1-my-contig,HIV1-B-FR-K03455-seed,1,101,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,2,102,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,3,103,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,4,104,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,5,105,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,6,106,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,7,107,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,8,108,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,9,109,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,10,110,0,1,M
1-my-contig,HIV1-B-FR-K03455-seed,11,,0,1,I
1-my-contig,HIV1-B-FR-K03455-seed,12,,0,1,I
1-my-contig,HIV1-B-FR-K03455-seed,13,,0,4,I
1-my-contig,HIV1-B-FR-K03455-seed,14,,0,1,I
1-my-contig,HIV1-B-FR-K03455-seed,15,,0,1,I
1-my-contig,HIV1-B-FR-K03455-seed,16,111,0,1,M
"""

    report_file = StringIO()
    sequence_report.write_genome_coverage_header(report_file)
    SequenceReport.write_sequence_coverage_counts(
        projects,
        sequence_report.genome_coverage_writer,
        '1-my-contig',
        landmarks,
        seq,
        seed_nucs=seed_nucs)

    report_text = report_file.getvalue()
    head = report_text[:len(expected_head)]
    assert head == expected_head


# noinspection DuplicatedCode
def test_write_sequence_coverage_counts_with_unaligned_middle(projects,
                                                              sequence_report):
    """ The middle 100 bases are from a different reference.

    They get reported with query positions, but no reference positions.
    """
    hxb2_name = 'HIV1-B-FR-K03455-seed'
    ref = projects.getReference(hxb2_name)
    hcv_ref = projects.getReference('HCV-1a')
    landmarks = dict(coordinates=hxb2_name)
    seq = ref[:100] + hcv_ref[1000:1100] + ref[1000:1100]
    expected_ref_positions = (list(range(1, 101)) +
                              list(range(501, 601)) +
                              list(range(1001, 1101)))
    expected_query_positions = list(range(1, 301))

    report_file = StringIO()
    sequence_report.write_genome_coverage_header(report_file)
    SequenceReport.write_sequence_coverage_counts(
        projects,
        sequence_report.genome_coverage_writer,
        '1-my-contig',
        landmarks,
        seq)

    report_file.seek(0)
    rows = list(DictReader(report_file))
    ref_positions = [int(row['refseq_nuc_pos'])
                     for row in rows
                     if row['refseq_nuc_pos']]
    query_positions = [int(row['query_nuc_pos'])
                       for row in rows
                       if row['query_nuc_pos']]
    assert query_positions == expected_query_positions
    assert ref_positions == expected_ref_positions


# noinspection DuplicatedCode
def test_write_sequence_coverage_counts_with_double_mapped_edges(
        projects,
        sequence_report):
    """ There's a large deletion, but the junction maps to both edges.

    The first 8 bases of the second section (8188-8195) also map to the 7 bases
    after the first section (2909-2915). We arbitrarily display them with the
    first section.
    """
    hxb2_name = 'HIV1-B-FR-K03455-seed'
    ref = projects.getReference(hxb2_name)
    landmarks = dict(coordinates=hxb2_name)
    seq = ref[2858:2908] + ref[8187:8237]
    expected_ref_positions = (list(range(2859, 2916)) + list(range(8196, 8238)))
    expected_query_positions = list(range(1, len(seq)+1))

    report_file = StringIO()
    sequence_report.write_genome_coverage_header(report_file)
    SequenceReport.write_sequence_coverage_counts(
        projects,
        sequence_report.genome_coverage_writer,
        '1-my-contig',
        landmarks,
        seq)

    report_file.seek(0)
    rows = list(DictReader(report_file))
    ref_positions = [int(row['refseq_nuc_pos'])
                     for row in rows
                     if row['refseq_nuc_pos']]
    query_positions = [int(row['query_nuc_pos'])
                       for row in rows
                       if row['query_nuc_pos']]
    assert query_positions == expected_query_positions
    assert ref_positions == expected_ref_positions
