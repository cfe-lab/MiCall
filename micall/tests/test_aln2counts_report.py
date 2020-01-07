import csv
import sys
from io import StringIO

import pytest

from micall.core import project_config
from micall.core.aln2counts import InsertionWriter, SequenceReport


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
