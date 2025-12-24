import csv
import logging
import sys
from collections import Counter
from csv import DictReader
from io import StringIO

import pytest
import yaml
from mappy import revcomp

from micall.core import project_config
from micall.core.aln2counts import InsertionWriter, SequenceReport, combine_region_nucleotides, \
    combine_region_insertions
from micall.core.project_config import ProjectConfig
from micall.utils.report_amino import SeedNucleotide, ReportNucleotide

# noinspection PyUnresolvedReferences
from micall.tests.test_remap import load_projects


def prepare_reads(aligned_reads_text):
    full_text = "refname,qcut,rank,count,offset,seq\n" + aligned_reads_text
    dummy_file = StringIO(full_text)
    return csv.DictReader(dummy_file)


@pytest.fixture
def sequence_report():
    return create_sequence_report()


@pytest.fixture
def default_sequence_report():
    """ Sequence report with stubbed files, but all default projects. """
    conseq_mixture_cutoffs = [0.1]
    return StubbedSequenceReport(InsertionWriter(StringIO()),
                                 ProjectConfig.loadDefault(),
                                 conseq_mixture_cutoffs)


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
"R1A": {
  "max_variants": 10,
  "regions": [
    {
      "coordinate_region": "R1A",
      "seed_region_names": ["R1A-seed"]
    },
    {
      "coordinate_region": "R1A_second",
      "seed_region_names": ["R1A-seed"]
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
"R1A-seed": {
  "is_nucleotide": true,
  "reference": [
    "AAATTTAGGTAGCCGAGACCTCAGGTCACTCTTTGGTAG"
  ]
},
"R1A": {
  "is_nucleotide": false,
  "reference": [
    "KFR*"
  ]
},
"R1A_second": {
  "is_nucleotide": false,
  "reference": [
    "PRPQVTLW*"
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
    "ATGGCAAACTCAATCAATGGG"
  ]
},
"R4": {
  "is_nucleotide": false,
  "reference": [
    "SING"
  ]
},
"R5-seed": {
  "comment": "Coord has G that's not in seed.",
  "is_nucleotide": true,
  "reference": [
    "AAATTTCCGAGACCTCAGGTCACTCTTTGG"
  ]
},
"R5": {
  "is_nucleotide": false,
  "reference": [
    "KFGPRPQVTLW"
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
    # Extra 3 positions for stop codon to get dropped.
    - {name: R1, start: 1, end: 12, colour: steelblue}
- seed_pattern: R1A-.*
  coordinates: R1A-seed
  landmarks:
    - {name: R1A, start: 1, end: 12, colour: steelblue}
    - {name: R1A_second, start: 13, end: 39, colour: pink}
- seed_pattern: R2-.*
  coordinates: R2-seed
  landmarks:
    # Extra 3 positions for stop codon to get dropped.
    - {name: R2, start: 1, end: 18, colour: steelblue}
- seed_pattern: R3-.*
  coordinates: R3-seed
  landmarks:
    # Extra 3 positions for stop codons to get dropped, one codon overlaps.
    - {name: R3, start: 1, end: 15, colour: lightblue}
    - {name: R3-extra, start: 13, end: 27, colour: steelblue}
- seed_pattern: R4-.*
  coordinates: R4-seed
  landmarks:
    # Extra 3 positions for stop codons to get dropped.
    - {name: R4, start: 10, end: 24}
- seed_pattern: R5-.*
  coordinates: R5-seed
  landmarks:
    # Extra 3 positions for stop codons to get dropped.
    - {name: R5, start: 1, end: 33}
- seed_pattern: R7-.*
  coordinates: R7-seed
  landmarks:
    # Extra 3 positions for stop codons to get dropped.
    - {name: R7a, start: 1, end: 12}
    - {name: R7b, start: 16, end: 27}
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
        "GGCCCCAAATTTAGGGAGCAC"
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
    sequence_report.landmarks = yaml.safe_load("""\
- seed_pattern: R1-seed
  coordinates: R1-seed
  landmarks:
    # Extra 3 nucleotides at end, because stop codons will get dropped.
    - {name: R1, start: 7, end: 18, frame: 0}
    - {name: R1-expanded, start: 1, end: 24, frame: 0}
""")
    return sequence_report


class StubbedSequenceReport(SequenceReport):
    def __init__(self, *args, **kwargs):
        SequenceReport.__init__(self, *args, **kwargs)
        self.overrides = {}

    def _pair_align(self, reference, query, *args, **kwargs):
        override = self.overrides.get((reference, query))
        return (override
                if override is not None
                else SequenceReport._pair_align(self,
                                                reference,
                                                query,
                                                *args,
                                                **kwargs))

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
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,genome.pos,\
A,C,G,T,N,del,ins,clip,v3_overlap,coverage,exact_coverage
R1-seed,R1,15,1,1,1,9,0,0,0,0,0,0,0,0,9,
R1-seed,R1,15,2,2,2,9,0,0,0,0,0,0,0,0,9,18
R1-seed,R1,15,3,3,3,9,0,0,0,0,0,0,0,0,9,18
R1-seed,R1,15,4,4,4,0,0,0,9,0,0,0,0,0,9,18
R1-seed,R1,15,5,5,5,0,0,0,9,0,0,0,0,0,9,18
R1-seed,R1,15,6,6,6,0,0,0,9,0,0,0,0,0,9,
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
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,genome.pos,\
A,C,G,T,N,del,ins,clip,v3_overlap,coverage,exact_coverage
R1-seed,R1,15,1,1,7,5,0,0,0,0,0,0,0,0,5,
R1-seed,R1,15,2,2,8,5,0,0,0,0,0,0,0,0,5,
R1-seed,R1,15,3,3,9,5,0,0,0,0,0,0,0,0,5,
R1-seed,R1,15,,4,10,0,0,0,7,0,0,0,0,0,7,
R1-seed,R1,15,,5,11,0,0,0,7,0,0,0,0,0,7,
R1-seed,R1,15,,6,12,0,0,0,7,0,0,0,0,0,7,
R1-seed,R1,15,4,7,13,2,0,0,0,0,0,0,0,0,2,
R1-seed,R1,15,5,8,14,0,0,2,0,0,0,0,0,0,2,
R1-seed,R1,15,6,9,15,0,0,2,0,0,0,0,0,0,2,
R1-seed,R1-expanded,15,1,7,7,5,0,0,0,0,0,0,0,0,5,
R1-seed,R1-expanded,15,2,8,8,5,0,0,0,0,0,0,0,0,5,
R1-seed,R1-expanded,15,3,9,9,5,0,0,0,0,0,0,0,0,5,
R1-seed,R1-expanded,15,,10,10,0,0,0,7,0,0,0,0,0,7,
R1-seed,R1-expanded,15,,11,11,0,0,0,7,0,0,0,0,0,7,
R1-seed,R1-expanded,15,,12,12,0,0,0,7,0,0,0,0,0,7,
R1-seed,R1-expanded,15,4,13,13,2,0,0,0,0,0,0,0,0,2,
R1-seed,R1-expanded,15,5,14,14,0,0,2,0,0,0,0,0,0,2,
R1-seed,R1-expanded,15,6,15,15,0,0,2,0,0,0,0,0,0,2,
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
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,genome.pos,\
A,C,G,T,N,del,ins,clip,v3_overlap,coverage,exact_coverage
R1-seed,R1,15,1,1,7,5,0,0,0,0,0,0,0,0,5,
R1-seed,R1,15,2,2,8,5,0,0,0,0,0,0,0,0,5,
R1-seed,R1,15,3,3,9,5,0,0,0,0,0,0,0,0,5,
R1-seed,R1,15,4,4,10,0,0,0,5,0,0,0,0,0,5,
R1-seed,R1,15,5,5,11,0,0,0,5,0,0,0,0,0,5,
R1-seed,R1,15,6,6,12,0,0,0,5,0,0,0,0,0,5,
R1-seed,R1,15,7,7,13,5,0,0,0,0,0,0,0,0,5,
R1-seed,R1,15,8,8,14,0,0,5,0,0,0,0,0,0,5,
R1-seed,R1,15,9,9,15,0,0,5,0,0,0,0,0,0,5,5
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
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,genome.pos,\
A,C,G,T,N,del,ins,clip,v3_overlap,coverage,exact_coverage
R1-seed,R1-expanded,15,1,7,7,5,0,0,0,0,0,0,0,0,5,
R1-seed,R1-expanded,15,2,8,8,5,0,0,0,0,0,0,0,0,5,
R1-seed,R1-expanded,15,3,9,9,5,0,0,0,0,0,0,0,0,5,
R1-seed,R1-expanded,15,4,10,10,0,0,0,5,0,0,0,0,0,5,
R1-seed,R1-expanded,15,5,11,11,0,0,0,5,0,0,0,0,0,5,
R1-seed,R1-expanded,15,6,12,12,0,0,0,5,0,0,0,0,0,5,
R1-seed,R1-expanded,15,7,13,13,5,0,0,0,0,0,0,0,0,5,
R1-seed,R1-expanded,15,8,14,14,0,0,5,0,0,0,0,0,0,5,
R1-seed,R1-expanded,15,9,15,15,0,0,5,0,0,0,0,0,0,5,5
"""

    report = sequence_report_overlapping_regions
    report_file = StringIO()
    report.write_nuc_header(report_file)
    report.read(aligned_reads, included_regions={'R1-expanded'})
    report.write_nuc_counts()

    assert report_file.getvalue() == expected_text


# noinspection DuplicatedCode
def test_duplicated_sars_base_amino(default_sequence_report):
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
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,1,4396,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,4,4397,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,7,4398,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,10,4399,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,13,4400,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,16,4401,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,18,4402,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,21,4403,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,24,4404,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,27,4405,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9"""

    report_file = StringIO()
    default_sequence_report.write_amino_header(report_file)
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_amino_counts()

    report = report_file.getvalue()
    report_lines = report.splitlines()
    expected_size = 45
    if len(report_lines) != expected_size:
        assert (len(report_lines), report) == (expected_size, '')

    key_lines = report_lines[1:11]
    key_report = '\n'.join(key_lines)
    assert key_report == expected_text


# noinspection DuplicatedCode
def test_duplicated_sars_base_amino_offset10(default_sequence_report):
    """ Special case for duplicated base in SARS orf1ab with offset.

    Expect amino sequence AQSFLNRVCG.
    """

    # refname,qcut,rank,count,offset,seq
    aligned_reads = prepare_reads("""\
SARS-CoV-2-seed,15,0,9,10,GCACAATCGTTTTTAAACGGGTTTGCGGTGTAAGTGCAGCCCGTCTTA
""")

    #                                        A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,...,coverage
    expected_text = """\
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,11,4396,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,14,4397,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,17,4398,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,20,4399,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,23,4400,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,26,4401,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,28,4402,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,31,4403,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,34,4404,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,37,4405,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9"""
    report_file = StringIO()
    default_sequence_report.write_amino_header(report_file)
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_amino_counts()

    report = report_file.getvalue()
    report_lines = report.splitlines()
    expected_size = 43
    if len(report_lines) != expected_size:
        assert (len(report_lines), report) == (expected_size, '')

    key_lines = report_lines[1:11]
    key_report = '\n'.join(key_lines)
    assert key_report == expected_text


# noinspection DuplicatedCode
def test_duplicated_sars_base_amino_offset11(default_sequence_report):
    """ Special case for duplicated base in SARS orf1ab (reading frame 1).

    Expect amino sequence AQSFLNRVCG.
    """

    # refname,qcut,rank,count,offset,seq
    aligned_reads = prepare_reads("""\
SARS-CoV-2-seed,15,0,9,11,GCACAATCGTTTTTAAACGGGTTTGCGGTGTAAGTGCAGCCCGTCTTA
""")

    #                                        A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,...,coverage
    expected_text = """\
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,12,4396,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,15,4397,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,18,4398,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,21,4399,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,24,4400,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,27,4401,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,29,4402,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,32,4403,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,35,4404,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,38,4405,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9"""
    report_file = StringIO()
    default_sequence_report.write_amino_header(report_file)
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_amino_counts()

    report = report_file.getvalue()
    report_lines = report.splitlines()
    expected_size = 43
    if len(report_lines) != expected_size:
        assert (len(report_lines), report) == (expected_size, '')

    key_lines = report_lines[1:11]
    key_report = '\n'.join(key_lines)
    assert key_report == expected_text


# noinspection DuplicatedCode
def test_duplicated_sars_base_nuc(default_sequence_report):
    """ Make sure duplicated base in SARS isn't duplicated in nuc.csv.

    Duplicated base is position 13468 in the whole genome, or 13203 in ORF1a.
    """

    # refname,qcut,rank,count,offset,seq
    aligned_reads = prepare_reads("""\
SARS-CoV-2-seed,15,0,9,10,ACAATCGTTTTTAAACGGGTTTGCGGTGTAAGTGCAGCCCGTCTTACACCG
""")
    #                                    ^ Duplicated base

    #                                               A,C,G,T,N,...,coverage
    expected_section = """\
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,21,13198,13463,0,0,0,9,0,0,0,0,0,9,
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,22,13199,13464,0,0,0,9,0,0,0,0,0,9,
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,23,13200,13465,9,0,0,0,0,0,0,0,0,9,
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,24,13201,13466,9,0,0,0,0,0,0,0,0,9,
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,25,13202,13467,9,0,0,0,0,0,0,0,0,9,
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,26,13203,13468,0,9,0,0,0,0,0,0,0,9,
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,27,13204,13469,0,0,9,0,0,0,0,0,0,9,
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,28,13205,13470,0,0,9,0,0,0,0,0,0,9,
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,29,13206,13471,0,0,9,0,0,0,0,0,0,9,
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,30,13207,13472,0,0,0,9,0,0,0,0,0,9,
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,31,13208,13473,0,0,0,9,0,0,0,0,0,9,"""

    report_file = StringIO()
    default_sequence_report.write_nuc_header(report_file)
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_nuc_counts()

    report = report_file.getvalue()
    report_lines = report.splitlines()
    expected_size = 131
    if len(report_lines) != expected_size:
        assert (len(report_lines), report) == (expected_size, '')

    key_lines = report_lines[11:22]
    key_report = '\n'.join(key_lines)
    assert key_report == expected_section


# noinspection DuplicatedCode
def test_nsp11_duplicated_sars_base(default_sequence_report):
    """ Make sure that the duplicated position is NOT duplicated in nsp11"""

    # refname,qcut,rank,count,offset,seq
    aligned_reads = prepare_reads("""\
SARS-CoV-2-seed,15,0,9,0,GCACAATCGTTTTTAAACGGGTTTGCGGTGTAAGTGCAGCCCGTCTTACAC
""")
    # Repeat is here:                     ^
    # Expected amino sequence: AQSFLNGFAV

    #                                   A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,...,coverage
    expected_text = """\
SARS-CoV-2-seed,SARS-CoV-2-nsp11,15,1,4,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-nsp11,15,4,5,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-nsp11,15,7,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-nsp11,15,10,7,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-nsp11,15,13,8,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-nsp11,15,16,9,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-nsp11,15,19,10,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-nsp11,15,22,11,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-nsp11,15,25,12,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
SARS-CoV-2-seed,SARS-CoV-2-nsp11,15,28,13,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,9"""

    report_file = StringIO()
    default_sequence_report.write_amino_header(report_file)
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_amino_counts()

    report = report_file.getvalue()
    report_lines = report.splitlines()
    expected_size = 45
    if len(report_lines) != expected_size:
        assert (len(report_lines), report) == (expected_size, '')

    key_lines = report_lines[18:28]
    key_report = '\n'.join(key_lines)
    assert key_report == expected_text


# noinspection DuplicatedCode
def test_duplicated_sars_base_last_region_nuc(default_sequence_report):
    """ Make sure that the last nucleotide of the region with the duplicated position is included in the report"""

    # refname,qcut,rank,count,offset,seq
    aligned_reads = prepare_reads("""\
SARS-CoV-2-seed,15,0,9,10,ACAATCGTTTTTAAACGGGTTTGCGGTGTAAGTGCAGCCCGTCTTACACCG
""")
    #                                ^ Duplicated base

    #                                               A,C,G,T,N,...,coverage
    expected_section = """\
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,34,13211,13476,0,9,0,0,0,0,0,0,0,9,
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,35,13212,13477,0,0,9,0,0,0,0,0,0,9,
SARS-CoV-2-seed,SARS-CoV-2-ORF1ab,15,36,13213,13478,0,0,9,0,0,0,0,0,0,9,"""

    report_file = StringIO()
    default_sequence_report.write_nuc_header(report_file)
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_nuc_counts()

    report = report_file.getvalue()
    report_lines = report.splitlines()
    expected_size = 131
    if len(report_lines) != expected_size:
        assert (len(report_lines), report) == (expected_size, '')

    key_lines = report_lines[24:27]
    key_report = '\n'.join(key_lines)
    assert key_report == expected_section


# noinspection DuplicatedCode
def test_duplicated_sars_base_last_contig_nuc(default_sequence_report):
    """ Make sure that the last nucleotide of the contig is included in the report,
    if the contig ends in a section with the duplicated nucleotide"""

    # refname,qcut,rank,count,offset,seq
    aligned_reads = prepare_reads("""\
SARS-CoV-2-seed,15,0,9,10,ACAATCGTTTTTAAACGGGTTTGCGGTGTAAGTGCAGCCCGTCTTACACCG
""")
    #                                ^ Duplicated base

    #                                           A,C,G,T,N,...,coverage
    expected_section = """\
SARS-CoV-2-seed,SARS-CoV-2-nsp12,15,58,59,13500,9,0,0,0,0,0,0,0,0,9,
SARS-CoV-2-seed,SARS-CoV-2-nsp12,15,59,60,13501,0,9,0,0,0,0,0,0,0,9,
SARS-CoV-2-seed,SARS-CoV-2-nsp12,15,60,61,13502,0,9,0,0,0,0,0,0,0,9,
SARS-CoV-2-seed,SARS-CoV-2-nsp12,15,61,62,13503,0,0,9,0,0,0,0,0,0,9,"""

    report_file = StringIO()
    default_sequence_report.write_nuc_header(report_file)
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_nuc_counts()

    report = report_file.getvalue()
    report_lines = report.splitlines()
    expected_size = 131
    if len(report_lines) != expected_size:
        assert (len(report_lines), report) == (expected_size, '')

    key_lines = report_lines[127:131]
    key_report = '\n'.join(key_lines)
    assert key_report == expected_section


# noinspection DuplicatedCode
def test_skipped_nucleotide_amino(default_sequence_report):
    """ Ensure that the aminos are aligned and output correctly for the skipped nucleotide in HXB2"""

    # refname,qcut,rank,count,offset,seq
    aligned_reads = prepare_reads("""\
HIV1-B-FR-K03455-seed,15,0,9,1,CAACAACTGCTGTTTATCCATTTTCAGAATTGGGTGTCGACATAGCAGAA
""")
    # skipped pos is here:                          ^
    # expected amino sequence: QQLLFIHFRIGCRHSR
    #                                    A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,...,coverage
    expected_section = """\
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,14,69,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,17,70,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,20,71,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,24,72,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,27,73,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,9
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,30,74,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,33,75,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9"""

    report_file = StringIO()
    default_sequence_report.write_amino_header(report_file)
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_amino_counts()

    report = report_file.getvalue()
    report_lines = report.splitlines()
    expected_size = 17
    if len(report_lines) != expected_size:
        assert (len(report_lines), report) == (expected_size, '')
    key_lines = report_lines[5:12]
    key_report = '\n'.join(key_lines)
    assert key_report == expected_section


# noinspection DuplicatedCode
def test_no_skipped_nucleotide_amino(default_sequence_report):
    """ Ensure that the aminos are aligned and output correctly if there is no skipped nucleotide in HXB2"""

    # refname,qcut,rank,count,offset,seq
    aligned_reads = prepare_reads("""\
HIV1-B-FR-K03455-seed,15,0,9,1,CAACAACTGCTGTTTATCCATTTCAGAATTGGGTGTCGACATAGCAGAA
""")
    # expected amino sequence: QQLLFIHFRIGCRHSR
    #                                    A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,...,coverage
    expected_section = """\
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,14,69,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,17,70,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,20,71,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,23,72,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,26,73,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,9
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,29,74,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,32,75,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9"""

    report_file = StringIO()
    default_sequence_report.write_amino_header(report_file)
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_amino_counts()

    report = report_file.getvalue()
    report_lines = report.splitlines()
    expected_size = 17
    if len(report_lines) != expected_size:
        assert (len(report_lines), report) == (expected_size, '')
    key_lines = report_lines[5:12]
    key_report = '\n'.join(key_lines)
    assert key_report == expected_section


# noinspection DuplicatedCode
def test_skipped_nucleotide_nuc(default_sequence_report):
    """ Ensure that the nucleotides are output correctly for the skipped nucleotide in HXB2"""

    # refname,qcut,rank,count,offset,seq
    aligned_reads = prepare_reads("""\
HIV1-B-FR-K03455-seed,15,0,9,1,CAACAACTGCTGTTTATCCATTTTCAGAATTGGGTGTCGACATAGCAGAA
""")
    # skipped pos is here:                          ^
    # skipped pos is 5772 in the genome, and 21 within this read

    expected_section = """\
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,21,212,5770,9,0,0,0,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,22,213,5771,0,0,0,9,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,23,214,5772,0,0,0,9,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,24,215,5773,0,0,0,9,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,25,216,5774,0,0,0,9,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,26,217,5775,0,9,0,0,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,27,218,5776,9,0,0,0,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,28,219,5777,0,0,9,0,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,29,220,5778,9,0,0,0,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,30,221,5779,9,0,0,0,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,31,222,5780,0,0,0,9,0,0,0,0,0,9,"""

    report_file = StringIO()
    default_sequence_report.write_nuc_header(report_file)
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_nuc_counts()
    report = report_file.getvalue()
    report_lines = report.splitlines()
    expected_size = 51
    if len(report_lines) != expected_size:
        assert (len(report_lines), report) == (expected_size, '')
    key_lines = report_lines[20:31]
    key_report = '\n'.join(key_lines)
    assert key_report == expected_section


# noinspection DuplicatedCode
def test_no_skipped_nucleotide_nuc(default_sequence_report):
    """ Ensure that the nucleotides are output correctly if the skipped nucleotide in HXB2 is not present"""

    # refname,qcut,rank,count,offset,seq
    aligned_reads = prepare_reads("""\
HIV1-B-FR-K03455-seed,15,0,9,1,CAACAACTGCTGTTTATCCATTTCAGAATTGGGTGTCGACATAGCAGAA
""")
    # skipped pos is 5772 in the genome

    expected_section = """\
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,21,212,5770,9,0,0,0,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,22,213,5771,0,0,0,9,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,,214,5772,0,0,0,0,0,9,0,0,0,9,
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,23,215,5773,0,0,0,9,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,24,216,5774,0,0,0,9,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,25,217,5775,0,9,0,0,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,26,218,5776,9,0,0,0,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,27,219,5777,0,0,9,0,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,28,220,5778,9,0,0,0,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,29,221,5779,9,0,0,0,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,30,222,5780,0,0,0,9,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,31,223,5781,0,0,0,9,0,0,0,0,0,9,"""

    report_file = StringIO()
    default_sequence_report.write_nuc_header(report_file)
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_nuc_counts()
    report = report_file.getvalue()
    report_lines = report.splitlines()
    expected_size = 51
    if len(report_lines) != expected_size:
        assert (len(report_lines), report) == (expected_size, '')
    key_lines = report_lines[20:32]
    key_report = '\n'.join(key_lines)
    assert key_report == expected_section


def test_whole_genome_consensus(default_sequence_report):
    """ Check that the whole genome consensus is correctly added from different regions.

    The given read spans more than 1 region.
    """
    aligned_reads = prepare_reads("""\
HIV1-B-FR-K03455-seed,15,0,9,1,AGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCAGTGGCGCCCGAACA\
GGGACCTGAAAGCGAAAGGGAAACCAGAGGAGCTCTCTCGACGCAGGACTCG
""")

    expected_text = """\
seed,q-cutoff,consensus-percent-cutoff,offset,sequence
HIV1-B-FR-K03455-seed,15,MAX,600,\
AGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCAGTGGCGCCCGAACA\
GGGACCTGAAAGCGAAAGGGAAACCAGAGGAGCTCTCTCGACGCAGGACTCG
HIV1-B-FR-K03455-seed,15,0.100,600,\
AGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCAGTGGCGCCCGAACA\
GGGACCTGAAAGCGAAAGGGAAACCAGAGGAGCTCTCTCGACGCAGGACTCG
"""

    report_file = StringIO()
    default_sequence_report.write_consensus_stitched_header(report_file)
    default_sequence_report.read(aligned_reads)
    default_sequence_report.combine_reports()
    default_sequence_report.write_whole_genome_consensus_from_nuc()
    assert report_file.getvalue() == expected_text


# noinspection DuplicatedCode
def test_whole_genome_consensus_amino_insertions(default_sequence_report):
    """ Check that insertions are correctly inserted into the whole genome consensus from a translated region"""
    aligned_reads = prepare_reads("""\
HIV1-B-FR-K03455-seed,15,0,9,1,\
CAGAAAAATTGTGGGTCACAGTCTATTATGGGAAAAAAAAAGTACCTGTGTGGAAGGAAGCAACCACCACTCTATTTTGTGCATCAGATGCTAA
""")
    # this is the ref genome from pos 6316 to 6400 (1 based)
    # plus an insertion of 'AAAAAAAAA' here
    #                           ^^^^^^^^^
    # the insertion is between nucleotide 6347 and 6348 (1 based)

    expected_text = """\
seed,q-cutoff,consensus-percent-cutoff,offset,sequence
HIV1-B-FR-K03455-seed,15,MAX,6315,\
CAGAAAAATTGTGGGTCACAGTCTATTATGGGAAAAAAAAAGTACCTGTGTGGAAGGAAGCAACCACCACTCTATTTTGTGCATCAGATGCTAA
HIV1-B-FR-K03455-seed,15,0.100,6315,\
CAGAAAAATTGTGGGTCACAGTCTATTATGGGAAAAAAAAAGTACCTGTGTGGAAGGAAGCAACCACCACTCTATTTTGTGCATCAGATGCTAA
"""
    expected_regions = """\
seed,region,q-cutoff,consensus-percent-cutoff,offset,sequence
HIV1-B-FR-K03455-seed,GP120,15,MAX,91,\
CAGAAAAATTGTGGGTCACAGTCTATTATGGGAAAAAAAAAGTACCTGTGTGGAAGGAAGCAACCACCACTCTATTTTGTGCATCAGATGCTAA
HIV1-B-FR-K03455-seed,GP120,15,0.100,91,\
CAGAAAAATTGTGGGTCACAGTCTATTATGGGAAAAAAAAAGTACCTGTGTGGAAGGAAGCAACCACCACTCTATTTTGTGCATCAGATGCTAA
"""

    report_file = StringIO()
    region_file = StringIO()
    default_sequence_report.write_consensus_stitched_header(report_file)
    default_sequence_report.write_consensus_regions_header(region_file)
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_insertions()
    default_sequence_report.combine_reports()
    default_sequence_report.write_whole_genome_consensus_from_nuc()
    default_sequence_report.write_consensus_regions()
    assert report_file.getvalue() == expected_text
    assert region_file.getvalue() == expected_regions


# noinspection DuplicatedCode
def test_whole_genome_consensus_nuc_insertions(default_sequence_report):
    """ Check that insertions are correctly inserted into the whole genome consensus from a non-translated region"""
    aligned_reads = prepare_reads("""\
HIV1-B-FR-K03455-seed,15,0,9,1,\
GGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGAAAAAAAAAGATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGG
""")
    # this is the ref genome from pos 1 to 99 (0 based) plus an insertion of 'AAAAAAAAA'
    #                                           ^^^^^^^^^
    # the insertion is between nucleotide 48 and 49

    expected_text = """\
seed,q-cutoff,consensus-percent-cutoff,offset,sequence
HIV1-B-FR-K03455-seed,15,MAX,1,\
GGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGAAAAAAAAAGATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGG
HIV1-B-FR-K03455-seed,15,0.100,1,\
GGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGAAAAAAAAAGATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGG
"""

    report_file = StringIO()
    default_sequence_report.write_consensus_stitched_header(report_file)
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_insertions()
    default_sequence_report.combine_reports()
    default_sequence_report.write_whole_genome_consensus_from_nuc()
    assert report_file.getvalue() == expected_text


# noinspection DuplicatedCode
def test_whole_genome_consensus_different_insertions(default_sequence_report):
    """ Check that mixed insertions are correctly inserted into the whole genome consensus"""
    aligned_reads = prepare_reads("""\
HIV1-B-FR-K03455-seed,15,0,12,1,\
ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGAAAAAAAAAGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
HIV1-B-FR-K03455-seed,15,0,9,1,\
ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGAAAGGGAAAGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
""")
    # this is the ref genome from pos 789 to 889 (0 based) plus two different insertions
    #                                           ^^^^^^^^^
    # the insertions are between nucleotide 48 and 49

    expected_text = """\
seed,q-cutoff,consensus-percent-cutoff,offset,sequence
HIV1-B-FR-K03455-seed,15,MAX,789,\
ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGAAAAAAAAAGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
HIV1-B-FR-K03455-seed,15,0.100,789,\
ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGAAARRRAAAGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
"""
    expected_regions = """\
seed,region,q-cutoff,consensus-percent-cutoff,offset,sequence
HIV1-B-FR-K03455-seed,HIV1B-gag,15,MAX,0,\
ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGAAAAAAAAAGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
HIV1-B-FR-K03455-seed,HIV1B-gag,15,0.100,0,\
ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGAAARRRAAAGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
"""

    report_file = StringIO()
    region_file = StringIO()
    default_sequence_report.write_consensus_stitched_header(report_file)
    default_sequence_report.write_consensus_regions_header(region_file)
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_insertions()
    default_sequence_report.combine_reports()
    default_sequence_report.write_whole_genome_consensus_from_nuc()
    default_sequence_report.write_consensus_regions()
    assert report_file.getvalue() == expected_text
    assert region_file.getvalue() == expected_regions


# noinspection DuplicatedCode
def test_whole_genome_consensus_half_insertions(default_sequence_report):
    """ Check that mixed insertions are correctly inserted into the whole genome consensus"""
    aligned_reads = prepare_reads("""\
HIV1-B-FR-K03455-seed,15,0,15,1,\
ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGAAAAAAAAAGAAAAAATTCGGTTAAGGCCA---GGGGGAAAGAAAAAATATAAATTAAAACATAT
HIV1-B-FR-K03455-seed,15,0,10,1,\
ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGG---------GAAAAAATTCGGTTAAGGCCACCCGGGGGAAAGAAAAAATATAAATTAAAACATAT
""")
    # this is the ref genome from pos 789 to 889 (0 based) plus different insertions in each of the reads
    # here and here:                            ^^^^^^^^^                     ^^^
    # In the MAX, only the insertion from the majority of reads will show up.
    # In the 0.1 mixture sequence, both insertions will show up as lower case.

    expected_text = """\
seed,q-cutoff,consensus-percent-cutoff,offset,sequence
HIV1-B-FR-K03455-seed,15,MAX,789,\
ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGAAAAAAAAAGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
HIV1-B-FR-K03455-seed,15,0.100,789,\
ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGaaaaaaaaaGAAAAAATTCGGTTAAGGCCAcccGGGGGAAAGAAAAAATATAAATTAAAACATAT
"""
    expected_regions = """\
seed,region,q-cutoff,consensus-percent-cutoff,offset,sequence
HIV1-B-FR-K03455-seed,HIV1B-gag,15,MAX,0,\
ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGAAAAAAAAAGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
HIV1-B-FR-K03455-seed,HIV1B-gag,15,0.100,0,\
ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGaaaaaaaaaGAAAAAATTCGGTTAAGGCCAcccGGGGGAAAGAAAAAATATAAATTAAAACATAT
"""

    report_file = StringIO()
    region_file = StringIO()
    default_sequence_report.write_consensus_stitched_header(report_file)
    default_sequence_report.write_consensus_regions_header(region_file)
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_insertions()
    default_sequence_report.combine_reports()
    default_sequence_report.write_whole_genome_consensus_from_nuc()
    default_sequence_report.write_consensus_regions()
    assert report_file.getvalue() == expected_text
    assert region_file.getvalue() == expected_regions


# noinspection DuplicatedCode
def test_whole_genome_consensus_minority_insertions(default_sequence_report):
    """ Check that insertions relative to the consensus are correctly inserted into the whole genome consensus"""
    aligned_reads = prepare_reads("""\
HIV1-B-FR-K03455-seed,15,0,10,1,\
ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
""")
    # this is the ref genome from pos 789 to 889 (0 based)
    conseq_ins_csv = StringIO("""\
qname,fwd_rev,refname,pos,insert,qual
Example_read_1,F,HIV1-B-FR-K03455-seed,20,AAC,AAA
Example_read_2,F,HIV1-B-FR-K03455-seed,20,AAC,AAA
""")

    expected_text = """\
seed,q-cutoff,consensus-percent-cutoff,offset,sequence
HIV1-B-FR-K03455-seed,15,MAX,789,\
ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
HIV1-B-FR-K03455-seed,15,0.100,789,\
ATGGGTGCGAGAGCGTCAGaacTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
"""
    expected_regions = """\
seed,region,q-cutoff,consensus-percent-cutoff,offset,sequence
HIV1-B-FR-K03455-seed,HIV1B-gag,15,MAX,0,\
ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
HIV1-B-FR-K03455-seed,HIV1B-gag,15,0.100,0,\
ATGGGTGCGAGAGCGTCAGaacTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
"""

    report_file = StringIO()
    region_file = StringIO()
    default_sequence_report.read_insertions(conseq_ins_csv)
    default_sequence_report.write_consensus_stitched_header(report_file)
    default_sequence_report.write_consensus_regions_header(region_file)
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_amino_header(StringIO())
    default_sequence_report.write_nuc_header(StringIO())
    default_sequence_report.write_nuc_counts()  # calculates ins counts
    default_sequence_report.write_amino_counts()
    default_sequence_report.write_insertions()
    default_sequence_report.combine_reports()
    default_sequence_report.write_whole_genome_consensus_from_nuc()
    default_sequence_report.write_consensus_regions()
    assert report_file.getvalue() == expected_text
    assert region_file.getvalue() == expected_regions


# noinspection DuplicatedCode
def test_nuc_minority_insertions(default_sequence_report):
    """ Check that insertions relative to the consensus are correctly inserted into the nuc.csv file for translated
    and untranslated region"""
    aligned_reads = prepare_reads("""\
HIV1-B-FR-K03455-seed,15,0,10,1,\
ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
""")
    # this is the ref genome from pos 789 to 889 (0 based)
    conseq_ins_csv = StringIO("""\
qname,fwd_rev,refname,pos,insert,qual
Example_read_1,F,HIV1-B-FR-K03455-seed,10,AAC,AAA
Example_read_2,F,HIV1-B-FR-K03455-seed,10,AAC,AAA
""")

    expected_text_untranslated = """\
HIV1-B-FR-K03455-seed,HIV1B-sl4,15,8,4,796,0,0,10,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-sl4,15,9,5,797,0,10,0,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-sl4,15,10,6,798,0,0,10,0,0,0,2,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-sl4,15,11,7,799,10,0,0,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-sl4,15,12,8,800,0,0,10,0,0,0,0,0,0,10,"""
    expected_text_translated = """\
HIV1-B-FR-K03455-seed,HIV1B-gag,15,8,7,796,0,0,10,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,9,8,797,0,10,0,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,10,9,798,0,0,10,0,0,0,2,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,11,10,799,10,0,0,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,12,11,800,0,0,10,0,0,0,0,0,0,10,"""

    nuc_file = StringIO()
    default_sequence_report.read_insertions(conseq_ins_csv)
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_nuc_header(nuc_file)
    default_sequence_report.write_nuc_counts()  # calculates ins counts
    report = nuc_file.getvalue()
    report_lines = report.splitlines()
    key_lines1 = report_lines[105:110]
    key_report1 = '\n'.join(key_lines1)
    key_lines2 = report_lines[7:12]
    key_report2 = '\n'.join(key_lines2)
    assert key_report1 == expected_text_untranslated
    assert key_report2 == expected_text_translated


# noinspection DuplicatedCode
def test_nuc_small_majority_insertion(default_sequence_report):
    """ Check that a small insertion relative to the reference is correctly inserted into the nuc.csv file
    and the insertions.csv file for a translated region"""
    aligned_reads = prepare_reads("""\
HIV1-B-FR-K03455-seed,15,0,10,1,\
ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAATTTTTTTTTATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
""")
    # this is the ref genome from pos 789 to 889 (0 based), with an insertion here:
    #                                                 ^^^^^^^^^

    expected_text = """\
HIV1-B-FR-K03455-seed,HIV1B-gag,15,53,52,841,10,0,0,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,54,53,842,10,0,0,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,55,54,843,10,0,0,0,0,0,10,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,65,55,844,10,0,0,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,66,56,845,0,0,0,10,0,0,0,0,0,10,"""

    expected_insertions = """\
seed,mixture_cutoff,region,ref_region_pos,ref_genome_pos,query_pos,insertion
HIV1-B-FR-K03455-seed,MAX,HIV1B-gag,54,843,55,TTTTTTTTT
HIV1-B-FR-K03455-seed,0.100,HIV1B-gag,54,843,55,TTTTTTTTT
"""

    nuc_file = StringIO()
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_nuc_header(nuc_file)
    default_sequence_report.write_insertions(default_sequence_report.insert_writer)
    default_sequence_report.write_nuc_counts()  # calculates ins counts
    report = nuc_file.getvalue()
    report_lines = report.splitlines()
    key_lines = report_lines[52:57]
    key_report = '\n'.join(key_lines)
    assert key_report == expected_text
    assert default_sequence_report.insert_writer.insert_file.getvalue() == expected_insertions


# noinspection DuplicatedCode
def test_nuc_large_majority_insertion(default_sequence_report):
    """ Check that a large (>30) insertion relative to the reference is correctly inserted into the nuc.csv file
    for translated regions"""
    aligned_reads = prepare_reads("""\
HIV1-B-FR-K03455-seed,15,0,10,1,\
ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATTCGGTTAAGG\
CCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
""")
    # this is the ref genome from pos 789 to 889 (0 based), with an insertion here:
    #                                                 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    expected_text = """\
HIV1-B-FR-K03455-seed,HIV1B-gag,15,53,52,841,10,0,0,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,54,53,842,10,0,0,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,55,54,843,10,0,0,0,0,0,10,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,89,55,844,10,0,0,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,90,56,845,0,0,0,10,0,0,0,0,0,10,"""

    nuc_file = StringIO()
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_nuc_header(nuc_file)
    default_sequence_report.write_insertions(default_sequence_report.insert_writer)
    default_sequence_report.write_nuc_counts()  # calculates ins counts
    report = nuc_file.getvalue()
    report_lines = report.splitlines()
    key_lines = report_lines[52:57]
    key_report = '\n'.join(key_lines)
    assert key_report == expected_text


# noinspection DuplicatedCode
def test_nuc_large_majority_insertion_offset(default_sequence_report):
    """ Check that a large (>30) insertion relative to the reference is correctly inserted into the nuc.csv file
    for translated regions, with a different consensus offset"""
    aligned_reads = prepare_reads("""\
HIV1-B-FR-K03455-seed,15,0,10,30,\
ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATTCGGTTAAGG\
CCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
""")
    # this is the ref genome from pos 789 to 889 (0 based), with an insertion here:
    #                                                 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    expected_text = """\
HIV1-B-FR-K03455-seed,HIV1B-gag,15,82,52,841,10,0,0,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,83,53,842,10,0,0,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,84,54,843,10,0,0,0,0,0,10,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,118,55,844,10,0,0,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,119,56,845,0,0,0,10,0,0,0,0,0,10,"""

    nuc_file = StringIO()
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_nuc_header(nuc_file)
    default_sequence_report.write_insertions(default_sequence_report.insert_writer)
    default_sequence_report.write_nuc_counts()  # calculates ins counts
    report = nuc_file.getvalue()
    report_lines = report.splitlines()
    key_lines = report_lines[52:57]
    key_report = '\n'.join(key_lines)
    assert key_report == expected_text


# noinspection DuplicatedCode
def test_nuc_large_majority_insertion_frameshift(default_sequence_report):
    """ Check that a large (>30) insertion relative to the reference is correctly inserted into the nuc.csv file
    for translated regions, even if it is out of frame"""
    aligned_reads = prepare_reads("""\
HIV1-B-FR-K03455-seed,15,0,10,1,\
ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAATTCGGTTAAGG\
CCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
""")
    # this is the ref genome from pos 789 to 889 (0 based), with an insertion here:
    #                                                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    expected_text = """\
HIV1-B-FR-K03455-seed,HIV1B-gag,15,51,50,839,10,0,0,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,52,51,840,10,0,0,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,53,52,841,10,0,0,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,54,53,842,10,0,0,0,0,0,10,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,88,54,843,10,0,0,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,89,55,844,10,0,0,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,90,56,845,0,0,0,10,0,0,0,0,0,10,"""

    nuc_file = StringIO()
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_nuc_header(nuc_file)
    default_sequence_report.write_insertions(default_sequence_report.insert_writer)
    default_sequence_report.write_nuc_counts()  # calculates ins counts
    report = nuc_file.getvalue()
    report_lines = report.splitlines()
    key_lines = report_lines[50:57]
    key_report = '\n'.join(key_lines)
    assert key_report == expected_text


# noinspection DuplicatedCode
def test_nuc_large_insertion_not_multiple_of_three(default_sequence_report):
    """ Check that a large (>30) insertion relative to the reference is correctly inserted into the nuc.csv file
    for translated regions, even if it not a multiple of three long"""
    aligned_reads = prepare_reads("""\
HIV1-B-FR-K03455-seed,15,0,10,1,\
ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAATTCGGTTAAGG\
CCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
""")
    # this is the ref genome from pos 789 to 889 (0 based), with an insertion here:
    #                                                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    expected_text = """\
HIV1-B-FR-K03455-seed,HIV1B-gag,15,51,50,839,10,0,0,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,52,51,840,10,0,0,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,53,52,841,10,0,0,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,54,53,842,10,0,0,0,0,0,10,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,89,54,843,10,0,0,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,90,55,844,10,0,0,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,91,56,845,0,0,0,10,0,0,0,0,0,10,"""

    nuc_file = StringIO()
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_nuc_header(nuc_file)
    default_sequence_report.write_insertions(default_sequence_report.insert_writer)
    default_sequence_report.write_nuc_counts()  # calculates ins counts
    report = nuc_file.getvalue()
    report_lines = report.splitlines()
    key_lines = report_lines[50:57]
    key_report = '\n'.join(key_lines)
    assert key_report == expected_text


def test_merge_extra_counts_insertion(projects, default_sequence_report):
    """ Check that an insertion after a deletion is correctly handled by merge_extra_counts """
    seed_name = 'HIV1-B-FR-K03455-seed'
    seed_seq = projects.getReference(seed_name)
    read_seq = seed_seq[789:888] + seed_seq[900:999] + 'TTTTTTTTTTTT' + seed_seq[999:1200]
    aligned_reads = prepare_reads(f"""\
HIV1-B-FR-K03455-seed,15,0,10,1,{read_seq}
""")
    # seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,
    # partial,del,ins,clip,v3_overlap,coverage
    expected_deletion = """\
HIV1-B-FR-K03455-seed,HIV1B-gag,15,98,33,0,0,0,0,0,0,10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10
HIV1-B-FR-K03455-seed,HIV1B-gag,15,,34,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10,0,0,0,10
HIV1-B-FR-K03455-seed,HIV1B-gag,15,,35,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10,0,0,0,10
HIV1-B-FR-K03455-seed,HIV1B-gag,15,,36,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10,0,0,0,10
HIV1-B-FR-K03455-seed,HIV1B-gag,15,,37,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10,0,0,0,10
HIV1-B-FR-K03455-seed,HIV1B-gag,15,101,38,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10,0,0,0,0,0,0,0,0,0,0,0,10"""
    expected_insertion_amino = """\
HIV1-B-FR-K03455-seed,HIV1B-gag,15,194,69,0,0,0,0,0,0,0,0,0,0,0,0,0,10,0,0,0,0,0,0,0,0,0,0,0,0,0,10
HIV1-B-FR-K03455-seed,HIV1B-gag,15,197,70,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10,0,0,0,0,0,0,0,10,0,0,10
HIV1-B-FR-K03455-seed,HIV1B-gag,15,212,71,0,0,0,0,0,10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10"""
    # seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,genome.pos,
    # A,C,G,T,N,del,ins,clip,v3_overlap,coverage,exact_coverage
    expected_insertion = """\
HIV1-B-FR-K03455-seed,HIV1B-gag,15,197,208,997,10,0,0,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,198,209,998,0,10,0,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,199,210,999,10,0,0,0,0,0,10,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,212,211,1000,0,0,10,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,213,212,1001,0,0,10,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,214,213,1002,10,0,0,0,0,0,0,0,0,10,"""

    nuc_csv = StringIO()
    amino_csv = StringIO()
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_nuc_header(nuc_csv)
    default_sequence_report.write_amino_header(amino_csv)
    default_sequence_report.write_insertions(default_sequence_report.insert_writer)
    default_sequence_report.write_nuc_counts()  # calls merge_extra_counts
    default_sequence_report.write_amino_counts()
    report_lines = nuc_csv.getvalue().splitlines()
    amino_lines = amino_csv.getvalue().splitlines()
    deletion_lines = amino_lines[33:39]
    amino_insertion_lines = amino_lines[69:72]
    insertion_lines = report_lines[208:214]
    deletion_report = '\n'.join(deletion_lines)
    amino_insertion_report = '\n'.join(amino_insertion_lines)
    insertion_report = '\n'.join(insertion_lines)
    assert deletion_report == expected_deletion
    assert insertion_report == expected_insertion
    assert amino_insertion_report == expected_insertion_amino


def test_merge_extra_counts_insertion_vpr(projects, default_sequence_report):
    """ Check that an insertion in vpr is correctly handled by merge_extra_counts """
    seed_name = 'HIV1-B-FR-K03455-seed'
    seed_seq = projects.getReference(seed_name)
    read_seq = seed_seq[5558:5800] + 'AAAAAAAAA' + seed_seq[5800:5850]
    # vpr, insertion behind skip pos (5772)
    aligned_reads = prepare_reads(f"""\
HIV1-B-FR-K03455-seed,15,0,10,0,{read_seq}
""")
# seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,genome.pos,\
# A,C,G,T,N,del,ins,clip,v3_overlap,coverage,exact_coverage
    expected_insertion = """\
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,239,239,5797,10,0,0,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,240,240,5798,0,0,10,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,241,241,5799,10,0,0,0,0,0,10,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,251,242,5800,10,0,0,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,252,243,5801,0,0,0,10,0,0,0,0,0,10,"""
    expected_amino_insertion = """\
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,236,79,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10,0,0,0,0,0,0,0,0,0,0,0,10
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,239,80,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10,0,0,0,0,0,0,0,0,0,10,0,0,10
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,251,81,0,0,0,0,0,0,0,10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10"""

    nuc_csv = StringIO()
    amino_csv = StringIO()
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_nuc_header(nuc_csv)
    default_sequence_report.write_amino_header(amino_csv)
    default_sequence_report.write_insertions(default_sequence_report.insert_writer)
    default_sequence_report.write_nuc_counts()  # calls merge_extra_counts
    default_sequence_report.write_amino_counts()
    report_lines = nuc_csv.getvalue().splitlines()
    insertion_report = '\n'.join(report_lines[320:325])
    amino_report_lines = amino_csv.getvalue().splitlines()
    amino_report = '\n'.join(amino_report_lines[105:108])
    assert insertion_report == expected_insertion
    assert amino_report == expected_amino_insertion


def test_merge_extra_counts_insertion_vpr_noskip(projects, default_sequence_report):
    """ Check that an insertion in vpr without the skip pos is correctly handled by merge_extra_counts """
    seed_name = 'HIV1-B-FR-K03455-seed'
    seed_seq = projects.getReference(seed_name)
    read_seq = seed_seq[5558:5771] + seed_seq[5772:5800] + 'AAAAAAAAA' + seed_seq[5800:5850]
    # vpr, insertion behind skip pos (5772)
    aligned_reads = prepare_reads(f"""\
HIV1-B-FR-K03455-seed,15,0,10,0,{read_seq}
""")
# seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,genome.pos,\
# A,C,G,T,N,del,ins,clip,v3_overlap,coverage,exact_coverage
    expected_insertion = """\
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,238,239,5797,10,0,0,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,239,240,5798,0,0,10,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,240,241,5799,10,0,0,0,0,0,10,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,250,242,5800,10,0,0,0,0,0,0,0,0,10,
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,251,243,5801,0,0,0,10,0,0,0,0,0,10,"""
    expected_amino_insertion = """\
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,235,79,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10,0,0,0,0,0,0,0,0,0,0,0,10
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,238,80,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10,0,0,0,0,0,0,0,0,0,10,0,0,10
HIV1-B-FR-K03455-seed,HIV1B-vpr,15,250,81,0,0,0,0,0,0,0,10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10"""

    nuc_csv = StringIO()
    amino_csv = StringIO()
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_nuc_header(nuc_csv)
    default_sequence_report.write_amino_header(amino_csv)
    default_sequence_report.write_insertions(default_sequence_report.insert_writer)
    default_sequence_report.write_nuc_counts()  # calls merge_extra_counts
    default_sequence_report.write_amino_counts()
    report_lines = nuc_csv.getvalue().splitlines()
    insertion_report = '\n'.join(report_lines[320:325])
    amino_report_lines = amino_csv.getvalue().splitlines()
    amino_report = '\n'.join(amino_report_lines[105:108])
    assert insertion_report == expected_insertion
    assert amino_report == expected_amino_insertion


def test_merge_extra_counts_insertion_nsp12(projects, default_sequence_report):
    """ Check that an insertion after the duplicated position is correctly handled by merge_extra_counts """
    seed_name = 'SARS-CoV-2-seed'
    seed_seq = projects.getReference(seed_name)
    read_seq = seed_seq[13441:13600] + 'AAAAAAAAA' + seed_seq[13600:13700]
    # nsp12, insertion behind repeated pos (13468)
    aligned_reads = prepare_reads(f"""\
SARS-CoV-2-seed,15,0,10,0,{read_seq}
""")
# seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,genome.pos,\
# A,C,G,T,N,del,ins,clip,v3_overlap,coverage,exact_coverage
    expected_insertion = """\
SARS-CoV-2-seed,SARS-CoV-2-nsp12,15,157,157,13598,0,0,10,0,0,0,0,0,0,10,
SARS-CoV-2-seed,SARS-CoV-2-nsp12,15,158,158,13599,0,0,0,10,0,0,10,0,0,10,
SARS-CoV-2-seed,SARS-CoV-2-nsp12,15,168,159,13600,10,0,0,0,0,0,0,0,0,10,"""
    expected_amino_insertion = """\
SARS-CoV-2-seed,SARS-CoV-2-nsp12,15,153,52,0,0,0,0,0,0,0,0,0,0,0,10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10
SARS-CoV-2-seed,SARS-CoV-2-nsp12,15,156,53,0,10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10,0,0,10
SARS-CoV-2-seed,SARS-CoV-2-nsp12,15,168,54,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10,0,0,0,0,0,0,0,0,0,0,0,10"""

    nuc_csv = StringIO()
    amino_csv = StringIO()
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_nuc_header(nuc_csv)
    default_sequence_report.write_amino_header(amino_csv)
    default_sequence_report.write_insertions(default_sequence_report.insert_writer)
    default_sequence_report.write_nuc_counts()  # calls merge_extra_counts
    default_sequence_report.write_amino_counts()
    report_lines = nuc_csv.getvalue().splitlines()
    insertion_report = '\n'.join(report_lines[455:458])
    amino_report_lines = amino_csv.getvalue().splitlines()
    amino_report = '\n'.join(amino_report_lines[151:154])
    assert insertion_report == expected_insertion
    assert amino_report == expected_amino_insertion


# noinspection DuplicatedCode
def test_majority_and_minority_insertion(default_sequence_report):
    """ Check that an insertion relative to the reference is correctly inserted into conseq_regions
    for translated regions, even if there are also insertions relative to the consensus at that position"""
    aligned_reads = prepare_reads("""\
1-HIV1-B-FR-K03455-seed,15,0,10,1,\
ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAATTCGGTTAAGG\
CCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
""")
    # this is the ref genome from pos 789 to 889 (0 based), with an insertion here:
    #                                                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    conseq_ins_csv = StringIO("""\
qname,fwd_rev,refname,pos,insert,qual
Example_read_1,F,1-HIV1-B-FR-K03455-seed,54,AAC,AAA
Example_read_2,F,1-HIV1-B-FR-K03455-seed,54,AAC,AAA
""")

    expected_regions = """\
seed,region,q-cutoff,consensus-percent-cutoff,offset,sequence
HIV1-B-FR-K03455-seed,HIV1B-gag,15,MAX,0,\
ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAATTCGGTTAAGG\
CCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
HIV1-B-FR-K03455-seed,HIV1B-gag,15,0.100,0,\
ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAaacTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAATTCGGTTAAGG\
CCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
"""

    regions_file = StringIO()
    default_sequence_report.read_insertions(conseq_ins_csv)
    default_sequence_report.write_consensus_stitched_header(StringIO())
    default_sequence_report.write_consensus_regions_header(regions_file)
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_nuc_header(StringIO())
    default_sequence_report.write_insertions(default_sequence_report.insert_writer)
    default_sequence_report.write_nuc_counts()  # calculates ins counts
    default_sequence_report.combine_reports()
    default_sequence_report.write_whole_genome_consensus_from_nuc()
    default_sequence_report.write_consensus_regions()
    assert regions_file.getvalue() == expected_regions


# noinspection DuplicatedCode
def test_insertion_end_of_region(default_sequence_report):
    """ Check that an insertion relative to the consensus is not inserted into conseq_regions
    if it is located at the very end of the region"""
    aligned_reads = prepare_reads("""\
1-HIV1-B-FR-K03455-seed,15,0,10,0,\
CTCCCCCTCAGAAGCAGGAGCCGATAGACAAGGAACTGTATCCTTTAACTTCCCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGG
""")
    # this is the ref genome from pos 2200 to 2300 (0 based)

    conseq_ins_csv = StringIO("""\
qname,fwd_rev,refname,pos,insert,qual
Example_read_1,F,1-HIV1-B-FR-K03455-seed,92,AAC,AAA
Example_read_2,F,1-HIV1-B-FR-K03455-seed,92,AAC,AAA
Example_read_3,F,1-HIV1-B-FR-K03455-seed,92,AAC,AAA
""")

    expected_regions = """\
seed,region,q-cutoff,consensus-percent-cutoff,offset,sequence
HIV1-B-FR-K03455-seed,HIV1B-gag,15,MAX,1411,\
CTCCCCCTCAGAAGCAGGAGCCGATAGACAAGGAACTGTATCCTTTAACTTCCCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAA
HIV1-B-FR-K03455-seed,HIV1B-gag,15,0.100,1411,\
CTCCCCCTCAGAAGCAGGAGCCGATAGACAAGGAACTGTATCCTTTAACTTCCCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAA"""

    regions_file = StringIO()
    default_sequence_report.read_insertions(conseq_ins_csv)
    default_sequence_report.write_consensus_stitched_header(StringIO())
    default_sequence_report.write_consensus_regions_header(regions_file)
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_nuc_header(StringIO())
    default_sequence_report.write_insertions(default_sequence_report.insert_writer)
    default_sequence_report.write_nuc_counts()  # calculates ins counts
    default_sequence_report.combine_reports()
    default_sequence_report.write_whole_genome_consensus_from_nuc()
    default_sequence_report.write_consensus_regions()
    report = regions_file.getvalue()
    report_lines = report.splitlines()
    key_lines = report_lines[0:3]
    key_report = '\n'.join(key_lines)
    assert key_report == expected_regions


# noinspection DuplicatedCode
def test_insertion_10percent(default_sequence_report):
    """ Check that an insertion relative to the consensus is inserted into conseq_regions
    if its coverage is exactly 10%"""
    aligned_reads = prepare_reads("""\
1-HIV1-B-FR-K03455-seed,15,0,9,0,\
CTCCCCCTCAGAAGCAGGAGCCGATAGACAAGGAACTGTATCCTTTAACTTCCCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGG
""")
    # this is the ref genome from pos 2200 to 2300 (0 based)

    conseq_ins_csv = StringIO("""\
qname,fwd_rev,refname,pos,insert,qual
Example_read_1,F,1-HIV1-B-FR-K03455-seed,91,AAC,AAA
""")

    expected_regions = """\
seed,region,q-cutoff,consensus-percent-cutoff,offset,sequence
HIV1-B-FR-K03455-seed,HIV1B-gag,15,MAX,1411,\
CTCCCCCTCAGAAGCAGGAGCCGATAGACAAGGAACTGTATCCTTTAACTTCCCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAA
HIV1-B-FR-K03455-seed,HIV1B-gag,15,0.100,1411,\
CTCCCCCTCAGAAGCAGGAGCCGATAGACAAGGAACTGTATCCTTTAACTTCCCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAaacA"""

    regions_file = StringIO()
    default_sequence_report.read_insertions(conseq_ins_csv)
    default_sequence_report.write_consensus_stitched_header(StringIO())
    default_sequence_report.write_consensus_regions_header(regions_file)
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_nuc_header(StringIO())
    default_sequence_report.write_insertions(default_sequence_report.insert_writer)
    default_sequence_report.write_nuc_counts()  # calculates ins counts
    default_sequence_report.combine_reports()
    default_sequence_report.write_whole_genome_consensus_from_nuc()
    default_sequence_report.write_consensus_regions()
    report = regions_file.getvalue()
    report_lines = report.splitlines()
    key_lines = report_lines[0:3]
    key_report = '\n'.join(key_lines)
    assert key_report == expected_regions


# noinspection DuplicatedCode
def test_whole_genome_consensus_different_minority_insertions(default_sequence_report):
    """ Check that different insertions relative to the consensus are correctly inserted
    into the whole genome consensus"""
    aligned_reads1 = prepare_reads("""\
1-HIV1-B-FR-K03455-seed,15,0,10,1,\
ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
""")
# this is the ref genome from pos 789 to 889 (0 based)
    aligned_reads2 = prepare_reads("""\
2-HIV1-B-FR-K03455-seed,15,0,10,1,\
GGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
""")
# this is the ref genome from pos 792 to 889 (0 based)
    conseq_ins_csv = StringIO("""\
qname,fwd_rev,refname,pos,insert,qual
Example_read_1,F,1-HIV1-B-FR-K03455-seed,20,AAC,AAA
Example_read_2,F,1-HIV1-B-FR-K03455-seed,20,AAC,AAA
Example_read_1,F,2-HIV1-B-FR-K03455-seed,23,TTT,AAA
Example_read_2,F,2-HIV1-B-FR-K03455-seed,23,TTT,AAA
""")

    expected_text = """\
seed,q-cutoff,consensus-percent-cutoff,offset,sequence
HIV1-B-FR-K03455-seed,15,MAX,789,\
ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
HIV1-B-FR-K03455-seed,15,0.100,789,\
ATGGGTGCGAGAGCGTCAGaacTATTAAtttGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
"""
    expected_regions = """\
seed,region,q-cutoff,consensus-percent-cutoff,offset,sequence
HIV1-B-FR-K03455-seed,HIV1B-gag,15,MAX,0,\
ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
HIV1-B-FR-K03455-seed,HIV1B-gag,15,0.100,0,\
ATGGGTGCGAGAGCGTCAGaacTATTAAtttGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
"""

    report_file = StringIO()
    region_file = StringIO()
    default_sequence_report.read_insertions(conseq_ins_csv)
    default_sequence_report.write_consensus_stitched_header(report_file)
    default_sequence_report.write_consensus_regions_header(region_file)
    default_sequence_report.read(aligned_reads1)
    default_sequence_report.write_amino_header(StringIO())
    default_sequence_report.write_nuc_header(StringIO())
    default_sequence_report.write_nuc_counts()  # calculates ins counts
    default_sequence_report.write_amino_counts()
    default_sequence_report.write_insertions()
    default_sequence_report.combine_reports()
    default_sequence_report.read(aligned_reads2)
    default_sequence_report.write_amino_header(StringIO())
    default_sequence_report.write_nuc_header(StringIO())
    default_sequence_report.write_nuc_counts()  # calculates ins counts
    default_sequence_report.write_amino_counts()
    default_sequence_report.write_insertions()
    default_sequence_report.combine_reports()
    default_sequence_report.write_whole_genome_consensus_from_nuc()
    default_sequence_report.write_consensus_regions()
    assert report_file.getvalue() == expected_text
    assert region_file.getvalue() == expected_regions


# noinspection DuplicatedCode
def test_whole_genome_consensus_different_contig_insertions(default_sequence_report):
    """ Check that different insertions relative to the consensus are correctly inserted
    into the whole genome consensus"""
    aligned_reads1 = prepare_reads("""\
1-HIV1-B-FR-K03455-seed,15,0,10,1,\
ATGGGTGCGAGAGCGTCAGTACCCTTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
""")
# this is the ref genome from pos 789 to 889 (0 based) plus an insertion here:
#                    ^^^
    aligned_reads2 = prepare_reads("""\
2-HIV1-B-FR-K03455-seed,15,0,5,1,\
ATGGGTGCGAGAGCGTCAGTATTATTTAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
""")
# this is the ref genome from pos 789 to 889 (0 based) plus an insertion here:
#                       ^^^

    expected_text = """\
seed,q-cutoff,consensus-percent-cutoff,offset,sequence
HIV1-B-FR-K03455-seed,15,MAX,789,\
ATGGGTGCGAGAGCGTCAGTACCCTTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
HIV1-B-FR-K03455-seed,15,0.100,789,\
ATGGGTGCGAGAGCGTCAGTAcccTTAtttAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
"""
    expected_regions = """\
seed,region,q-cutoff,consensus-percent-cutoff,offset,sequence
HIV1-B-FR-K03455-seed,HIV1B-gag,15,MAX,0,\
ATGGGTGCGAGAGCGTCAGTACCCTTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
HIV1-B-FR-K03455-seed,HIV1B-gag,15,0.100,0,\
ATGGGTGCGAGAGCGTCAGTAcccTTAtttAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATAT
"""

    report_file = StringIO()
    region_file = StringIO()
    default_sequence_report.write_consensus_stitched_header(report_file)
    default_sequence_report.write_consensus_regions_header(region_file)
    default_sequence_report.read(aligned_reads1)
    default_sequence_report.write_insertions()
    default_sequence_report.combine_reports()
    default_sequence_report.read(aligned_reads2)
    default_sequence_report.write_insertions()
    default_sequence_report.combine_reports()
    default_sequence_report.write_whole_genome_consensus_from_nuc()
    default_sequence_report.write_consensus_regions()
    assert report_file.getvalue() == expected_text
    assert region_file.getvalue() == expected_regions


# noinspection DuplicatedCode
def test_whole_genome_consensus_insertions_overlap(default_sequence_report):
    """ Check that insertions in overlapping regions are correctly inserted into the whole genome consensus"""
    aligned_reads = prepare_reads("""\
HIV1-B-FR-K03455-seed,15,0,9,1,\
TGAGAGTGAAGGAGAAATATCAGCACTTGTGGAGATGGGGGTGGAGATGGGGCACCATGAAAAAAAAACTCCTTGGGATGTTGATGATCTGTAGTGCTA
""")
    # this is the ref genome from pos 6225 to 6315 (0 based), plus an insertion of 'AAAAAAAAA' here
    #                                                      ^^^^^^^^^
    # the insertion is between nucleotide 6283 and 6284 (0 based)

    expected_text = """\
seed,q-cutoff,consensus-percent-cutoff,offset,sequence
HIV1-B-FR-K03455-seed,15,MAX,6225,\
TGAGAGTGAAGGAGAAATATCAGCACTTGTGGAGATGGGGGTGGAGATGGGGCACCATGAAAAAAAAACTCCTTGGGATGTTGATGATCTGTAGTGCTA
HIV1-B-FR-K03455-seed,15,0.100,6225,\
TGAGAGTGAAGGAGAAATATCAGCACTTGTGGAGATGGGGGTGGAGATGGGGCACCATGAAAAAAAAACTCCTTGGGATGTTGATGATCTGTAGTGCTA
"""
    expected_regions = """\
seed,region,q-cutoff,consensus-percent-cutoff,offset,sequence
HIV1-B-FR-K03455-seed,GP120,15,MAX,1,\
TGAGAGTGAAGGAGAAATATCAGCACTTGTGGAGATGGGGGTGGAGATGGGGCACCATGAAAAAAAAACTCCTTGGGATGTTGATGATCTGTAGTGCTA
HIV1-B-FR-K03455-seed,GP120,15,0.100,1,\
TGAGAGTGAAGGAGAAATATCAGCACTTGTGGAGATGGGGGTGGAGATGGGGCACCATGAAAAAAAAACTCCTTGGGATGTTGATGATCTGTAGTGCTA
HIV1-B-FR-K03455-seed,HIV1B-vpu,15,MAX,164,\
TGAGAGTGAAGGAGAAATATCAGCACTTGTGGAGATGGGGGTGGAGATGGGGCACCATGAAAAAAAAACTCCTTGGGATGTTGATGATCTGTAG
HIV1-B-FR-K03455-seed,HIV1B-vpu,15,0.100,164,\
TGAGAGTGAAGGAGAAATATCAGCACTTGTGGAGATGGGGGTGGAGATGGGGCACCATGAAAAAAAAACTCCTTGGGATGTTGATGATCTGTAG
"""

    report_file = StringIO()
    region_file = StringIO()
    default_sequence_report.write_consensus_stitched_header(report_file)
    default_sequence_report.write_consensus_regions_header(region_file)
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_insertions()
    default_sequence_report.combine_reports()
    default_sequence_report.write_whole_genome_consensus_from_nuc()
    default_sequence_report.write_consensus_regions()
    assert report_file.getvalue() == expected_text
    assert region_file.getvalue() == expected_regions


# noinspection DuplicatedCode
def test_deleted_nucleotide(default_sequence_report):
    """ Make sure that there is no lost nucleotide after a single amino acid that aligns
    in a different frame"""

    aligned_reads = prepare_reads("""\
HIV1-B-FR-K03455-seed,15,0,10,1,ATTCATCTGTATTACTTTGACTGTTTTTCAGACTCTGCTATAAGAA\
AGGCCTTATTAGGACACATAGTTAGCCCTAGGTGTGAATATCAAGCAGGACATAACAAGGTAGGATCTCTACAATACT\
TGGCACTGCATGCATTAATAACACCAAAAAAGATAAAGCCACCTTTGCCTAGTGTTACGAAACTGACAGAGGATAGAT\
GGAACAAGCCCCAGAAGACCAAGGGCCACAGAGGGAGCCACACAATGAATGGACACTAG
""")
    #  ^ lost nucleotide is G on the second line from the end, in TGC
    expected_text = """\
seed,q-cutoff,consensus-percent-cutoff,offset,sequence
HIV1-B-FR-K03455-seed,15,MAX,5358,\
ATTCATCTGTATTACTTTGACTGTTTTTCAGACTCTGCTATAAGAA\
AGGCCTTATTAGGACACATAGTTAGCCCTAGGTGTGAATATCAAGCAGGACATAACAAGGTAGGATCTCTACAATACT\
TGGCACTGCATGCATTAATAACACCAAAAAAGATAAAGCCACCTTTGCCTAGTGTTACGAAACTGACAGAGGATAGAT\
GGAACAAGCCCCAGAAGACCAAGGGCCACAGAGGGAGCCACACAATGAATGGACACTAG
HIV1-B-FR-K03455-seed,15,0.100,5358,\
ATTCATCTGTATTACTTTGACTGTTTTTCAGACTCTGCTATAAGAA\
AGGCCTTATTAGGACACATAGTTAGCCCTAGGTGTGAATATCAAGCAGGACATAACAAGGTAGGATCTCTACAATACT\
TGGCACTGCATGCATTAATAACACCAAAAAAGATAAAGCCACCTTTGCCTAGTGTTACGAAACTGACAGAGGATAGAT\
GGAACAAGCCCCAGAAGACCAAGGGCCACAGAGGGAGCCACACAATGAATGGACACTAG
"""
    report_file = StringIO()
    default_sequence_report.write_consensus_stitched_header(report_file)
    default_sequence_report.read(aligned_reads)
    default_sequence_report.combine_reports()
    default_sequence_report.write_whole_genome_consensus_from_nuc()
    assert report_file.getvalue() == expected_text


# noinspection DuplicatedCode
def test_deleted_nucleotide_vpr(default_sequence_report):
    """ Make sure that there is no lost nucleotide after a single amino acid that aligns
    in a different frame, specifically if this happens before the skip position in vpr"""

    aligned_reads = prepare_reads("""\
HIV1-B-FR-K03455-seed,15,0,10,1,\
ATGGAACAAGCCCCAGAAGACCAAGGGCCACAGAGGGAGCCACACAATGAATGGACACTAGAGCTTTTAGAGGAGCTT\
AAGAATGAAGCTGTTAGACATTTTCCTAGGATTTGGCTCCATGGCTTAGGGCAACATATCTATGAAACTTATGGGGAT\
ACTTGGGCAGGAGTGGAAGCCATAATAAGAATTCTGCAACAACTGCTGTTTATTCATTTCAGAATTGGGTGTCGACAT\
AGCAGAATAGGCGTTACTCGACAGAGGAGAGCAAGAAATGGAGCCAGTAGATCCTAG
""")
    #  ^ lost nucleotide is C towards the end of the 3rd line, in TTTCAGAA
    expected_text = """\
seed,q-cutoff,consensus-percent-cutoff,offset,sequence
HIV1-B-FR-K03455-seed,15,MAX,5558,\
ATGGAACAAGCCCCAGAAGACCAAGGGCCACAGAGGGAGCCACACAATGAATGGACACTAGAGCTTTTAGAGGAGCTT\
AAGAATGAAGCTGTTAGACATTTTCCTAGGATTTGGCTCCATGGCTTAGGGCAACATATCTATGAAACTTATGGGGAT\
ACTTGGGCAGGAGTGGAAGCCATAATAAGAATTCTGCAACAACTGCTGTTTATTCATTTCAGAATTGGGTGTCGACAT\
AGCAGAATAGGCGTTACTCGACAGAGGAGAGCAAGAAATGGAGCCAGTAGATCCTAG
HIV1-B-FR-K03455-seed,15,0.100,5558,\
ATGGAACAAGCCCCAGAAGACCAAGGGCCACAGAGGGAGCCACACAATGAATGGACACTAGAGCTTTTAGAGGAGCTT\
AAGAATGAAGCTGTTAGACATTTTCCTAGGATTTGGCTCCATGGCTTAGGGCAACATATCTATGAAACTTATGGGGAT\
ACTTGGGCAGGAGTGGAAGCCATAATAAGAATTCTGCAACAACTGCTGTTTATTCATTTCAGAATTGGGTGTCGACAT\
AGCAGAATAGGCGTTACTCGACAGAGGAGAGCAAGAAATGGAGCCAGTAGATCCTAG
"""
    report_file = StringIO()
    default_sequence_report.write_consensus_stitched_header(report_file)
    default_sequence_report.read(aligned_reads)
    default_sequence_report.combine_reports()
    default_sequence_report.write_whole_genome_consensus_from_nuc()
    assert report_file.getvalue() == expected_text


# noinspection DuplicatedCode
def test_consensus_translated_region_differences(caplog):
    """ Test that the consensus is stitched together correctly between two translated regions """
    nucleotide1 = SeedNucleotide()
    nucleotide2 = SeedNucleotide()
    nucleotide3 = SeedNucleotide()
    nucleotide4 = SeedNucleotide()
    nucleotide5 = SeedNucleotide()
    nucleotide6 = SeedNucleotide()
    nucleotide1.count_nucleotides('A', 3)
    nucleotide2.count_nucleotides('C', 3)
    nucleotide3.count_nucleotides('G', 3)
    nucleotide4.count_nucleotides('T', 3)
    nucleotide5.count_nucleotides('T', 6)
    nucleotide6.count_nucleotides('A', 6)
    nuc_dict = {0: nucleotide1, 1: nucleotide2, 2: nucleotide3, 3: nucleotide1, 4: nucleotide1}

    region_nucleotides = [ReportNucleotide(1, seed_nucleotide=nucleotide4),
                          ReportNucleotide(2, seed_nucleotide=nucleotide4),
                          ReportNucleotide(3, seed_nucleotide=nucleotide4),
                          ReportNucleotide(4, seed_nucleotide=nucleotide5),
                          ReportNucleotide(5, seed_nucleotide=nucleotide6)]

    region_start = 2
    prev_region_end = 5
    is_amino = True
    region_name = 'test-region'

    expected_counts = {0: nucleotide1, 1: nucleotide2, 2: nucleotide3, 3: nucleotide4, 4: nucleotide5, 5: nucleotide6}
    expected_log = [
        ('micall.core.aln2counts',
         logging.DEBUG,
         'Disagreement in counts for region test-region')]
    with caplog.at_level(logging.DEBUG):
        combine_region_nucleotides(nuc_dict, region_nucleotides, region_start, prev_region_end, is_amino, region_name)
    assert nuc_dict == expected_counts
    assert caplog.record_tuples == expected_log


# noinspection DuplicatedCode
def test_consensus_untranslated_region_differences(caplog):
    """ Test that the consensus is stitched together correctly for an untranslated region """
    nucleotide1 = SeedNucleotide()
    nucleotide2 = SeedNucleotide()
    nucleotide3 = SeedNucleotide()
    nucleotide4 = SeedNucleotide()
    nucleotide5 = SeedNucleotide()
    nucleotide6 = SeedNucleotide()
    nucleotide1.count_nucleotides('A', 3)
    nucleotide2.count_nucleotides('C', 3)
    nucleotide3.count_nucleotides('G', 3)
    nucleotide4.count_nucleotides('T', 3)
    nucleotide5.count_nucleotides('T', 6)
    nucleotide6.count_nucleotides('A', 6)
    nuc_dict = {0: nucleotide1, 1: nucleotide2, 2: nucleotide3, 3: nucleotide1, 4: nucleotide1}

    region_nucleotides = [ReportNucleotide(1, seed_nucleotide=nucleotide4),
                          ReportNucleotide(2, seed_nucleotide=nucleotide4),
                          ReportNucleotide(3, seed_nucleotide=nucleotide4),
                          ReportNucleotide(4, seed_nucleotide=nucleotide5),
                          ReportNucleotide(5, seed_nucleotide=nucleotide6)]

    region_start = 2
    prev_region_end = 5
    is_amino = False
    region_name = 'test-region'

    expected_counts = {0: nucleotide1, 1: nucleotide2, 2: nucleotide3, 3: nucleotide1, 4: nucleotide1, 5: nucleotide6}
    expected_log = [
        ('micall.core.aln2counts',
         logging.DEBUG,
         'Disagreement in counts for region test-region')]
    with caplog.at_level(logging.DEBUG):
        combine_region_nucleotides(nuc_dict, region_nucleotides, region_start, prev_region_end, is_amino, region_name)
    assert nuc_dict == expected_counts
    assert caplog.record_tuples == expected_log


# noinspection DuplicatedCode
def test_consensus_insertion_translated_differences():
    """ Test that the consensus insertions are handled correctly for translated regions """
    insertion1 = SeedNucleotide()
    insertion2 = SeedNucleotide()
    insertion3 = SeedNucleotide()
    insertion4 = SeedNucleotide()
    insertion1.count_nucleotides('A', 10)
    insertion2.count_nucleotides('C', 10)
    insertion3.count_nucleotides('G', 10)
    insertion4.count_nucleotides('T', 10)
    insertions_dict = {10: {1: insertion1}, 22: {1: insertion2}, 28: {1: insertion1}}
    region_insertions = {1: {1: insertion3}, 7: {1: insertion3}, 12: {1: insertion4}}
    region_start = 20
    prev_region_end = 29
    consumed_positions = set(i for i in range(0, 28))
    is_amino = True

    combine_region_insertions(insertions_dict,
                              region_insertions,
                              region_start,
                              prev_region_end,
                              is_amino,
                              consumed_positions)
    assert insertions_dict == {10: {1: insertion1}, 22: {1: insertion2}, 26: {1: insertion3}, 31: {1: insertion4}}


# noinspection DuplicatedCode
def test_consensus_insertion_untranslated_differences():
    """ Test that the consensus insertions are handled correctly for untranslated regions """
    insertion1 = SeedNucleotide()
    insertion2 = SeedNucleotide()
    insertion3 = SeedNucleotide()
    insertion4 = SeedNucleotide()
    insertion1.count_nucleotides('A', 10)
    insertion2.count_nucleotides('C', 10)
    insertion3.count_nucleotides('G', 10)
    insertion4.count_nucleotides('T', 10)
    insertions_dict = {10: {1: insertion1}, 22: {1: insertion2}, 28: {1: insertion1}}
    region_insertions = {1: {1: insertion3}, 7: {1: insertion3}, 12: {1: insertion4}}
    region_start = 20
    prev_region_end = 29
    consumed_positions = set(i for i in range(0, 28))
    is_amino = False

    combine_region_insertions(insertions_dict,
                              region_insertions,
                              region_start,
                              prev_region_end,
                              is_amino,
                              consumed_positions)
    assert insertions_dict == {10: {1: insertion1}, 22: {1: insertion2}, 28: {1: insertion1}, 31: {1: insertion4}}


def test_nucleotide_coordinates(default_sequence_report):
    # refname,qcut,rank,count,offset,seq
    aligned_reads = prepare_reads("""\
SARS-CoV-2-seed,15,0,9,0,ACGAACAAACT
""")

    expected_report = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,genome.pos,\
A,C,G,T,N,del,ins,clip,v3_overlap,coverage,exact_coverage
SARS-CoV-2-seed,SARS-CoV-2-TRS-B-8,15,1,1,28260,9,0,0,0,0,0,0,0,0,9,
SARS-CoV-2-seed,SARS-CoV-2-TRS-B-8,15,2,2,28261,0,9,0,0,0,0,0,0,0,9,
SARS-CoV-2-seed,SARS-CoV-2-TRS-B-8,15,3,3,28262,0,0,9,0,0,0,0,0,0,9,
SARS-CoV-2-seed,SARS-CoV-2-TRS-B-8,15,4,4,28263,9,0,0,0,0,0,0,0,0,9,
SARS-CoV-2-seed,SARS-CoV-2-TRS-B-8,15,5,5,28264,9,0,0,0,0,0,0,0,0,9,
SARS-CoV-2-seed,SARS-CoV-2-TRS-B-8,15,6,6,28265,0,9,0,0,0,0,0,0,0,9,
SARS-CoV-2-seed,SARS-CoV-2-TRS-B-8,15,7,7,28266,9,0,0,0,0,0,0,0,0,9,
SARS-CoV-2-seed,SARS-CoV-2-TRS-B-8,15,8,8,28267,9,0,0,0,0,0,0,0,0,9,
SARS-CoV-2-seed,SARS-CoV-2-TRS-B-8,15,9,9,28268,9,0,0,0,0,0,0,0,0,9,
SARS-CoV-2-seed,SARS-CoV-2-TRS-B-8,15,10,10,28269,0,9,0,0,0,0,0,0,0,9,
SARS-CoV-2-seed,SARS-CoV-2-TRS-B-8,15,11,11,28270,0,0,0,9,0,0,0,0,0,9,
"""

    report_file = StringIO()
    default_sequence_report.write_nuc_header(report_file)
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_nuc_counts()

    assert report_file.getvalue() == expected_report


def test_minimap_overlap(default_sequence_report, projects):
    """ Similar sections can fool minimap2 into reporting a section twice.

    In this example, positions 1-104 of the read map to pos 4401-4504 of the
    reference. Positions 101-200 of the read map to pos 3001-3100 of the ref.
    Even though positions 101-104 are in both alignments, only report them once.
    """
    seed_name = 'HIV1-B-FR-K03455-seed'
    seed_seq = projects.getReference(seed_name)
    read_seq = seed_seq[4440:4500] + seed_seq[3000:3060]

    # refname,qcut,rank,count,offset,seq
    aligned_reads = prepare_reads(f"""\
HIV1-B-FR-K03455-seed,15,0,9,0,{read_seq}
""")

    #                                    A,C,G,T
    expected_text = """\
HIV1-B-FR-K03455-seed,INT,15,51,262,4491,0,0,9,0,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,INT,15,52,263,4492,0,0,0,9,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,INT,15,53,264,4493,0,0,0,9,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,INT,15,54,265,4494,9,0,0,0,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,INT,15,55,266,4495,0,0,0,9,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,INT,15,56,267,4496,0,0,0,9,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,INT,15,57,268,4497,0,9,0,0,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,INT,15,58,269,4498,0,9,0,0,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,INT,15,59,270,4499,9,0,0,0,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,INT,15,60,271,4500,0,0,9,0,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,RT,15,61,452,3001,9,0,0,0,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,RT,15,62,453,3002,0,0,9,0,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,RT,15,63,454,3003,0,0,9,0,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,RT,15,64,455,3004,0,0,9,0,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,RT,15,65,456,3005,9,0,0,0,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,RT,15,66,457,3006,0,0,0,9,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,RT,15,67,458,3007,0,0,9,0,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,RT,15,68,459,3008,0,0,9,0,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,RT,15,69,460,3009,9,0,0,0,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,RT,15,70,461,3010,9,0,0,0,0,0,0,0,0,9,"""
    report_file = StringIO()
    default_sequence_report.write_nuc_header(report_file)
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_nuc_counts()

    report = report_file.getvalue()
    report_lines = report.splitlines()
    expected_size = 121
    if len(report_lines) != expected_size:
        assert (len(report_lines), report) == (expected_size, '')

    key_lines = report_lines[51:71]
    key_report = '\n'.join(key_lines)
    assert key_report == expected_text


def test_minimap_overlap_at_start(default_sequence_report, projects):
    """ Actual overlaps cause blank query position. Check consensus offset.

    In this example, the start of PR appears twice, so the consensus index
    gets blanked. Make sure that the PR consensus has the correct offset and
    doesn't crash.
    """
    seed_name = 'HIV1-B-FR-K03455-seed'
    seed_seq = projects.getReference(seed_name)
    read_seq = seed_seq[2252:2400] + seed_seq[2000:2258]

    # refname,qcut,rank,count,offset,seq
    aligned_reads = prepare_reads(f"""\
HIV1-B-FR-K03455-seed,15,0,9,0,{read_seq}
""")

    expected_text = f"""\
seed,region,q-cutoff,consensus-percent-cutoff,seed-offset,region-offset,sequence
HIV1-B-FR-K03455-seed,,15,MAX,0,,{read_seq}
HIV1-B-FR-K03455-seed,HIV1B-gag,15,MAX,0,1463,{read_seq}
HIV1-B-FR-K03455-seed,PR,15,MAX,0,0,{read_seq}
"""
    report_file = StringIO()
    default_sequence_report.write_consensus_all_header(report_file)
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_consensus_all()

    assert report_file.getvalue() == expected_text


def test_minimap_gap(default_sequence_report, projects):
    """ Large gap should still have coverage on either side, plus deletions. """
    seed_name = 'HIV1-B-FR-K03455-seed'
    seed_seq = projects.getReference(seed_name)
    read_seq = seed_seq[789:1282] + seed_seq[1861:2292]

    # refname,qcut,rank,count,offset,seq
    aligned_reads = prepare_reads(f"""\
HIV1-B-FR-K03455-seed,15,0,9,0,{read_seq}
""")
    #                                           A,C,G,T
    expected_text = """\
HIV1-B-FR-K03455-seed,HIV1B-gag,15,493,493,1282,9,0,0,0,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,,494,1283,0,0,0,0,0,9,0,0,0,9,
...
HIV1-B-FR-K03455-seed,HIV1B-gag,15,,1072,1861,0,0,0,0,0,9,0,0,0,9,
HIV1-B-FR-K03455-seed,HIV1B-gag,15,494,1073,1862,9,0,0,0,0,0,0,0,0,9,"""
    report_file = StringIO()
    default_sequence_report.write_nuc_header(report_file)
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_nuc_counts()

    report = report_file.getvalue()
    report_lines = report.splitlines()
    expected_size = 1558
    if len(report_lines) != expected_size:
        assert (len(report_lines), report) == (expected_size, '')

    key_lines = report_lines[493:495] + ['...'] + report_lines[1072:1074]
    key_report = '\n'.join(key_lines)
    assert key_report == expected_text


def test_minimap_gap_around_start(default_sequence_report, projects):
    """ Large gap surrounding GP41 end, nucleotide gap, and nef start.

    GP41 ends at 8795 gap is at 8796, and nef starts at 8797.
    """
    seed_name = 'HIV1-B-FR-K03455-seed'
    seed_seq = projects.getReference(seed_name)
    read_seq = seed_seq[8500:8700] + seed_seq[8900:9050]

    # refname,qcut,rank,count,offset,seq
    aligned_reads = prepare_reads(f"""\
HIV1-B-FR-K03455-seed,15,0,9,0,{read_seq}
""")
    expected_text = """\
HIV1-B-FR-K03455-seed,GP41,15,,1037,8794,0,0,0,0,0,9,0,0,0,9,
HIV1-B-FR-K03455-seed,GP41,15,,1038,8795,0,0,0,0,0,9,0,0,0,9,
HIV1-B-FR-K03455-seed,HIV1B-8796,15,,1,8796,0,0,0,0,0,9,0,0,0,9,
HIV1-B-FR-K03455-seed,HIV1B-nef,15,,1,8797,0,0,0,0,0,9,0,0,0,9,
HIV1-B-FR-K03455-seed,HIV1B-nef,15,,2,8798,0,0,0,0,0,9,0,0,0,9,"""
    report_file = StringIO()
    default_sequence_report.write_nuc_header(report_file)
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_nuc_counts()

    report = report_file.getvalue()
    report_lines = report.splitlines()
    expected_size = 549
    if len(report_lines) != expected_size:
        assert (len(report_lines), report) == (expected_size, '')

    key_lines = report_lines[294:299]
    key_report = '\n'.join(key_lines)
    assert key_report == expected_text


def test_minimap_reading_frame(default_sequence_report, projects):
    """ Deletion before PR shouldn't throw off the reading frame.

    PR runs from 2253 to 2549, so put a 2-base deletion before the start.
    """
    seed_name = 'HIV1-B-FR-K03455-seed'
    seed_seq = projects.getReference(seed_name)
    read_seq = seed_seq[2100:2200] + seed_seq[2202:2400]

    # refname,qcut,rank,count,offset,seq
    aligned_reads = prepare_reads(f"""\
HIV1-B-FR-K03455-seed,15,0,9,0,{read_seq}
""")
    #                                     A,C,G,T
    expected_text = """\
HIV1-B-FR-K03455-seed,HIV1B-gag,15,190,1503,2292,9,0,0,0,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,PR,15,151,1,2253,0,9,0,0,0,0,0,0,0,9,
HIV1-B-FR-K03455-seed,PR,15,152,2,2254,0,9,0,0,0,0,0,0,0,9,"""
    report_file = StringIO()
    default_sequence_report.write_nuc_header(report_file)
    default_sequence_report.read(aligned_reads)
    default_sequence_report.write_nuc_counts()

    report = report_file.getvalue()
    report_lines = report.splitlines()
    expected_size = 341
    if len(report_lines) != expected_size:
        assert (len(report_lines), report) == (expected_size, '')

    key_lines = report_lines[192:195]
    key_report = '\n'.join(key_lines)
    assert key_report == expected_text


# noinspection DuplicatedCode
def test_contig_coverage_report_huge_gap(default_sequence_report):
    """ A gap so big that Gotoh can't bridge it, but minimap2 can. """
    ref = default_sequence_report.projects.getReference('HIV1-B-FR-K03455-seed')
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
    default_sequence_report.read_remap_conseqs(remap_conseq_csv)
    default_sequence_report.write_amino_header(StringIO())
    default_sequence_report.write_genome_coverage_header(report_file)
    default_sequence_report.read(aligned_reads1)
    default_sequence_report.write_genome_coverage_counts()
    default_sequence_report.write_amino_counts()

    report_file.seek(0)
    covered_positions = [int(row['refseq_nuc_pos'])
                         for row in DictReader(report_file)
                         if row['refseq_nuc_pos']]
    assert covered_positions == expected_positions


# noinspection DuplicatedCode
def test_contig_coverage_report_past_reference_end(projects, sequence_report):
    hxb2_name = 'HIV1-B-FR-K03455-seed'
    ref = projects.getReference(hxb2_name)
    assert len(ref) == 9719
    seq = ref[-100:] + 'CGTAC'
    seed_nucs = [SeedNucleotide(Counter({'C': 1}))] * len(seq)
    expected_tail = """\
1-my-contig,HIV1-B-FR-K03455-seed,99,9718,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,100,9719,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,101,9720,0,1,,U
1-my-contig,HIV1-B-FR-K03455-seed,102,9721,0,1,,U
1-my-contig,HIV1-B-FR-K03455-seed,103,9722,0,1,,U
1-my-contig,HIV1-B-FR-K03455-seed,104,9723,0,1,,U
1-my-contig,HIV1-B-FR-K03455-seed,105,9724,0,1,,U
"""

    report_file = StringIO()
    sequence_report.projects = projects
    sequence_report.write_genome_coverage_header(report_file)
    sequence_report.write_sequence_coverage_counts('1-my-contig',
                                                   hxb2_name,
                                                   seq,
                                                   seed_nucs=seed_nucs)

    report_text = report_file.getvalue()
    tail = report_text[-len(expected_tail):]
    assert tail == expected_tail


# noinspection DuplicatedCode
def test_contig_coverage_report_past_reference_start(projects, sequence_report):
    hxb2_name = 'HIV1-B-FR-K03455-seed'
    ref = projects.getReference(hxb2_name)
    seq = 'CGTAC' + ref[:100]
    seed_nucs = [SeedNucleotide(Counter({'C': 1}))] * len(seq)
    # link is (M)apped, (U)nmapped, or (I)nserted
    expected_head = """\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage,concordance,link
1-my-contig,HIV1-B-FR-K03455-seed,1,-4,0,1,,U
1-my-contig,HIV1-B-FR-K03455-seed,2,-3,0,1,,U
1-my-contig,HIV1-B-FR-K03455-seed,3,-2,0,1,,U
1-my-contig,HIV1-B-FR-K03455-seed,4,-1,0,1,,U
1-my-contig,HIV1-B-FR-K03455-seed,5,0,0,1,,U
1-my-contig,HIV1-B-FR-K03455-seed,6,1,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,7,2,0,1,,M
"""

    report_file = StringIO()
    sequence_report.projects = projects
    sequence_report.write_genome_coverage_header(report_file)
    sequence_report.write_sequence_coverage_counts('1-my-contig',
                                                   hxb2_name,
                                                   seq,
                                                   seed_nucs=seed_nucs)

    report_text = report_file.getvalue()
    head = report_text[:len(expected_head)]
    assert head == expected_head


# noinspection DuplicatedCode
def test_contig_coverage_report_offset_reads(projects, sequence_report):
    hxb2_name = 'HIV1-B-FR-K03455-seed'
    ref = projects.getReference(hxb2_name)
    seq = ref[50:150]
    seed_nucs = ([SeedNucleotide()] * 50 +
                 [SeedNucleotide(Counter({'C': 1}))] * len(seq))
    expected_head = """\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage,concordance,link
1-my-contig,HIV1-B-FR-K03455-seed,51,51,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,52,52,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,53,53,0,1,,M
"""

    report_file = StringIO()
    sequence_report.projects = projects
    sequence_report.write_genome_coverage_header(report_file)
    sequence_report.write_sequence_coverage_counts('1-my-contig',
                                                   hxb2_name,
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
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage,concordance,link
1-R1-seed-partial,,1,1,0,5,,M
1-R1-seed-partial,,2,2,0,5,,M
1-R1-seed-partial,,3,3,0,5,,M
1-R1-seed-partial,,4,4,0,5,,M
1-R1-seed-partial,,5,5,0,5,,M
1-R1-seed-partial,,6,6,0,5,,M
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
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage,concordance,link
1-R1-seed-reversed,,1,1,0,5,,M
1-R1-seed-reversed,,2,2,0,5,,M
1-R1-seed-reversed,,3,3,0,5,,M
1-R1-seed-reversed,,4,4,0,5,,M
1-R1-seed-reversed,,5,5,0,5,,M
1-R1-seed-reversed,,6,6,0,5,,M
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
    R1 is KFR, and R2 is KFGPR.
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
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage,concordance,link
1_3-R1-seed,R1-seed,1,1,0,5,1.0,M
1_3-R1-seed,R1-seed,2,2,0,5,1.0,M
1_3-R1-seed,R1-seed,3,3,0,5,1.0,M
1_3-R1-seed,R1-seed,4,4,0,7,1.0,M
1_3-R1-seed,R1-seed,5,5,0,7,1.0,M
1_3-R1-seed,R1-seed,6,6,0,7,1.0,M
1_3-R1-seed,R1-seed,7,7,0,2,1.0,M
1_3-R1-seed,R1-seed,8,8,0,2,1.0,M
1_3-R1-seed,R1-seed,9,9,0,2,1.0,M
contig-1-R1-seed,R1-seed,1,1,,,1.0,M
contig-1-R1-seed,R1-seed,2,2,,,1.0,M
contig-1-R1-seed,R1-seed,3,3,,,1.0,M
contig-1-R1-seed,R1-seed,4,4,,,1.0,M
contig-1-R1-seed,R1-seed,5,5,,,1.0,M
contig-1-R1-seed,R1-seed,6,6,,,1.0,M
contig-3-R1-seed,R1-seed,1,4,,,1.0,M
contig-3-R1-seed,R1-seed,2,5,,,1.0,M
contig-3-R1-seed,R1-seed,3,6,,,1.0,M
contig-3-R1-seed,R1-seed,4,7,,,1.0,M
contig-3-R1-seed,R1-seed,5,8,,,1.0,M
contig-3-R1-seed,R1-seed,6,9,,,1.0,M
2-R2-seed,R2-seed,1,7,0,4,1.0,M
2-R2-seed,R2-seed,2,8,0,4,1.0,M
2-R2-seed,R2-seed,3,9,0,4,1.0,M
2-R2-seed,R2-seed,4,10,0,4,1.0,M
2-R2-seed,R2-seed,5,11,0,4,1.0,M
2-R2-seed,R2-seed,6,12,0,4,1.0,M
contig-2-R2-seed,R2-seed,1,7,,,1.0,M
contig-2-R2-seed,R2-seed,2,8,,,1.0,M
contig-2-R2-seed,R2-seed,3,9,,,1.0,M
contig-2-R2-seed,R2-seed,4,10,,,1.0,M
contig-2-R2-seed,R2-seed,5,11,,,1.0,M
contig-2-R2-seed,R2-seed,6,12,,,1.0,M
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
    seq = ref[100:150] + ref[1000:1050]
    expected_positions = list(range(101, 151)) + list(range(1001, 1051))

    report_file = StringIO()
    sequence_report.projects = projects
    sequence_report.write_genome_coverage_header(report_file)
    sequence_report.write_sequence_coverage_counts('1-my-contig',
                                                   hxb2_name,
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
    seq = ref[100:150] + ref[1000:1050]
    seed_nucs = [SeedNucleotide(Counter({'C': 1}))] * 100
    seed_nucs[2] = SeedNucleotide(Counter({'G': 4}))
    seed_nucs[98] = SeedNucleotide(Counter({'T': 5}))
    expected_head = """\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage,concordance,link
1-my-contig,HIV1-B-FR-K03455-seed,1,101,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,2,102,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,3,103,0,4,,M
1-my-contig,HIV1-B-FR-K03455-seed,4,104,0,1,,M
"""
    expected_tail = """\
1-my-contig,HIV1-B-FR-K03455-seed,98,1048,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,99,1049,0,5,,M
1-my-contig,HIV1-B-FR-K03455-seed,100,1050,0,1,,M
"""

    report_file = StringIO()
    sequence_report.projects = projects
    sequence_report.write_genome_coverage_header(report_file)
    sequence_report.write_sequence_coverage_counts('1-my-contig',
                                                   hxb2_name,
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
    seq = ref[100:110] + ref[115:160]
    seed_nucs = [SeedNucleotide(Counter({'C': 1}))] * len(seq)
    seed_nucs[12] = SeedNucleotide(Counter({'G': 4}))
    expected_head = """\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage,concordance,link
1-my-contig,HIV1-B-FR-K03455-seed,1,101,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,2,102,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,3,103,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,4,104,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,5,105,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,6,106,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,7,107,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,8,108,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,9,109,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,10,110,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,,111,,,,D
1-my-contig,HIV1-B-FR-K03455-seed,,112,,,,D
1-my-contig,HIV1-B-FR-K03455-seed,,113,,,,D
1-my-contig,HIV1-B-FR-K03455-seed,,114,,,,D
1-my-contig,HIV1-B-FR-K03455-seed,,115,,,,D
1-my-contig,HIV1-B-FR-K03455-seed,11,116,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,12,117,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,13,118,0,4,,M
1-my-contig,HIV1-B-FR-K03455-seed,14,119,0,1,,M
"""

    report_file = StringIO()
    sequence_report.projects = projects
    sequence_report.write_genome_coverage_header(report_file)
    sequence_report.write_sequence_coverage_counts('1-my-contig',
                                                   hxb2_name,
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
    seq = ref[100:150]
    seed_nucs = [SeedNucleotide(Counter({'C': 1}))] * len(seq)
    seed_nucs[5] = SeedNucleotide(Counter({'G': 4, '-': 2}))
    expected_head = """\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage,concordance,link
1-my-contig,HIV1-B-FR-K03455-seed,1,101,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,2,102,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,3,103,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,4,104,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,5,105,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,6,106,2,6,,M
1-my-contig,HIV1-B-FR-K03455-seed,7,107,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,8,108,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,9,109,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,10,110,0,1,,M
"""

    report_file = StringIO()
    sequence_report.projects = projects
    sequence_report.write_genome_coverage_header(report_file)
    sequence_report.write_sequence_coverage_counts('1-my-contig',
                                                   hxb2_name,
                                                   seq,
                                                   seed_nucs=seed_nucs)

    report_text = report_file.getvalue()
    head = report_text[:len(expected_head)]
    assert head == expected_head


# noinspection DuplicatedCode
def test_write_sequence_coverage_counts_with_insert(projects, sequence_report):
    hxb2_name = 'HIV1-B-FR-K03455-seed'
    ref = projects.getReference(hxb2_name)
    seq = ref[100:110] + 'ACTGA' + ref[110:160]
    seed_nucs = [SeedNucleotide(Counter({'C': 1}))] * len(seq)
    seed_nucs[12] = SeedNucleotide(Counter({'T': 4}))
    expected_head = """\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage,concordance,link
1-my-contig,HIV1-B-FR-K03455-seed,1,101,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,2,102,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,3,103,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,4,104,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,5,105,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,6,106,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,7,107,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,8,108,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,9,109,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,10,110,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,11,,0,1,,I
1-my-contig,HIV1-B-FR-K03455-seed,12,,0,1,,I
1-my-contig,HIV1-B-FR-K03455-seed,13,,0,4,,I
1-my-contig,HIV1-B-FR-K03455-seed,14,,0,1,,I
1-my-contig,HIV1-B-FR-K03455-seed,15,,0,1,,I
1-my-contig,HIV1-B-FR-K03455-seed,16,111,0,1,,M
"""

    report_file = StringIO()
    sequence_report.projects = projects
    sequence_report.write_genome_coverage_header(report_file)
    sequence_report.write_sequence_coverage_counts('1-my-contig',
                                                   hxb2_name,
                                                   seq,
                                                   seed_nucs=seed_nucs)

    report_text = report_file.getvalue()
    head = report_text[:len(expected_head)]
    assert head == expected_head


# noinspection DuplicatedCode
def test_write_sequence_coverage_counts_with_concordance(projects, sequence_report):
    hxb2_name = 'HIV1-B-FR-K03455-seed'
    ref = projects.getReference(hxb2_name)
    seq = ref[100:110] + 'ACTGA' + ref[110:160] + ref[180:250]
    seed_nucs = [SeedNucleotide(Counter({'C': 1}))] * len(seq)
    seed_nucs[12] = SeedNucleotide(Counter({'T': 4}))
    expected_head = """\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage,concordance,link
1-my-contig,HIV1-B-FR-K03455-seed,1,101,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,2,102,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,3,103,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,4,104,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,5,105,0,1,,M
1-my-contig,HIV1-B-FR-K03455-seed,6,106,0,1,0.95,M
1-my-contig,HIV1-B-FR-K03455-seed,7,107,0,1,0.9,M
1-my-contig,HIV1-B-FR-K03455-seed,8,108,0,1,0.85,M
1-my-contig,HIV1-B-FR-K03455-seed,9,109,0,1,0.8,M
1-my-contig,HIV1-B-FR-K03455-seed,10,110,0,1,0.7,M
1-my-contig,HIV1-B-FR-K03455-seed,11,,0,1,,I
1-my-contig,HIV1-B-FR-K03455-seed,12,,0,1,,I
1-my-contig,HIV1-B-FR-K03455-seed,13,,0,4,,I
1-my-contig,HIV1-B-FR-K03455-seed,14,,0,1,,I
1-my-contig,HIV1-B-FR-K03455-seed,15,,0,1,,I
1-my-contig,HIV1-B-FR-K03455-seed,16,111,0,1,0.7,M
1-my-contig,HIV1-B-FR-K03455-seed,17,112,0,1,0.7,M
1-my-contig,HIV1-B-FR-K03455-seed,18,113,0,1,0.7,M
1-my-contig,HIV1-B-FR-K03455-seed,19,114,0,1,0.7,M
1-my-contig,HIV1-B-FR-K03455-seed,20,115,0,1,0.7,M
"""

    expected_middle = """\
1-my-contig,HIV1-B-FR-K03455-seed,64,159,0,1,0.6,M
1-my-contig,HIV1-B-FR-K03455-seed,65,160,0,1,0.6,M
1-my-contig,HIV1-B-FR-K03455-seed,,161,,,,D
1-my-contig,HIV1-B-FR-K03455-seed,,162,,,,D
1-my-contig,HIV1-B-FR-K03455-seed,,163,,,,D
1-my-contig,HIV1-B-FR-K03455-seed,,164,,,,D
1-my-contig,HIV1-B-FR-K03455-seed,,165,,,,D
"""

    report_file = StringIO()
    sequence_report.projects = projects
    sequence_report.coord_concordance = [None] * 5 + [0.95, 0.9, 0.85, 0.8, 0.7] + [None] * 5 + [0.7] * 45 + [0.6] * 5
    sequence_report.write_genome_coverage_header(report_file)
    sequence_report.write_sequence_coverage_counts('1-my-contig',
                                                   hxb2_name,
                                                   seq,
                                                   seed_nucs=seed_nucs)

    report_text = report_file.getvalue()
    head = report_text[:len(expected_head)]
    middle = "".join(report_text.splitlines(keepends=True)[64:71])
    assert head == expected_head
    assert middle == expected_middle


# noinspection DuplicatedCode
def test_write_sequence_coverage_counts_with_unaligned_middle(projects,
                                                              sequence_report):
    """ The middle 100 bases are from a different reference.

    They get reported with query positions, but no reference positions.
    """
    hxb2_name = 'HIV1-B-FR-K03455-seed'
    ref = projects.getReference(hxb2_name)
    hcv_ref = projects.getReference('HCV-1a')
    seq = ref[:100] + hcv_ref[1000:1100] + ref[1000:1100]
    expected_ref_positions = (list(range(1, 101)) +
                              list(range(1001, 1101)))
    expected_query_positions = list(range(1, 301))

    report_file = StringIO()
    sequence_report.projects = projects
    sequence_report.write_genome_coverage_header(report_file)
    sequence_report.write_sequence_coverage_counts('1-my-contig',
                                                   hxb2_name,
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
    seq = ref[2858:2908] + ref[8187:8237]
    expected_ref_positions = (list(range(2859, 2916)) + list(range(8196, 8238)))
    expected_query_positions = list(range(1, len(seq)+1))

    report_file = StringIO()
    sequence_report.projects = projects
    sequence_report.write_genome_coverage_header(report_file)
    sequence_report.write_sequence_coverage_counts('1-my-contig',
                                                   hxb2_name,
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
def test_write_sequence_coverage_minimap_hits(projects, sequence_report):
    hxb2_name = 'HIV1-B-FR-K03455-seed'
    ref = projects.getReference(hxb2_name)
    seq = ref[1000:1100] + ref[2000:2100]
    expected_minimap_hits = """\
contig,ref_name,start,end,ref_start,ref_end
1-my-contig,HIV1-B-FR-K03455-seed,1,100,1001,1100
1-my-contig,HIV1-B-FR-K03455-seed,101,200,2001,2100
"""
    report_file = StringIO()
    sequence_report.projects = projects
    sequence_report.write_genome_coverage_header(StringIO())
    sequence_report.write_minimap_hits_header(report_file)
    sequence_report.write_sequence_coverage_counts('1-my-contig', hxb2_name, seq)

    assert report_file.getvalue() == expected_minimap_hits


# noinspection DuplicatedCode
def test_write_sequence_coverage_minimap_hits_reversed(projects, sequence_report):
    hxb2_name = 'HIV1-B-FR-K03455-seed'
    ref = projects.getReference(hxb2_name)
    seq = ref[1000:1100] + revcomp(ref[2000:2100])
    expected_minimap_hits = """\
contig,ref_name,start,end,ref_start,ref_end
1-my-contig,HIV1-B-FR-K03455-seed,1,100,1001,1100
1-my-contig,HIV1-B-FR-K03455-seed,101,200,2100,2001
"""
    report_file = StringIO()
    sequence_report.projects = projects
    sequence_report.write_genome_coverage_header(StringIO())
    sequence_report.write_minimap_hits_header(report_file)
    sequence_report.write_sequence_coverage_counts('1-my-contig', hxb2_name, seq)

    assert report_file.getvalue() == expected_minimap_hits
