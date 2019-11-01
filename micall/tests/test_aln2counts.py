import csv
from io import StringIO
import sys
import unittest

from micall.core.aln2counts import SequenceReport, InsertionWriter, SeedAmino, \
    ReportAmino
from micall.core import project_config


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


class SequenceReportTest(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        self.insertion_file = StringIO()
        insert_writer = InsertionWriter(
            insert_file=self.insertion_file)
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
        self.report = StubbedSequenceReport(insert_writer,
                                            projects,
                                            conseq_mixture_cutoffs,
                                            landmarks_yaml=landmarks_yaml)
        self.report_file = StringIO()
        self.detail_report_file = StringIO()

    @staticmethod
    def prepareReads(aligned_reads_text):
        full_text = "refname,qcut,rank,count,offset,seq\n" + aligned_reads_text
        dummy_file = StringIO(full_text)
        return csv.DictReader(dummy_file)

    def testEmptyAminoReport(self):
        expected_text = ""

        self.report.write_amino_counts(self.report_file)

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testEmptyNucReport(self):
        expected_text = ""

        self.report.write_nuc_counts(self.report_file)

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testConsensusFromSingleRead(self):
        """ In this sample, there is a single read with two codons.
        AAA -> K
        TTT -> F
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,9,0,AAATTT
""")
        expected_text = """\
region,q-cutoff,min_coverage,consensus-percent-cutoff,offset,sequence
R1-seed,15,0,MAX,0,AAATTT
R1-seed,15,0,0.100,0,AAATTT
"""

        self.report.write_consensus_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_consensus()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testConsensusFromTwoReads(self):
        """ The second read is out voted by the first one.
        CCC -> P
        GGG -> G
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,9,0,AAATTT
R1-seed,15,0,1,0,CCCGGG
""")
        expected_text = """\
region,q-cutoff,min_coverage,consensus-percent-cutoff,offset,sequence
R1-seed,15,0,MAX,0,AAATTT
R1-seed,15,0,0.100,0,MMMKKK
"""

        self.report.write_consensus_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_consensus()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testConsensusWithOffset(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,9,3,AAATTT
R1-seed,15,0,1,7,TTGGG
""")
        expected_text = """\
region,q-cutoff,min_coverage,consensus-percent-cutoff,offset,sequence
R1-seed,15,0,MAX,3,AAATTTGGG
R1-seed,15,0,0.100,3,AAATTTGGG
"""

        self.report.write_consensus_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_consensus()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testConsensusFromPartialContig(self):
        """ Contigs with the -partial suffix report consensus. """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
1-R2-seed-partial,15,0,9,0,AAATTT
""")
        expected_text = """\
region,q-cutoff,min_coverage,consensus-percent-cutoff,offset,sequence
1-R2-seed-partial,15,0,MAX,0,AAATTT
1-R2-seed-partial,15,0,0.100,0,AAATTT
"""

        self.report.write_consensus_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_consensus()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testConsensusLowQualitySections(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,9,3,NNNTTT
R1-seed,15,0,1,7,TTNGG
""")
        expected_text = """\
region,q-cutoff,min_coverage,consensus-percent-cutoff,offset,sequence
R1-seed,15,1,MAX,6,TTT-GG
R1-seed,15,1,0.100,6,TTT-GG
"""
        self.report.consensus_min_coverage = 1

        self.report.write_consensus_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_consensus()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testConsensusLowQuality(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,9,3,NNNNNN
R1-seed,15,0,1,7,NNNNN
""")
        expected_text = """\
region,q-cutoff,min_coverage,consensus-percent-cutoff,offset,sequence
"""
        self.report.consensus_min_coverage = 1

        self.report.write_consensus_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_consensus()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testConsensusLowCoverageInMiddle(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,9,0,AAATTTGGG
R1-seed,15,0,1,0,AAAT
R1-seed,15,0,1,6,GGG
""")
        expected_text = """\
region,q-cutoff,min_coverage,consensus-percent-cutoff,offset,sequence
R1-seed,15,10,MAX,0,AAAT--GGG
R1-seed,15,10,0.100,0,AAAT--GGG
"""
        self.report.consensus_min_coverage = 10

        self.report.write_consensus_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_consensus()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testConsensusLowCoverageAtStart(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,9,0,AAATTTGGG
R1-seed,15,0,1,4,TTGGG
""")
        expected_text = """\
region,q-cutoff,min_coverage,consensus-percent-cutoff,offset,sequence
R1-seed,15,10,MAX,4,TTGGG
R1-seed,15,10,0.100,4,TTGGG
"""
        self.report.consensus_min_coverage = 10

        self.report.write_consensus_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_consensus()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testConsensusLowCoverageAtEnd(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,9,0,AAATTTGGG
R1-seed,15,0,1,0,AAAT
""")
        expected_text = """\
region,q-cutoff,min_coverage,consensus-percent-cutoff,offset,sequence
R1-seed,15,10,MAX,0,AAAT
R1-seed,15,10,0.100,0,AAAT
"""
        self.report.consensus_min_coverage = 10

        self.report.write_consensus_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_consensus()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testConsensusMultipleCoverageLevels(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,9,0,AAATTTGGG
R1-seed,15,0,1,0,AAAT
""")
        expected_text = """\
region,q-cutoff,min_coverage,consensus-percent-cutoff,offset,sequence
R1-seed,15,10,MAX,0,AAAT
R1-seed,15,10,0.100,0,AAAT
R1-seed,15,1,MAX,0,AAATTTGGG
R1-seed,15,1,0.100,0,AAATTTGGG
"""
        self.report.consensus_min_coverages = (10, 1)

        self.report.write_consensus_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_consensus()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testSingleReadAminoReport(self):
        """ In this sample, there is a single read with two codons.
        AAA -> K
        TTT -> F
        The coordinate reference has three codons, so the third position is
        empty.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,9,0,AAATTT
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
"""

        self.report.write_amino_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_amino_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testMultiplePrefixAminoReport(self):
        """ Assemble counts from three contigs to two references.

        Contig 1-R1 AAATTT -> KF
        Contig 2-R2 GGCCCG -> GP
        Contig 3-R1 TTTAGG -> FR

        Contig 1 and 3 should combine into R1 with KFR.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads1 = self.prepareReads("1-R1-seed,15,0,5,0,AAATTT")
        aligned_reads2 = self.prepareReads("2-R2-seed,15,0,4,0,GGCCCG")
        aligned_reads3 = self.prepareReads("3-R1-seed,15,0,2,0,TTTAGG")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R1-seed,R1,15,,1,0,0,0,0,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5
R1-seed,R1,15,,2,0,0,0,0,7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,2
R2-seed,R2,15,,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,,3,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4
R2-seed,R2,15,,4,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4
R2-seed,R2,15,,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
"""

        expected_detail_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
1-R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5
1-R1-seed,R1,15,4,2,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5
1-R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
2-R2-seed,R2,15,,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
2-R2-seed,R2,15,,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
2-R2-seed,R2,15,1,3,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4
2-R2-seed,R2,15,4,4,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4
2-R2-seed,R2,15,,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
3-R1-seed,R1,15,,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
3-R1-seed,R1,15,1,2,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2
3-R1-seed,R1,15,4,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,2
"""

        self.report.write_amino_header(self.report_file)
        self.report.write_amino_detail_header(self.detail_report_file)
        self.report.read(aligned_reads1)
        self.report.write_amino_detail_counts()
        self.report.combine_reports()
        self.report.read(aligned_reads2)
        self.report.write_amino_detail_counts()
        self.report.combine_reports()
        self.report.read(aligned_reads3)
        self.report.write_amino_detail_counts()
        self.report.combine_reports()
        self.report.write_amino_counts()

        self.assertMultiLineEqual(expected_detail_text,
                                  self.detail_report_file.getvalue())
        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    # noinspection DuplicatedCode
    def testMultiplePrefixPartialDeletionAminoReport(self):
        """ Assemble counts from multiple contigs.

        Contig 1-R1 AAATTT -> KF
        Contig 2-R1 TT-AGG -> fR (partial deletion)
        Contig 3-R1 AAA---AGG -> K-R (full deletion)

        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads1 = self.prepareReads("1-R1-seed,15,0,5,0,AAATTT")
        aligned_reads2 = self.prepareReads("2-R1-seed,15,0,2,0,TT-AGG")
        aligned_reads3 = self.prepareReads("3-R1-seed,15,0,3,0,AAATTTAGG\n"
                                           "3-R1-seed,15,0,1,0,AAA---AGG")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R1-seed,R1,15,,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,,2,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,1,0,0,0,9
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,0,0,0,0,0,0,0,0,0,0,0,0,6
"""

        expected_detail_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
1-R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5
1-R1-seed,R1,15,4,2,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5
1-R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
2-R1-seed,R1,15,,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
2-R1-seed,R1,15,1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0
2-R1-seed,R1,15,4,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,2
3-R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4
3-R1-seed,R1,15,4,2,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,4
3-R1-seed,R1,15,7,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,4
"""

        self.report.write_amino_header(self.report_file)
        self.report.write_amino_detail_header(self.detail_report_file)
        self.report.read(aligned_reads1)
        self.report.write_amino_detail_counts()
        self.report.combine_reports()
        self.report.read(aligned_reads2)
        self.report.write_amino_detail_counts()
        self.report.combine_reports()
        self.report.read(aligned_reads3)
        self.report.write_amino_detail_counts()
        self.report.combine_reports()
        self.report.write_amino_counts()

        self.assertMultiLineEqual(expected_detail_text,
                                  self.detail_report_file.getvalue())
        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testMultiplePrefixNucleotideReport(self):
        """ Assemble counts from three contigs to two references.

        Contig 1-R1 AAATTT -> KF
        Contig 2-R2 GGCCCG -> GP
        Contig 3-R1 TTTAGG -> FR

        Contig 1 and 3 should combine into R1 with KFR.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads1 = self.prepareReads("1-R1-seed,15,0,5,0,AAATTT")
        aligned_reads2 = self.prepareReads("2-R2-seed,15,0,4,0,GGCCCG")
        aligned_reads3 = self.prepareReads("3-R1-seed,15,0,2,0,TTTAGG")

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
R2-seed,R2,15,,1,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,,2,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,,3,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,,4,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,,5,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,,6,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,,7,0,0,4,0,0,0,0,0,0,4
R2-seed,R2,15,,8,0,0,4,0,0,0,0,0,0,4
R2-seed,R2,15,,9,0,4,0,0,0,0,0,0,0,4
R2-seed,R2,15,,10,0,4,0,0,0,0,0,0,0,4
R2-seed,R2,15,,11,0,4,0,0,0,0,0,0,0,4
R2-seed,R2,15,,12,0,0,4,0,0,0,0,0,0,4
R2-seed,R2,15,,13,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,,14,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,,15,0,0,0,0,0,0,0,0,0,0
"""

        expected_detail_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,\
A,C,G,T,N,del,ins,clip,v3_overlap,coverage
1-R1-seed,R1,15,1,1,5,0,0,0,0,0,0,0,0,5
1-R1-seed,R1,15,2,2,5,0,0,0,0,0,0,0,0,5
1-R1-seed,R1,15,3,3,5,0,0,0,0,0,0,0,0,5
1-R1-seed,R1,15,4,4,0,0,0,5,0,0,0,0,0,5
1-R1-seed,R1,15,5,5,0,0,0,5,0,0,0,0,0,5
1-R1-seed,R1,15,6,6,0,0,0,5,0,0,0,0,0,5
1-R1-seed,R1,15,,7,0,0,0,0,0,0,0,0,0,0
1-R1-seed,R1,15,,8,0,0,0,0,0,0,0,0,0,0
1-R1-seed,R1,15,,9,0,0,0,0,0,0,0,0,0,0
2-R2-seed,R2,15,,1,0,0,0,0,0,0,0,0,0,0
2-R2-seed,R2,15,,2,0,0,0,0,0,0,0,0,0,0
2-R2-seed,R2,15,,3,0,0,0,0,0,0,0,0,0,0
2-R2-seed,R2,15,,4,0,0,0,0,0,0,0,0,0,0
2-R2-seed,R2,15,,5,0,0,0,0,0,0,0,0,0,0
2-R2-seed,R2,15,,6,0,0,0,0,0,0,0,0,0,0
2-R2-seed,R2,15,1,7,0,0,4,0,0,0,0,0,0,4
2-R2-seed,R2,15,2,8,0,0,4,0,0,0,0,0,0,4
2-R2-seed,R2,15,3,9,0,4,0,0,0,0,0,0,0,4
2-R2-seed,R2,15,4,10,0,4,0,0,0,0,0,0,0,4
2-R2-seed,R2,15,5,11,0,4,0,0,0,0,0,0,0,4
2-R2-seed,R2,15,6,12,0,0,4,0,0,0,0,0,0,4
2-R2-seed,R2,15,,13,0,0,0,0,0,0,0,0,0,0
2-R2-seed,R2,15,,14,0,0,0,0,0,0,0,0,0,0
2-R2-seed,R2,15,,15,0,0,0,0,0,0,0,0,0,0
3-R1-seed,R1,15,,1,0,0,0,0,0,0,0,0,0,0
3-R1-seed,R1,15,,2,0,0,0,0,0,0,0,0,0,0
3-R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0
3-R1-seed,R1,15,1,4,0,0,0,2,0,0,0,0,0,2
3-R1-seed,R1,15,2,5,0,0,0,2,0,0,0,0,0,2
3-R1-seed,R1,15,3,6,0,0,0,2,0,0,0,0,0,2
3-R1-seed,R1,15,4,7,2,0,0,0,0,0,0,0,0,2
3-R1-seed,R1,15,5,8,0,0,2,0,0,0,0,0,0,2
3-R1-seed,R1,15,6,9,0,0,2,0,0,0,0,0,0,2
"""

        self.report.write_nuc_header(self.report_file)
        self.report.write_nuc_detail_header(self.detail_report_file)
        self.report.read(aligned_reads1)
        self.report.write_nuc_detail_counts()
        self.report.combine_reports()
        self.report.read(aligned_reads2)
        self.report.write_nuc_detail_counts()
        self.report.combine_reports()
        self.report.read(aligned_reads3)
        self.report.write_nuc_detail_counts()
        self.report.combine_reports()
        self.report.write_nuc_counts()

        self.assertMultiLineEqual(expected_detail_text,
                                  self.detail_report_file.getvalue())
        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testContigCoverageReport(self):
        """ Assemble counts from three contigs to two references.

        Contig 1-R1 AAATTT -> KF
        Contig 2-R2 GGCCCG -> GP
        Contig 3-R1 TTTAGG -> FR

        Contig 1 and 3 should combine into R1 with KFR.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads1 = self.prepareReads("1-R1-seed,15,0,5,0,AAATTT")
        aligned_reads2 = self.prepareReads("2-R2-seed,15,0,4,0,GGCCCG")
        aligned_reads3 = self.prepareReads("3-R1-seed,15,0,2,0,TTTAGG")

        expected_text = """\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage
1-R1-seed,R1-seed,1,1,0,5
1-R1-seed,R1-seed,2,2,0,5
1-R1-seed,R1-seed,3,3,0,5
1-R1-seed,R1-seed,4,4,0,5
1-R1-seed,R1-seed,5,5,0,5
1-R1-seed,R1-seed,6,6,0,5
2-R2-seed,R2-seed,1,7,0,4
2-R2-seed,R2-seed,2,8,0,4
2-R2-seed,R2-seed,3,9,0,4
2-R2-seed,R2-seed,4,10,0,4
2-R2-seed,R2-seed,5,11,0,4
2-R2-seed,R2-seed,6,12,0,4
3-R1-seed,R1-seed,1,4,0,2
3-R1-seed,R1-seed,2,5,0,2
3-R1-seed,R1-seed,3,6,0,2
3-R1-seed,R1-seed,4,7,0,2
3-R1-seed,R1-seed,5,8,0,2
3-R1-seed,R1-seed,6,9,0,2
"""

        self.report.write_amino_header(StringIO())
        self.report.write_amino_detail_header(StringIO())
        self.report.write_genome_coverage_header(self.report_file)
        self.report.read(aligned_reads1)
        self.report.write_genome_coverage_counts()
        self.report.write_amino_detail_counts()
        self.report.read(aligned_reads2)
        self.report.write_genome_coverage_counts()
        self.report.write_amino_detail_counts()
        self.report.read(aligned_reads3)
        self.report.write_genome_coverage_counts()
        self.report.write_amino_detail_counts()
        self.report.write_amino_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    # noinspection DuplicatedCode
    def testContigCoverageReportDeletions(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads1 = self.prepareReads("1-R2-seed,15,0,4,0,GGC-CG")

        expected_text = """\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage
1-R2-seed,R2-seed,1,7,0,4
1-R2-seed,R2-seed,2,8,0,4
1-R2-seed,R2-seed,3,9,0,4
1-R2-seed,R2-seed,4,10,4,4
1-R2-seed,R2-seed,5,11,0,4
1-R2-seed,R2-seed,6,12,0,4
"""

        self.report.write_amino_header(StringIO())
        self.report.write_amino_detail_header(StringIO())
        self.report.write_genome_coverage_header(self.report_file)
        self.report.read(aligned_reads1)
        self.report.write_genome_coverage_counts()
        self.report.write_amino_detail_counts()
        self.report.write_amino_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    # noinspection DuplicatedCode
    def testContigCoverageReportGap(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads1 = self.prepareReads("1-R3-seed,15,0,4,0,AAATTTCAGCCACGAGAGCAT")
        expected_text = """\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage
1-R3-seed,R3-seed,1,1,0,4
1-R3-seed,R3-seed,2,2,0,4
1-R3-seed,R3-seed,3,3,0,4
1-R3-seed,R3-seed,4,4,0,4
1-R3-seed,R3-seed,5,5,0,4
1-R3-seed,R3-seed,6,6,0,4
1-R3-seed,R3-seed,7,7,0,4
1-R3-seed,R3-seed,8,8,0,4
1-R3-seed,R3-seed,9,9,0,4
1-R3-seed,R3-seed,10,13,0,4
1-R3-seed,R3-seed,11,14,0,4
1-R3-seed,R3-seed,12,15,0,4
1-R3-seed,R3-seed,13,16,0,4
1-R3-seed,R3-seed,14,17,0,4
1-R3-seed,R3-seed,15,18,0,4
1-R3-seed,R3-seed,16,19,0,4
1-R3-seed,R3-seed,17,20,0,4
1-R3-seed,R3-seed,18,21,0,4
1-R3-seed,R3-seed,19,22,0,4
1-R3-seed,R3-seed,20,23,0,4
1-R3-seed,R3-seed,21,24,0,4
"""

        self.report.write_amino_header(StringIO())
        self.report.write_amino_detail_header(StringIO())
        self.report.write_genome_coverage_header(self.report_file)
        self.report.read(aligned_reads1)
        self.report.write_genome_coverage_counts()
        self.report.write_amino_detail_counts()
        self.report.write_amino_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    # noinspection DuplicatedCode
    def testContigCoverageReportInsertions(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads1 = self.prepareReads(
            "1-R3-seed,15,0,4,0,AAATTTCAGACCACACCACGAGAGCAT")
        # insertion:                        ^^^

        expected_text = """\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage
1-R3-seed,R3-seed,1,1,0,4
1-R3-seed,R3-seed,2,2,0,4
1-R3-seed,R3-seed,3,3,0,4
1-R3-seed,R3-seed,4,4,0,4
1-R3-seed,R3-seed,5,5,0,4
1-R3-seed,R3-seed,6,6,0,4
1-R3-seed,R3-seed,7,7,0,4
1-R3-seed,R3-seed,8,8,0,4
1-R3-seed,R3-seed,9,9,0,4
1-R3-seed,R3-seed,10,10,0,4
1-R3-seed,R3-seed,11,11,0,4
1-R3-seed,R3-seed,12,12,0,4
1-R3-seed,R3-seed,13,,0,4
1-R3-seed,R3-seed,14,,0,4
1-R3-seed,R3-seed,15,,0,4
1-R3-seed,R3-seed,16,13,0,4
1-R3-seed,R3-seed,17,14,0,4
1-R3-seed,R3-seed,18,15,0,4
1-R3-seed,R3-seed,19,16,0,4
1-R3-seed,R3-seed,20,17,0,4
1-R3-seed,R3-seed,21,18,0,4
1-R3-seed,R3-seed,22,19,0,4
1-R3-seed,R3-seed,23,20,0,4
1-R3-seed,R3-seed,24,21,0,4
1-R3-seed,R3-seed,25,22,0,4
1-R3-seed,R3-seed,26,23,0,4
1-R3-seed,R3-seed,27,24,0,4
"""

        self.report.write_amino_header(StringIO())
        self.report.write_amino_detail_header(StringIO())
        self.report.write_genome_coverage_header(self.report_file)
        self.report.read(aligned_reads1)
        self.report.write_genome_coverage_counts()
        self.report.write_amino_detail_counts()
        self.report.write_amino_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    # noinspection DuplicatedCode
    def testContigCoverageReportPastReferenceEnd(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads1 = self.prepareReads("1-R1-seed,15,0,4,0,AAATTTAGGGAGCAT")
        expected_text = """\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage
1-R1-seed,R1-seed,1,1,0,4
1-R1-seed,R1-seed,2,2,0,4
1-R1-seed,R1-seed,3,3,0,4
1-R1-seed,R1-seed,4,4,0,4
1-R1-seed,R1-seed,5,5,0,4
1-R1-seed,R1-seed,6,6,0,4
1-R1-seed,R1-seed,7,7,0,4
1-R1-seed,R1-seed,8,8,0,4
1-R1-seed,R1-seed,9,9,0,4
1-R1-seed,R1-seed,10,10,0,4
1-R1-seed,R1-seed,11,11,0,4
1-R1-seed,R1-seed,12,12,0,4
1-R1-seed,R1-seed,13,13,0,4
1-R1-seed,R1-seed,14,14,0,4
1-R1-seed,R1-seed,15,15,0,4
"""

        self.report.write_amino_header(StringIO())
        self.report.write_amino_detail_header(StringIO())
        self.report.write_genome_coverage_header(self.report_file)
        self.report.read(aligned_reads1)
        self.report.write_genome_coverage_counts()
        self.report.write_amino_detail_counts()
        self.report.write_amino_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    # noinspection DuplicatedCode
    def testContigCoverageReportPastReferenceStart(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads1 = self.prepareReads("1-R1-seed,15,0,4,0,GAGCATAAATTTAGG")
        expected_text = """\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage
1-R1-seed,R1-seed,1,-5,0,4
1-R1-seed,R1-seed,2,-4,0,4
1-R1-seed,R1-seed,3,-3,0,4
1-R1-seed,R1-seed,4,-2,0,4
1-R1-seed,R1-seed,5,-1,0,4
1-R1-seed,R1-seed,6,0,0,4
1-R1-seed,R1-seed,7,1,0,4
1-R1-seed,R1-seed,8,2,0,4
1-R1-seed,R1-seed,9,3,0,4
1-R1-seed,R1-seed,10,4,0,4
1-R1-seed,R1-seed,11,5,0,4
1-R1-seed,R1-seed,12,6,0,4
1-R1-seed,R1-seed,13,7,0,4
1-R1-seed,R1-seed,14,8,0,4
1-R1-seed,R1-seed,15,9,0,4
"""

        self.report.write_amino_header(StringIO())
        self.report.write_amino_detail_header(StringIO())
        self.report.write_genome_coverage_header(self.report_file)
        self.report.read(aligned_reads1)
        self.report.write_genome_coverage_counts()
        self.report.write_amino_detail_counts()
        self.report.write_amino_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    # noinspection DuplicatedCode
    def testContigCoverageReportOffsetReads(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads1 = self.prepareReads("1-R1-seed,15,0,4,10,AAATTTAGG")
        expected_text = """\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage
1-R1-seed,R1-seed,11,1,0,4
1-R1-seed,R1-seed,12,2,0,4
1-R1-seed,R1-seed,13,3,0,4
1-R1-seed,R1-seed,14,4,0,4
1-R1-seed,R1-seed,15,5,0,4
1-R1-seed,R1-seed,16,6,0,4
1-R1-seed,R1-seed,17,7,0,4
1-R1-seed,R1-seed,18,8,0,4
1-R1-seed,R1-seed,19,9,0,4
"""

        self.report.write_amino_header(StringIO())
        self.report.write_amino_detail_header(StringIO())
        self.report.write_genome_coverage_header(self.report_file)
        self.report.read(aligned_reads1)
        self.report.write_genome_coverage_counts()
        self.report.write_amino_detail_counts()
        self.report.write_amino_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    # noinspection DuplicatedCode
    def testContigCoverageReportLargeGap(self):
        """ When mapping against known references, there can be large gaps. """
        remap_conseq_csv = StringIO("""\
region,sequence
R3-seed,AAATTTCAGACCCCACGAGAGCAT
""")
        # refname,qcut,rank,count,offset,seq
        aligned_reads1 = self.prepareReads("""\
R3-seed,15,0,4,3,TTT
R3-seed,15,0,5,21,CAT
""")
        expected_text = """\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage
R3-seed,R3-seed,4,4,0,4
R3-seed,R3-seed,5,5,0,4
R3-seed,R3-seed,6,6,0,4
R3-seed,R3-seed,22,22,0,5
R3-seed,R3-seed,23,23,0,5
R3-seed,R3-seed,24,24,0,5
"""

        self.report.read_remap_conseqs(remap_conseq_csv)
        self.report.write_amino_header(StringIO())
        self.report.write_genome_coverage_header(self.report_file)
        self.report.read(aligned_reads1)
        self.report.write_genome_coverage_counts()
        self.report.write_amino_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testContigCoverageReportForPartialContig(self):
        """ Contig coverage is reported for partial contigs.

        Reference columns are left blank, though, because they're not aligned.
        See blast.csv for best guess at alignment.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads1 = self.prepareReads("1-R1-seed-partial,15,0,5,0,CCCCCC")
        contigs_csv = StringIO("""\
ref,match,group_ref,contig
R1-seed,1,R1-seed,CCCCCC
""")

        expected_text = """\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage
1-R1-seed-partial,,1,,0,5
1-R1-seed-partial,,2,,0,5
1-R1-seed-partial,,3,,0,5
1-R1-seed-partial,,4,,0,5
1-R1-seed-partial,,5,,0,5
1-R1-seed-partial,,6,,0,5
"""

        self.report.write_amino_header(StringIO())
        self.report.write_amino_detail_header(StringIO())
        self.report.read_contigs(contigs_csv)
        self.report.write_genome_coverage_header(self.report_file)
        self.report.read(aligned_reads1)
        self.report.write_genome_coverage_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testContigCoverageReportForReversedContig(self):
        """ Contig coverage is reported for reversed contigs.

        Reference columns are left blank, though, because they're not aligned.
        See blast.csv for best guess at alignment.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads1 = self.prepareReads("1-R1-seed-reversed,15,0,5,0,CCCCCC")
        contigs_csv = StringIO("""\
ref,match,group_ref,contig
R1-seed,1,R1-seed,CCCCCC
""")

        expected_text = """\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage
1-R1-seed-reversed,,1,,0,5
1-R1-seed-reversed,,2,,0,5
1-R1-seed-reversed,,3,,0,5
1-R1-seed-reversed,,4,,0,5
1-R1-seed-reversed,,5,,0,5
1-R1-seed-reversed,,6,,0,5
"""

        self.report.read_contigs(contigs_csv)
        self.report.write_genome_coverage_header(self.report_file)
        self.report.read(aligned_reads1)
        self.report.write_genome_coverage_counts()
        self.report.write_amino_detail_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testContigCoverageReportMergedContigs(self):
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
        aligned_reads1 = self.prepareReads("1_3-R1-seed,15,0,5,0,AAATTT\n"
                                           "1_3-R1-seed,15,0,2,3,TTTAGG")
        aligned_reads2 = self.prepareReads("2-R2-seed,15,0,4,0,GGCCCG")
        contigs_csv = StringIO("""\
ref,match,group_ref,contig
R1-seed,1,R1-seed,AAATTT
R2-seed,1,R2-seed,GGCCCG
R1-seed,1,R1-seed,TTTAGG
""")

        expected_text = """\
contig,coordinates,query_nuc_pos,refseq_nuc_pos,dels,coverage
1_3-R1-seed,R1-seed,1,1,0,5
1_3-R1-seed,R1-seed,2,2,0,5
1_3-R1-seed,R1-seed,3,3,0,5
1_3-R1-seed,R1-seed,4,4,0,7
1_3-R1-seed,R1-seed,5,5,0,7
1_3-R1-seed,R1-seed,6,6,0,7
1_3-R1-seed,R1-seed,7,7,0,2
1_3-R1-seed,R1-seed,8,8,0,2
1_3-R1-seed,R1-seed,9,9,0,2
contig-1-R1-seed,R1-seed,1,1,,
contig-1-R1-seed,R1-seed,2,2,,
contig-1-R1-seed,R1-seed,3,3,,
contig-1-R1-seed,R1-seed,4,4,,
contig-1-R1-seed,R1-seed,5,5,,
contig-1-R1-seed,R1-seed,6,6,,
contig-3-R1-seed,R1-seed,1,4,,
contig-3-R1-seed,R1-seed,2,5,,
contig-3-R1-seed,R1-seed,3,6,,
contig-3-R1-seed,R1-seed,4,7,,
contig-3-R1-seed,R1-seed,5,8,,
contig-3-R1-seed,R1-seed,6,9,,
2-R2-seed,R2-seed,1,7,0,4
2-R2-seed,R2-seed,2,8,0,4
2-R2-seed,R2-seed,3,9,0,4
2-R2-seed,R2-seed,4,10,0,4
2-R2-seed,R2-seed,5,11,0,4
2-R2-seed,R2-seed,6,12,0,4
contig-2-R2-seed,R2-seed,1,7,,
contig-2-R2-seed,R2-seed,2,8,,
contig-2-R2-seed,R2-seed,3,9,,
contig-2-R2-seed,R2-seed,4,10,,
contig-2-R2-seed,R2-seed,5,11,,
contig-2-R2-seed,R2-seed,6,12,,
"""

        self.report.read_contigs(contigs_csv)
        self.report.write_amino_header(StringIO())
        self.report.write_amino_detail_header(StringIO())
        self.report.write_genome_coverage_header(self.report_file)
        self.report.read(aligned_reads1)
        self.report.write_genome_coverage_counts()
        self.report.write_amino_counts()
        self.report.read(aligned_reads2)
        self.report.write_genome_coverage_counts()
        self.report.write_amino_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testAminoReportForPartialContig(self):
        """ Contigs with the -partial suffix shouldn't be reported. """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
1-R1-seed-partial,15,0,9,0,AAATTT
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
"""

        self.report.write_amino_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_amino_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testSingleReadNucleotideReport(self):
        """ In this sample, there is a single read with two codons.
        AAA -> K
        TTT -> F
        The coordinate reference has three codons, so the third position is
        empty.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
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

        self.report.write_nuc_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_nuc_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testSoftClippingNucleotideReport(self):
        """ Combine the soft clipping data with the read counts.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,9,2,ATTTA
""")
        clipping = StringIO("""\
refname,pos,count
R1-seed,1,9
R1-seed,2,9
R1-seed,8,9
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T,N,del,ins,clip,v3_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,9,0,0
R1-seed,R1,15,2,2,0,0,0,0,0,0,0,9,0,0
R1-seed,R1,15,3,3,9,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,4,4,0,0,0,9,0,0,0,0,0,9
R1-seed,R1,15,5,5,0,0,0,9,0,0,0,0,0,9
R1-seed,R1,15,6,6,0,0,0,9,0,0,0,0,0,9
R1-seed,R1,15,7,7,9,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,8,8,0,0,0,0,0,0,0,9,0,0
R1-seed,R1,15,9,9,0,0,0,0,0,0,0,0,0,0
"""

        self.report.read_clipping(clipping)
        self.report.write_nuc_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_nuc_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testSoftClippingAminoReport(self):
        """ Combine the soft clipping data with the read counts.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,9,2,ATTTA
""")
        clipping = StringIO("""\
refname,pos,count
R1-seed,1,9
R1-seed,2,9
R1-seed,8,9
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0
R1-seed,R1,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,7,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0
"""

        self.report.read_clipping(clipping)
        self.report.write_amino_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_nuc_header(StringIO())
        self.report.write_nuc_counts()
        self.report.write_amino_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testMultiplePrefixSoftClippingAminoReport(self):
        """ Combine the soft clipping data with the read counts.
        """
        """ Assemble counts from three contigs to two references.
        Contig 1-R1 AAATTT -> KF
        Contig 2-R2 GGCCCG -> GP
        Contig 3-R1 TTTAGG -> FR
        Contig 1 and 3 should combine into R1 with KFR.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads1 = self.prepareReads("1-R1-seed,15,0,5,0,AAATTT")
        aligned_reads2 = self.prepareReads("2-R2-seed,15,0,4,0,GGCCCG")
        aligned_reads3 = self.prepareReads("3-R1-seed,15,0,2,0,TTTAGG")

        clipping = StringIO("""\
refname,pos,count
1-R1-seed,7,5
3-R1-seed,-1,2
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R1-seed,R1,15,,1,0,0,0,0,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,5
R1-seed,R1,15,,2,0,0,0,0,7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,5,0,2
R2-seed,R2,15,,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,,3,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4
R2-seed,R2,15,,4,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4
R2-seed,R2,15,,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
"""

        self.report.read_clipping(clipping)
        self.report.write_amino_header(self.report_file)
        self.report.write_amino_detail_header(self.detail_report_file)
        self.report.write_nuc_header(StringIO())
        self.report.read(aligned_reads1)
        self.report.write_nuc_counts()
        self.report.write_amino_detail_counts()
        self.report.combine_reports()
        self.report.read(aligned_reads2)
        self.report.write_nuc_counts()
        self.report.write_amino_detail_counts()
        self.report.combine_reports()
        self.report.read(aligned_reads3)
        self.report.write_nuc_counts()
        self.report.write_amino_detail_counts()
        self.report.combine_reports()
        self.report.write_amino_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testInsertionBetweenReadAndConsensusNucleotideReport(self):
        """ Combine the soft clipping data with the read counts.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,9,0,AAATTT
""")
        conseq_ins_csv = StringIO("""\
qname,fwd_rev,refname,pos,insert,qual
Example_read_1,F,R1-seed,3,AAC,AAA
Example_read_2,F,R1-seed,3,AAC,AAA
Example_read_2,R,R1-seed,3,AAC,AAA
Example_read_3,F,R2-seed,6,GTA,AAA
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T,N,del,ins,clip,v3_overlap,coverage
R1-seed,R1,15,1,1,9,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,2,2,9,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,3,3,9,0,0,0,0,0,2,0,0,9
R1-seed,R1,15,4,4,0,0,0,9,0,0,0,0,0,9
R1-seed,R1,15,5,5,0,0,0,9,0,0,0,0,0,9
R1-seed,R1,15,6,6,0,0,0,9,0,0,0,0,0,9
R1-seed,R1,15,,7,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,,8,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,,9,0,0,0,0,0,0,0,0,0,0
"""

        self.report.read_insertions(conseq_ins_csv)
        self.report.write_nuc_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_nuc_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    # noinspection DuplicatedCode
    def testInsertionBetweenReadAndConsensusAminoReport(self):
        """ Combine the soft clipping data with the read counts.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,9,0,AAATTT
""")
        conseq_ins_csv = StringIO("""\
qname,fwd_rev,refname,pos,insert,qual
Example_read_1,F,R1-seed,3,AAC,AAA
Example_read_2,F,R1-seed,3,AAC,AAA
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,9
R1-seed,R1,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
"""

        self.report.read_insertions(conseq_ins_csv)
        self.report.write_amino_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_nuc_header(StringIO())
        self.report.write_nuc_counts()  # calculates ins counts
        self.report.write_amino_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testSubstitutionAtBoundary(self):
        """ In this sample, there are nine identical reads with six codons.
        ATG -> M
        GCA -> A
        AAC -> N
        TGG -> W
        ATC -> I
        AAT -> N
        The R4 coordinate reference is SIN, so its first position will not map.
        However, the R4-seed reference matches the BAD and the IN, so the first
        position should get treated as a substitution.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R4-seed,15,0,9,0,ATGGCAAACTGGATCAAT
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R4-seed,R4,15,10,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,9
R4-seed,R4,15,13,2,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R4-seed,R4,15,16,3,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
"""

        self.report.write_amino_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_amino_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testCoverageSummary(self):
        """ R1 has coverage 9.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,9,0,AAATTTCGA
""")
        expected_summary = dict(avg_coverage=9.0,
                                coverage_region='R1',
                                region_width=3)

        summary = {}
        self.report.read(aligned_reads)
        self.report.write_amino_header(StringIO())
        self.report.write_amino_counts(coverage_summary=summary)

        self.assertEqual(expected_summary, summary)

    def testCoverageSummaryNotImproved(self):
        """ R2 has coverage 9, and R1 had coverage 50. Report R1.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R2-seed,15,0,9,0,AAATTTCGA
""")
        expected_summary = dict(avg_coverage=50.0,
                                coverage_region='R1',
                                region_width=3)

        summary = dict(expected_summary)
        self.report.read(aligned_reads)
        self.report.write_amino_header(StringIO())
        self.report.write_amino_counts(coverage_summary=summary)

        self.assertEqual(expected_summary, summary)

    def testCoverageSummaryNoCoverage(self):
        """ Stuff mapped to the seed, but didn't align with the coordinate region.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,9,10,TGGTGGTGG
""")
        expected_summary = {}

        summary = {}
        self.report.read(aligned_reads)
        self.report.write_amino_header(StringIO())
        self.report.write_amino_counts(coverage_summary=summary)

        self.assertEqual(expected_summary, summary)

    def testOffsetNucleotideReport(self):
        """ The first row provides alignment so the partial codon at the start
        of the second row will map to the reference.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,1,3,TTT
R1-seed,15,0,8,5,TCGA
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T,N,del,ins,clip,v3_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,2,2,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,3,3,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,4,4,0,0,0,1,0,0,0,0,0,1
R1-seed,R1,15,5,5,0,0,0,1,0,0,0,0,0,1
R1-seed,R1,15,6,6,0,0,0,9,0,0,0,0,0,9
R1-seed,R1,15,7,7,0,8,0,0,0,0,0,0,0,8
R1-seed,R1,15,8,8,0,0,8,0,0,0,0,0,0,8
R1-seed,R1,15,9,9,8,0,0,0,0,0,0,0,0,8
"""

        self.report.read(aligned_reads)
        self.report.write_nuc_header(self.report_file)
        self.report.write_nuc_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testPartialCodonNucleotideReport(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,9,0,AAATT
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T,N,del,ins,clip,v3_overlap,coverage
R1-seed,R1,15,1,1,9,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,2,2,9,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,3,3,9,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,4,4,0,0,0,9,0,0,0,0,0,9
R1-seed,R1,15,5,5,0,0,0,9,0,0,0,0,0,9
R1-seed,R1,15,6,6,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,,7,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,,8,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,,9,0,0,0,0,0,0,0,0,0,0
"""

        self.report.read(aligned_reads)
        self.report.write_nuc_header(self.report_file)
        self.report.write_nuc_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testPartialStartCodonNucleotideReport(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,9,0,TTAGG
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,\
A,C,G,T,N,del,ins,clip,v3_overlap,coverage
R1-seed,R1,15,,1,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,,2,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,0,4,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,1,5,0,0,0,9,0,0,0,0,0,9
R1-seed,R1,15,2,6,0,0,0,9,0,0,0,0,0,9
R1-seed,R1,15,3,7,9,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,4,8,0,0,9,0,0,0,0,0,0,9
R1-seed,R1,15,5,9,0,0,9,0,0,0,0,0,0,9
"""

        self.report.read(aligned_reads)
        self.report.write_nuc_header(self.report_file)
        self.report.write_nuc_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testLowQualityNucleotideReport(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,9,0,AAATNT
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T,N,del,ins,clip,v3_overlap,coverage
R1-seed,R1,15,1,1,9,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,2,2,9,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,3,3,9,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,4,4,0,0,0,9,0,0,0,0,0,9
R1-seed,R1,15,5,5,0,0,0,0,9,0,0,0,0,0
R1-seed,R1,15,6,6,0,0,0,9,0,0,0,0,0,9
R1-seed,R1,15,,7,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,,8,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,,9,0,0,0,0,0,0,0,0,0,0
"""

        self.report.read(aligned_reads)
        self.report.write_nuc_header(self.report_file)
        self.report.write_nuc_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testLowQualityAminoReport(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,9,0,AAATNT
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,4,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
"""

        self.report.read(aligned_reads)
        self.report.write_amino_header(self.report_file)
        self.report.write_amino_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testPartialDeletionAminoReport(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,9,0,AAAT-T
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,4,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
"""

        self.report.read(aligned_reads)
        self.report.write_amino_header(self.report_file)
        self.report.write_amino_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testShiftedReadingFrameAminoReport(self):
        """ The seed's reading frame doesn't match the coordinate reference's
        reading frame, so there is an extra nucleotide at the beginning of the
        reads.
        It will try padding the first codon to see which of the three possible
        reading frames gives the highest alignment score.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,9,0,GAAATTTCGA
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R1-seed,R1,15,2,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,5,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,8,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,9
"""

        self.report.read(aligned_reads)
        self.report.write_amino_header(self.report_file)
        self.report.write_amino_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testShiftedReadingFrameNucleotideReport(self):
        """ The seed's reading frame doesn't match the coordinate reference's
        reading frame, so there is an extra nucleotide at the beginning of the
        reads.
        It will try padding the first codon to see which of the three possible
        reading frames gives the highest alignment score.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,9,0,GAAATTTCGA
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T,N,del,ins,clip,v3_overlap,coverage
R1-seed,R1,15,2,1,9,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,3,2,9,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,4,3,9,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,5,4,0,0,0,9,0,0,0,0,0,9
R1-seed,R1,15,6,5,0,0,0,9,0,0,0,0,0,9
R1-seed,R1,15,7,6,0,0,0,9,0,0,0,0,0,9
R1-seed,R1,15,8,7,0,9,0,0,0,0,0,0,0,9
R1-seed,R1,15,9,8,0,0,9,0,0,0,0,0,0,9
R1-seed,R1,15,10,9,9,0,0,0,0,0,0,0,0,9
"""

        self.report.read(aligned_reads)
        self.report.write_nuc_header(self.report_file)
        self.report.write_nuc_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testDeletionNucleotideReport(self):
        """ Coordinate sequence is KFGPR, and this aligned read is KFPR.

        Must be a deletion in the seed reference with respect to the coordinate
        reference.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,9,0,AAA---AGG
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T,N,del,ins,clip,v3_overlap,coverage
R1-seed,R1,15,1,1,9,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,2,2,9,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,3,3,9,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,4,4,0,0,0,0,0,9,0,0,0,9
R1-seed,R1,15,5,5,0,0,0,0,0,9,0,0,0,9
R1-seed,R1,15,6,6,0,0,0,0,0,9,0,0,0,9
R1-seed,R1,15,7,7,9,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,8,8,0,0,9,0,0,0,0,0,0,9
R1-seed,R1,15,9,9,0,0,9,0,0,0,0,0,0,9
"""

        self.report.read(aligned_reads)
        self.report.write_nuc_header(self.report_file)
        self.report.write_nuc_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testDeletionBetweenSeedAndCoordinateNucleotideReport(self):
        """ Coordinate sequence is KFGPR, and this aligned read is KFPR.

        Must be a deletion in the seed reference with respect to the coordinate
        reference.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R2-seed,15,0,9,0,AAATTTCCCCGA
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T,N,del,ins,clip,v3_overlap,coverage
R2-seed,R2,15,1,1,9,0,0,0,0,0,0,0,0,9
R2-seed,R2,15,2,2,9,0,0,0,0,0,0,0,0,9
R2-seed,R2,15,3,3,9,0,0,0,0,0,0,0,0,9
R2-seed,R2,15,4,4,0,0,0,9,0,0,0,0,0,9
R2-seed,R2,15,5,5,0,0,0,9,0,0,0,0,0,9
R2-seed,R2,15,6,6,0,0,0,9,0,0,0,0,0,9
R2-seed,R2,15,,7,0,0,0,0,0,9,0,0,0,9
R2-seed,R2,15,,8,0,0,0,0,0,9,0,0,0,9
R2-seed,R2,15,,9,0,0,0,0,0,9,0,0,0,9
R2-seed,R2,15,7,10,0,9,0,0,0,0,0,0,0,9
R2-seed,R2,15,8,11,0,9,0,0,0,0,0,0,0,9
R2-seed,R2,15,9,12,0,9,0,0,0,0,0,0,0,9
R2-seed,R2,15,10,13,0,9,0,0,0,0,0,0,0,9
R2-seed,R2,15,11,14,0,0,9,0,0,0,0,0,0,9
R2-seed,R2,15,12,15,9,0,0,0,0,0,0,0,0,9
"""

        self.report.read(aligned_reads)
        self.report.write_nuc_header(self.report_file)
        self.report.write_nuc_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testDeletionBetweenSeedAndCoordinateAminoReport(self):
        """ Coordinate sequence is KFGPR, and this aligned read is KFPR.

        Must be a deletion in the seed reference with respect to the coordinate
        reference.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R2-seed,15,0,9,0,AAATTTCCCCGA
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R2-seed,R2,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R2-seed,R2,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R2-seed,R2,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,9
R2-seed,R2,15,7,4,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R2-seed,R2,15,10,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,9
"""

        self.report.read(aligned_reads)
        self.report.write_amino_header(self.report_file)
        self.report.write_amino_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testDeletionBetweenSeedAndConsensusAminoReport(self):
        """ Coordinate and consensus are KFGPR, but seed is KFPR.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R5-seed,15,0,9,0,AAATTTGGCCCCCGA
""")

        # seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,
        #         A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip
        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R5-seed,R5,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R5-seed,R5,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R5-seed,R5,15,7,3,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R5-seed,R5,15,10,4,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R5-seed,R5,15,13,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,9
"""

        self.report.write_amino_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_amino_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testDeletionWithMinorityVariant(self):
        """ Aligned reads are mostly K-R, but some are KFR.

        Must be a deletion in the sample with respect to the seed reference,
        but some variants in the sample do not have that deletion.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,5,0,AAA---CGA
R1-seed,15,0,2,0,AAATTTCGA
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7
R1-seed,R1,15,4,2,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,0,0,0,7
R1-seed,R1,15,7,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,0,0,0,0,0,0,0,0,0,0,0,0,7
"""

        self.report.read(aligned_reads)
        self.report.write_amino_header(self.report_file)
        self.report.write_amino_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testDeletionNotAlignedToCodons(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,5,0,AAAC---GA
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5
R1-seed,R1,15,4,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,0,0,0,5
R1-seed,R1,15,7,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,5
"""
        self.report.remap_conseqs = {'R1-seed': 'AAATTTAGG'}

        self.report.read(aligned_reads)
        self.report.write_amino_header(self.report_file)
        self.report.write_amino_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testInsertionBetweenSeedAndCoordinateAminoReport(self):
        """ Coordinate sequence is KFQTPREH, and this aligned read is HERKFQTGPREHQFK.

        The G must be an insertion in the seed reference with respect to the
        coordinate reference.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R3-seed,15,0,9,0,CATGAGCGAAAATTTCAGACTGGGCCCCGAGAGCATCAGTTTAAA
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R3-seed,R3,15,10,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R3-seed,R3,15,13,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R3-seed,R3,15,16,3,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R3-seed,R3,15,19,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,9,0,0,9
R3-seed,R3,15,25,5,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R3-seed,R3,15,28,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,9
R3-seed,R3,15,31,7,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R3-seed,R3,15,34,8,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
"""
        expected_insertions = """\
seed,region,qcut,left,insert,count,before
R3-seed,R3,15,22,G,9,5
"""

        self.report.read(aligned_reads)
        self.report.write_insertions()
        self.report.write_nuc_header(StringIO())
        self.report.write_nuc_counts()  # calculates insertion counts
        self.report.write_amino_header(self.report_file)
        self.report.write_amino_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())
        self.assertMultiLineEqual(expected_insertions,
                                  self.insertion_file.getvalue())

    def testInsertionBetweenSeedAndCoordinateNucleotideReport(self):
        """ Coordinate sequence is KFQTPREH, and this aligned read is HERKFQTGPREHQFK.

        The G must be an insertion in the seed reference with respect to the
        coordinate reference.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R3-seed,15,0,9,0,CATGAGCGAAAATTTCAGACTGGGCCCCGAGAGCATCAGTTTAAA
""")
        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T,N,del,ins,clip,v3_overlap,coverage
R3-seed,R3,15,10,1,9,0,0,0,0,0,0,0,0,9
R3-seed,R3,15,11,2,9,0,0,0,0,0,0,0,0,9
R3-seed,R3,15,12,3,9,0,0,0,0,0,0,0,0,9
R3-seed,R3,15,13,4,0,0,0,9,0,0,0,0,0,9
R3-seed,R3,15,14,5,0,0,0,9,0,0,0,0,0,9
R3-seed,R3,15,15,6,0,0,0,9,0,0,0,0,0,9
R3-seed,R3,15,16,7,0,9,0,0,0,0,0,0,0,9
R3-seed,R3,15,17,8,9,0,0,0,0,0,0,0,0,9
R3-seed,R3,15,18,9,0,0,9,0,0,0,0,0,0,9
R3-seed,R3,15,19,10,9,0,0,0,0,0,0,0,0,9
R3-seed,R3,15,20,11,0,9,0,0,0,0,0,0,0,9
R3-seed,R3,15,21,12,0,0,0,9,0,0,9,0,0,9
R3-seed,R3,15,25,13,0,9,0,0,0,0,0,0,0,9
R3-seed,R3,15,26,14,0,9,0,0,0,0,0,0,0,9
R3-seed,R3,15,27,15,0,9,0,0,0,0,0,0,0,9
R3-seed,R3,15,28,16,0,9,0,0,0,0,0,0,0,9
R3-seed,R3,15,29,17,0,0,9,0,0,0,0,0,0,9
R3-seed,R3,15,30,18,9,0,0,0,0,0,0,0,0,9
R3-seed,R3,15,31,19,0,0,9,0,0,0,0,0,0,9
R3-seed,R3,15,32,20,9,0,0,0,0,0,0,0,0,9
R3-seed,R3,15,33,21,0,0,9,0,0,0,0,0,0,9
R3-seed,R3,15,34,22,0,9,0,0,0,0,0,0,0,9
R3-seed,R3,15,35,23,9,0,0,0,0,0,0,0,0,9
R3-seed,R3,15,36,24,0,0,0,9,0,0,0,0,0,9
"""

        self.report.read(aligned_reads)
        self.report.write_insertions()
        self.report.write_nuc_header(self.report_file)
        self.report.write_nuc_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testInsertionsSortedByCount(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R3-seed,15,0,9,0,CATGAGCGAAAATTTCAGACTGGGCCCCGAGAGCATCAGTTTAAA
R3-seed,15,0,8,0,CATGAGCGAAAATTTCAGACTAAACCCCGAGAGCATCAGTTTAAA
""")
        expected_insertions = """\
seed,region,qcut,left,insert,count,before
R3-seed,R3,15,22,G,9,5
R3-seed,R3,15,22,K,8,5
"""

        self.report.read(aligned_reads)
        self.report.write_insertions()

        self.assertMultiLineEqual(expected_insertions,
                                  self.insertion_file.getvalue())

    def testInsertionsSortedByLeft(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R3-seed,15,0,9,0,CATGAGCGAAAATTTCAGACTGGGCCCCGAAAAGAGCATCAGTTTAAA
""")
        expected_insertions = """\
seed,region,qcut,left,insert,count,before
R3-seed,R3,15,22,G,9,5
R3-seed,R3,15,31,K,9,7
"""

        self.report.read(aligned_reads)
        self.report.write_insertions()

        self.assertMultiLineEqual(expected_insertions,
                                  self.insertion_file.getvalue())

    def testInsertionInDifferentReadingFrame(self):
        """ Delete part of the first codon to throw off the reading frame.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R3-seed,15,0,9,0,AATTTCAGACTGGGCCCCGAGAGCAT
""")

        expected_insertions = """\
seed,region,qcut,left,insert,count,before
R3-seed,R3,15,12,G,9,5
"""

        self.report.read(aligned_reads)
        self.report.write_amino_header(StringIO())
        self.report.write_amino_counts()
        self.report.write_insertions()

        self.assertMultiLineEqual(expected_insertions,
                                  self.insertion_file.getvalue())

    def testInsertionInSomeReads(self):
        """ Not all reads have the insertion, some end before it.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R3-seed,15,0,9,0,AAATTTCAGACTGGGCCCCGAGAGCAT
R3-seed,15,1,5,0,AAATTTCAG
R3-seed,15,2,4,0,AAATTTCAGACTG
""")

        expected_insertions = """\
seed,region,qcut,left,insert,count,before
R3-seed,R3,15,13,G,9,5
"""

        self.report.read(aligned_reads)
        self.report.write_amino_header(StringIO())
        self.report.write_amino_counts()
        self.report.write_insertions()

        self.assertMultiLineEqual(expected_insertions,
                                  self.insertion_file.getvalue())

    def testMultipleCoordinateInsertionReport(self):
        """ Two coordinate regions map the same seed region, the consensus
        has an insertion relative to only one of them.
        """
        self.report.projects.load(StringIO("""\
{
  "projects": {
    "R3": {
      "max_variants": 0,
      "regions": [
        {
          "coordinate_region": "R3a",
          "seed_region_names": ["R3-seed"]
        },
        {
          "coordinate_region": "R3b",
          "seed_region_names": ["R3-seed"]
        }
      ]
    }
  },
  "regions": {
    "R3-seed": {
      "is_nucleotide": true,
      "reference": [
        "AAATTTCAGACCGGGCCACGAGAGCAT"
      ]
    },
    "R3a": {
      "is_nucleotide": false,
      "reference": [
        "KFQTPREH"
      ]
    },
    "R3b": {
      "is_nucleotide": false,
      "reference": [
        "KFQTGPREH"
      ]
    }
  }
}
"""))
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R3-seed,15,0,9,0,AAATTTCAGACTGGGCCCCGAGAGCAT
""")

        expected_insertions = """\
seed,region,qcut,left,insert,count,before
R3-seed,R3a,15,13,G,9,5
"""

        self.report.read(aligned_reads)
        self.report.write_insertions()

        self.assertMultiLineEqual(expected_insertions,
                                  self.insertion_file.getvalue())

    def testGapBetweenForwardAndReverse(self):
        """ Lower-case n represents a gap between forward and reverse reads.

        Region R2 has sequence KFGPR, so this read has a gap at the end of G
        and beginning of P. Partial codons at the ends of a read or next to the
        gap are ignored, even though G is still unambiguous.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R2-seed,15,0,5,0,AAATTTGGnnCCCGA
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R2-seed,R2,15,1,1,0,0,0,0,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5
R2-seed,R2,15,4,2,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5
R2-seed,R2,15,7,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,10,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,13,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,5
"""

        self.report.read(aligned_reads)
        self.report.write_amino_header(self.report_file)
        self.report.write_amino_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testFailedAlignmentAminoReport(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,2,0,TTATCCTAC
""")

        expected_text = ""

        self.report.read(aligned_reads)
        self.report.write_amino_counts(self.report_file)

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testFailedAlignmentFailureReport(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,2,0,TTATCCTAC
""")

        expected_text = """\
seed,region,qcut,queryseq,refseq
R1-seed,R1,15,LSY,KFR
"""

        self.report.write_failure_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_failure()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testFailedAlignmentWithHeadToTailMatch(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R3-seed,15,0,2,0,TTATCCTACTTATCCTACTTATCCAAA
""")

        expected_text = """\
seed,region,qcut,queryseq,refseq
R3-seed,R3,15,LSYLSYLSK,KFQTPREH
"""

        self.report.read(aligned_reads)
        self.report.write_failure_header(self.report_file)
        self.report.write_failure()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testGoodAlignmentWithTinyCoordinateReference(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,2,0,AAATTTCGATTATCCTACTTATCCTACTTATCCTACTTATCCTACTTATCCTACTTATCCTACTTATCCTACTTATCCTAC
""")

        expected_text = """\
seed,region,qcut,queryseq,refseq
"""

        self.report.read(aligned_reads)
        self.report.write_failure_header(self.report_file)
        self.report.write_failure()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testGoodAlignmentWithGiantSeed(self):
        """ Short consensus with long seed and long coordinate reference.

        Even when the consensus maps to the end of the seed, it should still
        only require a low alignment score.
        """
        self.report.projects.load(StringIO("""\
{
  "projects": {
    "R3": {
      "max_variants": 0,
      "regions": [
        {
          "coordinate_region": "R3",
          "seed_region_names": ["R3-seed"]
        }
      ]
    }
  },
  "regions": {
    "R3-seed": {
      "is_nucleotide": true,
      "reference": [
        "TGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGG",
        "TGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGAAATTTAGG"
      ]
    },
    "R3": {
      "is_nucleotide": false,
      "reference": [
        "WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWKFR"
      ]
    }
  }
}
"""))
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R3-seed,15,0,9,123,AAATTTCGA
""")

        # seed,region,qcut,queryseq,refseq
        expected_text = """\
seed,region,qcut,queryseq,refseq
"""

        self.report.read(aligned_reads)
        self.report.write_failure_header(self.report_file)
        self.report.write_failure()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testFailedAlignmentWithOffset(self):
        """ Be careful that an offset from the seed reference doesn't match
        the dashes in the failed alignment.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,2,3,TTATCCTAC
""")

        expected_text = """\
seed,region,qcut,queryseq,refseq
R1-seed,R1,15,-LSY,KFR
"""

        self.report.read(aligned_reads)
        self.report.write_failure_header(self.report_file)
        self.report.write_failure()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testMultipleCoordinateRefsNoAlignment(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R7-seed,15,0,2,0,TTATCCTAC
""")

        expected_text = """\
seed,region,qcut,queryseq,refseq
R7-seed,R7a,15,LSY,KFQ
R7-seed,R7b,15,LSY,REH
"""

        self.report.read(aligned_reads)
        self.report.write_failure_header(self.report_file)
        self.report.write_failure()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testMultipleCoordinateRefsOneAlignment(self):
        """ If one coordinate aligns, don't complain about the others.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R7-seed,15,0,2,0,AAATTT
""")

        expected_text = """\
seed,region,qcut,queryseq,refseq
"""

        self.report.read(aligned_reads)
        self.report.write_failure_header(self.report_file)
        self.report.write_failure()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testNoFailureReport(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,9,0,AAATTT
""")

        expected_text = """\
seed,region,qcut,queryseq,refseq
"""

        self.report.read(aligned_reads)
        self.report.write_failure_header(self.report_file)
        self.report.write_failure()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testMultipleCoordinateAminoReport(self):
        """ Two coordinate regions map the same seed region, report both.
        """
        self.report.projects.load(StringIO("""\
{
  "projects": {
    "R1": {
      "max_variants": 10,
      "regions": [
        {
          "coordinate_region": "R1a",
          "seed_region_names": ["R1-seed"]
        },
        {
          "coordinate_region": "R1b",
          "seed_region_names": ["R1-seed"]
        }
      ]
    }
  },
  "regions": {
    "R1-seed": {
      "is_nucleotide": true,
      "reference": [
        "TGGAAATTTAGG"
      ]
    },
    "R1a": {
      "is_nucleotide": false,
      "reference": [
        "KFR"
      ]
    },
    "R1b": {
      "is_nucleotide": false,
      "reference": [
        "WKFR"
      ]
    }
  }
}
"""))
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,9,0,AAATTT
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R1-seed,R1a,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1a,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1a,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R1-seed,R1b,15,,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R1-seed,R1b,15,1,2,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1b,15,4,3,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1b,15,,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
"""

        self.report.read(aligned_reads)
        self.report.write_amino_header(self.report_file)
        self.report.write_amino_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testMultipleCoordinateConsensusRegionsReport(self):
        """ Two coordinate regions map the same seed region, report both.
        """
        self.report.projects.load(StringIO("""\
{
  "projects": {
    "R1": {
      "max_variants": 10,
      "regions": [
        {
          "coordinate_region": "R1a",
          "seed_region_names": ["R1-seed"]
        },
        {
          "coordinate_region": "R1b",
          "seed_region_names": ["R1-seed"]
        }
      ]
    }
  },
  "regions": {
    "R1-seed": {
      "is_nucleotide": true,
      "reference": [
        "TGGAAATTTAGG"
      ]
    },
    "R1a": {
      "is_nucleotide": false,
      "reference": [
        "KFR"
      ]
    },
    "R1b": {
      "is_nucleotide": false,
      "reference": [
        "WKFR"
      ]
    }
  }
}
"""))
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,9,0,AAATTT
R1-seed,15,0,1,0,CCCGGG
""")

        expected_text = """\
seed,region,q-cutoff,min_coverage,consensus-percent-cutoff,offset,sequence
R1-seed,R1a,15,0,MAX,0,AAATTT
R1-seed,R1a,15,0,0.100,0,MMMKKK
R1-seed,R1b,15,0,MAX,3,AAATTT
R1-seed,R1b,15,0,0.100,3,MMMKKK
"""

        self.report.read(aligned_reads)
        self.report.write_consensus_regions_header(self.report_file)
        self.report.write_consensus_regions()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testReadRemapConseqs(self):
        remap_conseqs_csv = StringIO("""\
region,sequence
R1,ACATAGCCCGGG
R2,GCCATTAAA
""")
        expected_conseqs = {'R1': 'ACATAGCCCGGG', 'R2': 'GCCATTAAA'}

        self.report.read_remap_conseqs(remap_conseqs_csv)

        self.assertEqual(expected_conseqs, self.report.remap_conseqs)

    def testAlignDeletionsWithoutDeletion(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("R1-seed,15,0,10,0,AAATTTAGG")
        self.report.remap_conseqs = {'R1-seed': 'AAATTTAGG'}
        expected_reads = [dict(refname='R1-seed',
                               qcut='15',
                               rank='0',
                               count='10',
                               offset='0',
                               seq='AAATTTAGG')]

        reads = list(self.report.align_deletions(aligned_reads))

        self.assertEqual(expected_reads, reads)

    def testAlignDeletionsNoChange(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("R1-seed,15,0,10,0,AAA---AGG")
        self.report.remap_conseqs = {'R1-seed': 'AAATTTAGG'}
        expected_reads = [dict(refname='R1-seed',
                               qcut='15',
                               rank='0',
                               count='10',
                               offset='0',
                               seq='AAA---AGG')]

        reads = list(self.report.align_deletions(aligned_reads))

        self.assertEqual(expected_reads, reads)

    def testAlignDeletionsShiftedRight(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("R1-seed,15,0,10,0,AAAA---GG")
        self.report.remap_conseqs = {'R1-seed': 'AAATTTAGG'}
        expected_reads = [dict(refname='R1-seed',
                               qcut='15',
                               rank='0',
                               count='10',
                               offset='0',
                               seq='AAA---AGG')]

        reads = list(self.report.align_deletions(aligned_reads))

        self.assertEqual(expected_reads, reads)

    def testAlignDeletionsShiftedLeft(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("R1-seed,15,0,10,0,AA---AAGG")
        self.report.remap_conseqs = {'R1-seed': 'AAATTTAGG'}
        expected_reads = [dict(refname='R1-seed',
                               qcut='15',
                               rank='0',
                               count='10',
                               offset='0',
                               seq='AAA---AGG')]

        reads = list(self.report.align_deletions(aligned_reads))

        self.assertEqual(expected_reads, reads)

    def testAlignDeletionsTwoCodons(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("R2-seed,15,0,10,0,AA------ACCGAGA")
        self.report.remap_conseqs = {'R2-seed': 'AAATTTGGCCCGAGA'}
        expected_reads = [dict(refname='R2-seed',
                               qcut='15',
                               rank='0',
                               count='10',
                               offset='0',
                               seq='AAA------CCGAGA')]

        reads = list(self.report.align_deletions(aligned_reads))

        self.assertEqual(expected_reads, reads)

    def testAlignDeletionsUsingOffset(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("R1-seed,15,0,10,1,AA---AGG")
        self.report.remap_conseqs = {'R1-seed': 'AAATTTAGG'}
        expected_reads = [dict(refname='R1-seed',
                               qcut='15',
                               rank='0',
                               count='10',
                               offset='1',
                               seq='AA---AGG')]

        reads = list(self.report.align_deletions(aligned_reads))

        self.assertEqual(expected_reads, reads)

    def testAlignDeletionsUsingReadingFrame1(self):
        self.report.projects.load(StringIO("""\
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
            }
          },
          "regions": {
            "R1-seed": {
              "is_nucleotide": true,
              "reference": ["CCAAATTTAGG"]
            },
            "R1": {
              "is_nucleotide": false,
              "reference": ["KFR"]
            }
          }
        }
        """))
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("R1-seed,15,0,10,2,AAA---AGG")
        self.report.remap_conseqs = {'R1-seed': 'CCAAATTTAGG'}
        expected_reads = [dict(refname='R1-seed',
                               qcut='15',
                               rank='0',
                               count='10',
                               offset='2',
                               seq='AAA---AGG')]

        reads = list(self.report.align_deletions(aligned_reads))

        self.assertEqual(expected_reads, reads)

    def testAlignDeletionsMultipleReadingFrames(self):
        self.report.projects.load(StringIO("""\
{
  "projects": {
    "Rs": {
      "max_variants": 10,
      "regions": [
        {
          "coordinate_region": "R1",
          "seed_region_names": ["R-seed"]
        },
        {
          "coordinate_region": "R2",
          "seed_region_names": ["R-seed"]
        }
      ]
    }
  },
  "regions": {
    "R-seed": {
      "is_nucleotide": true,
      "reference": ["GAAATTTCAGTTTTTTTTCGAGAGCAT"],
      "comment": "   ^KkkFffQqq^^^^^^^^RrrEeeHhh (two reading frames)"
    },
    "R1": {
      "is_nucleotide": false,
      "reference": ["KFQ"]
    },
    "R2": {
      "is_nucleotide": false,
      "reference": ["REH"]
    }
  }
}
"""))
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("R-seed,15,0,10,1,AAA---CAGTTTTTTTTC---AGCAT")
        self.report.remap_conseqs = {'R-seed': 'GAAATTTCAGTTTTTTTTCGAGAGCAT'}
        expected_reads = [dict(refname='R-seed',
                               qcut='15',
                               rank='0',
                               count='10',
                               offset='1',
                               seq='AAA---CAGTTTTTTTT---CAGCAT')]

        reads = list(self.report.align_deletions(aligned_reads))

        self.assertEqual(expected_reads, reads)

    def testCombineDeletions(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("R2-seed,15,0,10,0,AA-TCG--CCCGAGA")
        self.report.remap_conseqs = {'R2-seed': 'AAATTTGGCCCGAGA'}
        expected_reads = [dict(refname='R2-seed',
                               qcut='15',
                               rank='0',
                               count='10',
                               offset='0',
                               seq='AATCGC---CCGAGA')]

        reads = list(self.report.align_deletions(aligned_reads))

        self.assertEqual(expected_reads, reads)

    def testCombineDeletionsTwoCodons(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("R2-seed,15,0,10,0,AA--TC----CGAGA")
        self.report.remap_conseqs = {'R2-seed': 'AAATTTGGCCCGAGA'}
        expected_reads = [dict(refname='R2-seed',
                               qcut='15',
                               rank='0',
                               count='10',
                               offset='0',
                               seq='AAT------CCGAGA')]

        reads = list(self.report.align_deletions(aligned_reads))

        self.assertEqual(expected_reads, reads)

    def testCombineDeletionsMaxSpread(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("R3-seed,15,0,10,0,AA--TTCAGACCCC-CGAGAGCAT")
        self.report.remap_conseqs = {'R3-seed': 'AAATTTCAGACCCCACGAGAGCAT'}
        expected_reads = [dict(refname='R3-seed',
                               qcut='15',
                               rank='0',
                               count='10',
                               offset='0',
                               seq='AAT---TCAGACCCCCGAGAGCAT')]

        reads = list(self.report.align_deletions(aligned_reads))

        self.assertEqual(expected_reads, reads)

    def testCombineDeletionsBeyondMaxSpread(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("R3-seed,15,0,10,0,AA--TTCAGACCCCA-GAGAGCAT")
        self.report.remap_conseqs = {'R3-seed': 'AAATTTCAGACCCCACGAGAGCAT'}
        expected_reads = [dict(refname='R3-seed',
                               qcut='15',
                               rank='0',
                               count='10',
                               offset='0',
                               seq='AA--TTCAGACCCCA-GAGAGCAT')]

        reads = list(self.report.align_deletions(aligned_reads))

        self.assertEqual(expected_reads, reads)

    def testCombineDeletionsTooCrowded(self):
        """ There must be a buffer of 13 with no deletions.

        Otherwise, sweeping is not allowed.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("R3-seed,15,0,10,0,AA--TTCAGACCCC-CGA---CAT")
        self.report.remap_conseqs = {'R3-seed': 'AAATTTCAGACCCCACGAGAGCAT'}
        expected_reads = [dict(refname='R3-seed',
                               qcut='15',
                               rank='0',
                               count='10',
                               offset='0',
                               seq='AA--TTCAGACCCC-CGA---CAT')]

        reads = list(self.report.align_deletions(aligned_reads))

        self.assertEqual(expected_reads, reads)


class InsertionWriterTest(unittest.TestCase):
    def setUp(self):
        self.insert_file = StringIO()
        self.writer = InsertionWriter(self.insert_file)
        self.writer.start_group(seed='R1-seed', qcut=15)
        self.nuc_seq_acdef = 'GCTTGTGACGAGTTT'
        self.nuc_seq_afdef = 'GCTTTTGACGAGTTT'

    def testNoInserts(self):
        expected_text = """\
seed,region,qcut,left,insert,count,before
"""

        self.writer.add_nuc_read(offset_sequence=self.nuc_seq_acdef, count=1)
        self.writer.write(inserts=[], region='R1')

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())

    def testInsert(self):
        expected_text = """\
seed,region,qcut,left,insert,count,before
R1-seed,R1,15,7,D,1,
"""

        self.writer.add_nuc_read(offset_sequence=self.nuc_seq_acdef, count=1)
        self.writer.write(inserts=[6], region='R1')

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())

    def testInsertWithBefore(self):
        expected_text = """\
seed,region,qcut,left,insert,count,before
R1-seed,R1,15,10,E,1,2
"""
        expected_counts = {('R1-seed', 'R1'): {1: 1}}
        seed_amino_after = SeedAmino(12)
        report_amino_after = ReportAmino(seed_amino_after, position=2)

        self.writer.add_nuc_read(offset_sequence=self.nuc_seq_acdef, count=1)
        self.writer.write(inserts=[9],
                          region='R1',
                          report_aminos=[report_amino_after])

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())
        self.assertEqual(expected_counts, self.writer.insert_pos_counts)

    def testInsertDifferentReadingFrame(self):
        """ Add a partial codon at the start of the read to shift the reading
        frame.
        """
        expected_text = """\
seed,region,qcut,left,insert,count,before
R1-seed,R1,15,8,D,1,
"""

        self.writer.add_nuc_read(offset_sequence='A' + self.nuc_seq_acdef,
                                 count=1)
        self.writer.write(inserts=[7], region='R1')

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())

    def testInsertWithOffset(self):
        expected_text = """\
seed,region,qcut,left,insert,count,before
R1-seed,R1,15,7,D,1,
"""

        #                                            C  D  E  F
        self.writer.add_nuc_read(offset_sequence='---TGTGACGAGTTT', count=1)
        self.writer.write(inserts=[6], region='R1')

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())

    def testInsertWithDeletion(self):
        expected_text = """\
seed,region,qcut,left,insert,count,before
"""

        #                                         C  D     E  F
        self.writer.add_nuc_read(offset_sequence='TGTGAC---GAGTTT', count=1)
        self.writer.write(inserts=[3, 6], region='R1')

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())

    def testTwoInsertsWithOffset(self):
        expected_text = """\
seed,region,qcut,left,insert,count,before
R1-seed,R1,15,7,D,1,
R1-seed,R1,15,13,F,1,
"""

        #                                            C  D  E  F  G
        self.writer.add_nuc_read(offset_sequence='---TGTGACGAGTTTGGG', count=1)
        self.writer.write(inserts=[6, 12], region='R1')

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())

    def testInsertsWithVariants(self):
        expected_text = """\
seed,region,qcut,left,insert,count,before
R1-seed,R1,15,7,D,2,
"""

        self.writer.add_nuc_read(offset_sequence=self.nuc_seq_acdef, count=1)
        self.writer.add_nuc_read(offset_sequence=self.nuc_seq_afdef, count=1)
        self.writer.write(inserts=[6], region='R1')

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())

    def testDifferentInserts(self):
        expected_text = """\
seed,region,qcut,left,insert,count,before
R1-seed,R1,15,4,C,2,
R1-seed,R1,15,4,F,3,"""

        self.writer.add_nuc_read(offset_sequence=self.nuc_seq_acdef, count=2)
        self.writer.add_nuc_read(offset_sequence=self.nuc_seq_afdef, count=3)
        self.writer.write(inserts=[3], region='R1')

        lines = self.insert_file.getvalue().splitlines()
        lines.sort()
        lines.insert(0, lines.pop())
        text = '\n'.join(lines)
        self.assertMultiLineEqual(expected_text, text)

    def testMulticharacterInsert(self):
        expected_text = """\
seed,region,qcut,left,insert,count,before
R1-seed,R1,15,7,DE,1,
"""

        self.writer.add_nuc_read(offset_sequence=self.nuc_seq_acdef, count=1)
        self.writer.write(inserts=[6, 9], region='R1')

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())

    def testReadGapInInsert(self):
        nuc_seq = 'GCTCTnGACGAGTTT'

        expected_text = """\
seed,region,qcut,left,insert,count,before
"""

        self.writer.add_nuc_read(nuc_seq, count=1)
        self.writer.write(inserts=[3], region='R1')

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())

    def testUnsortedInserts(self):
        expected_text = """\
seed,region,qcut,left,insert,count,before
R1-seed,R1,15,7,DE,1,
"""

        self.writer.add_nuc_read(offset_sequence=self.nuc_seq_acdef, count=1)
        self.writer.write(inserts=(9, 6), region='R1')

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())
