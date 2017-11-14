import csv
from io import StringIO
import sys
import unittest

from micall.core.aln2counts import SequenceReport, SeedNucleotide,\
    InsertionWriter, MAX_CUTOFF, SeedAmino, ReportAmino
from micall.core import project_config


class StubbedSequenceReport(SequenceReport):
    def __init__(self,
                 insert_writer,
                 projects,
                 conseq_mixture_cutoffs):
        SequenceReport.__init__(self,
                                insert_writer,
                                projects,
                                conseq_mixture_cutoffs)
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
        conseq_mixture_cutoffs = [0.1]
        self.report = StubbedSequenceReport(insert_writer,
                                            projects,
                                            conseq_mixture_cutoffs)
        self.report_file = StringIO()

    @staticmethod
    def prepareReads(aligned_reads_text):
        full_text = "refname,qcut,rank,count,offset,seq\n" + aligned_reads_text
        dummy_file = StringIO(full_text)
        return csv.DictReader(dummy_file)

    def testEmptyAminoReport(self):
        expected_text = ""
        aligned_reads = []

        self.report.read(aligned_reads)
        self.report.write_amino_counts(self.report_file)

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
region,q-cutoff,consensus-percent-cutoff,offset,sequence
R1-seed,15,MAX,0,AAATTT
R1-seed,15,0.100,0,AAATTT
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
region,q-cutoff,consensus-percent-cutoff,offset,sequence
R1-seed,15,MAX,0,AAATTT
R1-seed,15,0.100,0,MMMKKK
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
region,q-cutoff,consensus-percent-cutoff,offset,sequence
R1-seed,15,MAX,3,AAATTTGGG
R1-seed,15,0.100,3,AAATTTGGG
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
region,q-cutoff,consensus-percent-cutoff,offset,sequence
R1-seed,15,MAX,6,TTTNGG
R1-seed,15,0.100,6,TTTNGG
"""

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
region,q-cutoff,consensus-percent-cutoff,offset,sequence
"""

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
region,q-cutoff,consensus-percent-cutoff,offset,sequence
R1-seed,15,MAX,0,AAAT--GGG
R1-seed,15,0.100,0,AAAT--GGG
"""
        self.report.consensus_min_coverage=10

        self.report.write_consensus_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_consensus()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testConsensusLowCoverageAtStart(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,9,0,AAATTTGGG
R1-seed,15,0,1,3,TTTGGG
""")
        expected_text = """\
region,q-cutoff,consensus-percent-cutoff,offset,sequence
R1-seed,15,MAX,3,TTTGGG
R1-seed,15,0.100,3,TTTGGG
"""
        self.report.consensus_min_coverage=10

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

    def testSecondSourceAminoReport(self):
        """ In this sample, there are two sequences, each with two codons.
        AAA -> K
        TTT -> F
        The coordinate references have three codons and five codons, so the
        later positions are empty.
        """
        g2p_aligned_csv = StringIO("""\
refname,qcut,rank,count,offset,seq
R1-seed,15,0,9,0,AAATTT
""")
        aligned_csv = StringIO("""\
refname,qcut,rank,count,offset,seq
R2-seed,15,0,8,0,AAATTT
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R2-seed,R2,15,1,1,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8
R2-seed,R2,15,4,2,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8
R2-seed,R2,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
"""

        self.report.write_amino_header(self.report_file)
        self.report.process_reads(g2p_aligned_csv, aligned_csv)

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testAminoReportWithG2pOverlap(self):
        """ If the same region appears in both sources, the second is overlap.
        """
        g2p_aligned_csv = StringIO("""\
refname,qcut,rank,count,offset,seq
R6a-seed,15,0,9,0,AAATTT
""")
        aligned_csv = StringIO("""\
refname,qcut,rank,count,offset,seq
R6b-seed,15,0,8,3,AAATTT
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R6a-seed,R6,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,9
R6a-seed,R6,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,9
R6a-seed,R6,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
"""

        self.report.write_amino_header(self.report_file)
        self.report.process_reads(g2p_aligned_csv, aligned_csv, g2p_region_name='R6')

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testAminoReportWithoutG2pOverlap(self):
        """ Bowtie2 aligned reads don't overlap the coordinate reference at all.
        """
        g2p_aligned_csv = StringIO("""\
refname,qcut,rank,count,offset,seq
R6a-seed,15,0,9,0,AAATTT
""")
        aligned_csv = StringIO("""\
refname,qcut,rank,count,offset,seq
R6b-seed,15,0,8,15,GGGGGGGGG
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R6a-seed,R6,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R6a-seed,R6,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R6a-seed,R6,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
"""

        self.report.write_amino_header(self.report_file)
        self.report.process_reads(g2p_aligned_csv, aligned_csv, g2p_region_name='R6')

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testAminoReportWithIndirectG2pOverlap(self):
        """ Bowtie2 aligned reads map to a different part of the coordinate ref.
        """
        g2p_aligned_csv = StringIO("""\
refname,qcut,rank,count,offset,seq
R6a-seed,15,0,9,0,AAATTT
""")
        aligned_csv = StringIO("""\
refname,qcut,rank,count,offset,seq
R6b-seed,15,0,8,9,AGG
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R6a-seed,R6,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R6a-seed,R6,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R6a-seed,R6,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0
"""

        self.report.write_amino_header(self.report_file)
        self.report.process_reads(g2p_aligned_csv, aligned_csv, g2p_region_name='R6')

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testInsertionReportWithG2pOverlap(self):
        """ If the same region appears in both sources, the second is overlap.
        """
        g2p_aligned_csv = StringIO("""\
refname,qcut,rank,count,offset,seq
R7-seed,15,0,9,0,AAATTTCAGACCCCACGAGAGCAT
""")
        aligned_csv = StringIO("""\
refname,qcut,rank,count,offset,seq
R7-seed,15,0,8,6,AAATTTCAGACCCCACGAGAGCAT
""")

        expected_text = """\
seed,region,qcut,left,insert,count,before
"""

        self.report.process_reads(g2p_aligned_csv, aligned_csv, g2p_region_name='R7a')

        self.assertMultiLineEqual(expected_text, self.insertion_file.getvalue())

    def testAminoReportWithDifferingConsensus(self):
        """ Some reads that didn't map in G2P, mapped to same seed with bowtie2.
        """
        self.report.projects.load(StringIO("""\
{
  "projects": {
    "HIV": {
      "max_variants": 0,
      "regions": [
        {
          "coordinate_region": "R1",
          "seed_region_names": ["HIV1-CON-XX-Consensus-seed"]
        }
      ]
    }
  },
  "regions": {
    "HIV1-CON-XX-Consensus-seed": {
      "is_nucleotide": true,
      "reference": [
        "GAAATTTGGCCCGAGA"
      ]
    },
    "R1": {
      "is_nucleotide": false,
      "reference": [
        "KFGPR"
      ]
    }
  }
}
"""))
        self.report.remap_conseqs = {'HIV1-CON-XX-Consensus-seed': "GAAATTTCAAGGCCCGAGA"}
        # Insert a Q in the middle                                      Q^^
        g2p_aligned_csv = StringIO("""\
refname,qcut,rank,count,offset,seq
HIV1-CON-XX-Consensus-seed,15,0,9,1,AAATTTGGC
""")
        aligned_csv = StringIO("""\
refname,qcut,rank,count,offset,seq
HIV1-CON-XX-Consensus-seed,15,0,6,1,AAATTTCAAGGCCCG
HIV1-CON-XX-Consensus-seed,15,0,2,1,AA---TCAAGGCCCG
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
HIV1-CON-XX-Consensus-seed,R1,15,2,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,9
HIV1-CON-XX-Consensus-seed,R1,15,5,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,9
HIV1-CON-XX-Consensus-seed,R1,15,8,3,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,9
HIV1-CON-XX-Consensus-seed,R1,15,,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0
HIV1-CON-XX-Consensus-seed,R1,15,,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
"""

        self.report.write_amino_header(self.report_file)
        self.report.process_reads(g2p_aligned_csv, aligned_csv, g2p_region_name='R1')

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testAminoReportOverlapWithDifferentSeed(self):
        """ Some reads that didn't map in G2P, mapped to different seed with bowtie2.
        """
        self.report.projects.load(StringIO("""\
{
  "projects": {
    "HIV": {
      "max_variants": 0,
      "regions": [
        {
          "coordinate_region": "R1",
          "seed_region_names": ["HIV1-CON-XX-Consensus-seed", "R1-seed"]
        }
      ]
    }
  },
  "regions": {
    "HIV1-CON-XX-Consensus-seed": {
      "is_nucleotide": true,
      "reference": [
        "GAAATTTGGCCCGAGA"
      ]
    },
    "R1-seed": {
      "is_nucleotide": true,
      "reference": [
        "AAATTTGGCCCGAGA"
      ]
    },
    "R1": {
      "is_nucleotide": false,
      "reference": [
        "KFGPR"
      ]
    }
  }
}
"""))
        self.report.remap_conseqs = {'R1-seed': "AAATTTGGCCCGAGA"}
        g2p_aligned_csv = StringIO("""\
refname,qcut,rank,count,offset,seq
HIV1-CON-XX-Consensus-seed,15,0,9,0,GAAATTTGGC
""")
        aligned_csv = StringIO("""\
refname,qcut,rank,count,offset,seq
R1-seed,15,0,8,0,AAATTTGGCCCG
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
HIV1-CON-XX-Consensus-seed,R1,15,2,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,9
HIV1-CON-XX-Consensus-seed,R1,15,5,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,9
HIV1-CON-XX-Consensus-seed,R1,15,8,3,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,9
HIV1-CON-XX-Consensus-seed,R1,15,,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0
HIV1-CON-XX-Consensus-seed,R1,15,,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
"""

        self.report.write_amino_header(self.report_file)
        self.report.process_reads(g2p_aligned_csv, aligned_csv, g2p_region_name='R1')

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testSecondSourceNucleotideReport(self):
        """ In this sample, there are two sequences, each with two codons.

        The coordinate references have three codons and five codons, so the
        later positions are empty.
        """
        g2p_aligned_csv = StringIO("""\
refname,qcut,rank,count,offset,seq
R1-seed,15,0,9,0,AAATTT
""")
        aligned_csv = StringIO("""\
refname,qcut,rank,count,offset,seq
R2-seed,15,0,8,0,AAATTT
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T,N,del,ins,clip,v3_overlap,coverage
R2-seed,R2,15,1,1,8,0,0,0,0,0,0,0,0,8
R2-seed,R2,15,2,2,8,0,0,0,0,0,0,0,0,8
R2-seed,R2,15,3,3,8,0,0,0,0,0,0,0,0,8
R2-seed,R2,15,4,4,0,0,0,8,0,0,0,0,0,8
R2-seed,R2,15,5,5,0,0,0,8,0,0,0,0,0,8
R2-seed,R2,15,6,6,0,0,0,8,0,0,0,0,0,8
R2-seed,R2,15,,7,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,,8,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,,9,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,,10,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,,11,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,,12,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,,13,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,,14,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,,15,0,0,0,0,0,0,0,0,0,0
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
        self.report.process_reads(g2p_aligned_csv, aligned_csv)

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testNucleotideReportWithG2pOverlap(self):
        g2p_aligned_csv = StringIO("""\
refname,qcut,rank,count,offset,seq
R1-seed,15,0,9,0,AAATTT
""")
        aligned_csv = StringIO("""\
refname,qcut,rank,count,offset,seq
R1-seed,15,0,8,0,AAATTT
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T,N,del,ins,clip,v3_overlap,coverage
R1-seed,R1,15,1,1,9,0,0,0,0,0,0,0,8,9
R1-seed,R1,15,2,2,9,0,0,0,0,0,0,0,8,9
R1-seed,R1,15,3,3,9,0,0,0,0,0,0,0,8,9
R1-seed,R1,15,4,4,0,0,0,9,0,0,0,0,8,9
R1-seed,R1,15,5,5,0,0,0,9,0,0,0,0,8,9
R1-seed,R1,15,6,6,0,0,0,9,0,0,0,0,8,9
R1-seed,R1,15,,7,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,,8,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,,9,0,0,0,0,0,0,0,0,0,0
"""

        self.report.write_nuc_header(self.report_file)
        self.report.process_reads(g2p_aligned_csv, aligned_csv, g2p_region_name='R1')

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

    def testRegionWithoutCoordinateReferenceNucleotideReport(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R-NO-COORD,15,0,9,0,AAATTT
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T,N,del,ins,clip,v3_overlap,coverage
R-NO-COORD,R-NO-COORD,15,1,,9,0,0,0,0,0,0,0,0,9
R-NO-COORD,R-NO-COORD,15,2,,9,0,0,0,0,0,0,0,0,9
R-NO-COORD,R-NO-COORD,15,3,,9,0,0,0,0,0,0,0,0,9
R-NO-COORD,R-NO-COORD,15,4,,0,0,0,9,0,0,0,0,0,9
R-NO-COORD,R-NO-COORD,15,5,,0,0,0,9,0,0,0,0,0,9
R-NO-COORD,R-NO-COORD,15,6,,0,0,0,9,0,0,0,0,0,9
R-NO-COORD,R-NO-COORD,15,7,,0,0,0,0,0,0,0,0,0,0
R-NO-COORD,R-NO-COORD,15,8,,0,0,0,0,0,0,0,0,0,0
R-NO-COORD,R-NO-COORD,15,9,,0,0,0,0,0,0,0,0,0,0
"""

        self.report.read(aligned_reads)
        self.report.write_nuc_header(self.report_file)
        self.report.write_nuc_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testRegionWithoutCoordinateReferenceFailureReport(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R-NO-COORD,15,0,9,0,AAATTT
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


class SeedAminoTest(unittest.TestCase):
    def setUp(self):
        self.amino = SeedAmino(None)

    def testSingleRead(self):
        """ Read a single codon, and report on counts.
        Columns are:       A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
        """
        nuc_seq = 'AAA'  # -> K
        expected_counts = '0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0'

        self.amino.count_aminos(nuc_seq, 8)
        counts = self.amino.get_report()

        self.assertSequenceEqual(expected_counts, counts)

    def testDifferentCodon(self):
        """ Read two different codons, and report on counts.
        Columns are:       A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
        """
        nuc_seq1 = 'AAA'  # -> K
        nuc_seq2 = 'GGG'  # -> G
        expected_counts = '0,0,0,0,0,5,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0'

        self.amino.count_aminos(nuc_seq1, 8)
        self.amino.count_aminos(nuc_seq2, 5)
        counts = self.amino.get_report()

        self.assertSequenceEqual(expected_counts, counts)

    def testSameAminoAcid(self):
        """ Read same codon twice, and report on counts.
        Columns are:       A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
        """
        nuc_seq1 = 'AAA'  # -> K
        nuc_seq2 = 'AAG'  # -> K
        expected_counts = '0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0'

        self.amino.count_aminos(nuc_seq1, 4)
        self.amino.count_aminos(nuc_seq2, 5)
        counts = self.amino.get_report()

        self.assertSequenceEqual(expected_counts, counts)

    def testNucleotides(self):
        nuc_seq1 = 'AAA'  # -> K
        nuc_seq2 = 'AAG'  # -> K
        expected_nuc_counts = '4,0,5,0'

        self.amino.count_aminos(nuc_seq1, 4)
        self.amino.count_aminos(nuc_seq2, 5)
        counts = self.amino.nucleotides[2].get_report()

        self.assertSequenceEqual(expected_nuc_counts, counts)

    def testConsensus(self):
        nuc_seq1 = 'AAA'  # -> K
        nuc_seq2 = 'GGG'  # -> G
        expected_consensus = 'G'

        self.amino.count_aminos(nuc_seq1, 4)
        self.amino.count_aminos(nuc_seq2, 5)
        consensus = self.amino.get_consensus()

        self.assertSequenceEqual(expected_consensus, consensus)

    def testConsensusMixture(self):
        nuc_seq1 = 'AAA'  # -> K
        nuc_seq2 = 'GGG'  # -> G
        nuc_seq3 = 'TTT'  # -> F
        allowed_consensus_values = ('G', 'K')

        self.amino.count_aminos(nuc_seq1, 4)
        self.amino.count_aminos(nuc_seq2, 4)
        self.amino.count_aminos(nuc_seq3, 3)
        consensus = self.amino.get_consensus()

        self.assertIn(consensus, allowed_consensus_values)

    def testConsensusWithNoReads(self):
        consensus = self.amino.get_consensus()

        self.assertEqual(consensus, '-')

    def testMissingData(self):
        """ Lower-case n represents a gap between the forward and reverse reads. """

        nuc_seq = 'CTn'
        expected_consensus = '-'

        self.amino.count_aminos(nuc_seq, 1)
        consensus = self.amino.get_consensus()

        self.assertEqual(expected_consensus, consensus)

    def testAmbiguousData(self):
        """If a read is ambiguous, don't count it toward consensus."""

        nuc_seq1 = 'Cnn'  # -> ?
        nuc_seq2 = 'AAA'  # -> K
        expected_consensus = 'K'

        self.amino.count_aminos(nuc_seq1, 9)
        self.amino.count_aminos(nuc_seq2, 1)
        consensus = self.amino.get_consensus()

        self.assertEqual(expected_consensus, consensus)

    def testOverlap(self):
        self.amino.count_aminos('GGG', 4)
        other = SeedAmino(consensus_nuc_index=7)
        other.count_aminos('TAG', 5)
        expected_counts = {'G': 4}
        expected_v3_overlap = 5

        self.amino.count_overlap(other)

        self.assertEqual(expected_counts, self.amino.counts)
        self.assertEqual(expected_v3_overlap, self.amino.v3_overlap)
        self.assertEqual(expected_v3_overlap,
                         self.amino.nucleotides[0].v3_overlap)

    def testOverlapPartialCodon(self):
        self.amino.count_aminos('GGG', 4)
        other = SeedAmino(consensus_nuc_index=7)
        other.count_aminos('TA', 5)
        expected_counts = {'G': 4}
        expected_v3_overlap = 5

        self.amino.count_overlap(other)

        self.assertEqual(expected_counts, self.amino.counts)
        self.assertEqual(expected_v3_overlap, self.amino.v3_overlap)
        self.assertEqual(expected_v3_overlap,
                         self.amino.nucleotides[0].v3_overlap)


class SeedNucleotideTest(unittest.TestCase):
    def setUp(self):
        self.nuc = SeedNucleotide()

    def testSingleRead(self):
        """ Read a single nucleotide, and report on counts.
        Columns are:       A,C,G,T
        """
        nuc_seq = 'C'
        expected_counts = '0,8,0,0'

        self.nuc.count_nucleotides(nuc_seq, 8)
        counts = self.nuc.get_report()

        self.assertSequenceEqual(expected_counts, counts)

    def testConsensusNoMixes(self):
        self.nuc.count_nucleotides('C', 1)
        consensus_max = self.nuc.get_consensus(MAX_CUTOFF)
        consensus_mix = self.nuc.get_consensus(0.1)

        expected_consensus = 'C'
        self.assertEqual(expected_consensus, consensus_max)
        self.assertEqual(expected_consensus, consensus_mix)

    def testConsensusMixed(self):
        self.nuc.count_nucleotides('C', 2)
        self.nuc.count_nucleotides('T', 1)
        consensus_max = self.nuc.get_consensus(MAX_CUTOFF)
        consensus_mix = self.nuc.get_consensus(0.1)

        expected_consensus_max = 'C'
        expected_consensus_mix = 'Y'
        self.assertEqual(expected_consensus_max, consensus_max)
        self.assertEqual(expected_consensus_mix, consensus_mix)

    def testConsensusMixedThree(self):
        self.nuc.count_nucleotides('C', 2)
        self.nuc.count_nucleotides('T', 1)
        self.nuc.count_nucleotides('G', 1)
        consensus_max = self.nuc.get_consensus(MAX_CUTOFF)
        consensus_mix = self.nuc.get_consensus(0.1)

        expected_consensus_max = 'C'
        expected_consensus_mix = 'B'  # B is a mix of T, G, and C
        self.assertEqual(expected_consensus_max, consensus_max)
        self.assertEqual(expected_consensus_mix, consensus_mix)

    def testConsensusMixedAll(self):
        self.nuc.count_nucleotides('C', 2)
        self.nuc.count_nucleotides('T', 1)
        self.nuc.count_nucleotides('G', 1)
        self.nuc.count_nucleotides('A', 1)
        consensus_max = self.nuc.get_consensus(MAX_CUTOFF)
        consensus_mix = self.nuc.get_consensus(0.1)

        expected_consensus_max = 'C'
        expected_consensus_mix = 'N'  # All four are reported as N
        self.assertEqual(expected_consensus_max, consensus_max)
        self.assertEqual(expected_consensus_mix, consensus_mix)

    def testConsensusMixedMax(self):
        self.nuc.count_nucleotides('C', 2)
        self.nuc.count_nucleotides('T', 2)
        self.nuc.count_nucleotides('G', 1)
        consensus_max = self.nuc.get_consensus(MAX_CUTOFF)
        consensus_mix = self.nuc.get_consensus(0.1)

        expected_consensus_max = 'Y'  # C and T tie for max, mix is Y
        expected_consensus_mix = 'B'  # C, T, and G mix is B
        self.assertEqual(expected_consensus_max, consensus_max)
        self.assertEqual(expected_consensus_mix, consensus_mix)

    def testConsensusCutoff(self):
        self.nuc.count_nucleotides('C', 2)
        self.nuc.count_nucleotides('T', 1)
        consensus_mix = self.nuc.get_consensus(0.5)

        expected_consensus = 'C'  # T was below the cutoff
        self.assertEqual(expected_consensus, consensus_mix)

    def testConsensusCutoffAtBoundary(self):
        self.nuc.count_nucleotides('C', 9000)
        self.nuc.count_nucleotides('T', 1000)
        consensus_mix = self.nuc.get_consensus(0.1)

        expected_consensus = 'Y'  # T was at the cutoff
        self.assertEqual(expected_consensus, consensus_mix)

    def testConsensusCutoffBelowBoundary(self):
        self.nuc.count_nucleotides('C', 9001)
        self.nuc.count_nucleotides('T', 999)
        consensus_mix = self.nuc.get_consensus(0.5)

        expected_consensus = 'C'  # T was below the cutoff
        self.assertEqual(expected_consensus, consensus_mix)

    def testConsensusMixedWithPoorQuality(self):
        self.nuc.count_nucleotides('N', 2)
        self.nuc.count_nucleotides('T', 1)
        consensus_max = self.nuc.get_consensus(MAX_CUTOFF)
        consensus_mix = self.nuc.get_consensus(0.1)

        expected_consensus_max = 'T'  # N always overruled
        expected_consensus_mix = 'T'
        self.assertEqual(expected_consensus_max, consensus_max)
        self.assertEqual(expected_consensus_mix, consensus_mix)

    def testConsensusMixedWithGap(self):
        self.nuc.count_nucleotides('-', 2)
        self.nuc.count_nucleotides('T', 1)
        consensus_max = self.nuc.get_consensus(MAX_CUTOFF)
        consensus_mix = self.nuc.get_consensus(0.1)

        expected_consensus_max = 'T'  # dash always overruled
        expected_consensus_mix = 'T'
        self.assertEqual(expected_consensus_max, consensus_max)
        self.assertEqual(expected_consensus_mix, consensus_mix)

    def testConsensusMixedWithGapAndPoorQuality(self):
        self.nuc.count_nucleotides('N', 3)
        self.nuc.count_nucleotides('-', 2)
        self.nuc.count_nucleotides('T', 1)
        consensus_max = self.nuc.get_consensus(MAX_CUTOFF)
        consensus_mix = self.nuc.get_consensus(0.1)

        expected_consensus_max = 'T'
        expected_consensus_mix = 'T'
        self.assertEqual(expected_consensus_max, consensus_max)
        self.assertEqual(expected_consensus_mix, consensus_mix)

    def testConsensusPoorQualityOnly(self):
        self.nuc.count_nucleotides('N', 1)
        consensus_max = self.nuc.get_consensus(MAX_CUTOFF)
        consensus_mix = self.nuc.get_consensus(0.1)

        expected_consensus_max = 'N'
        expected_consensus_mix = 'N'
        self.assertEqual(expected_consensus_max, consensus_max)
        self.assertEqual(expected_consensus_mix, consensus_mix)

    def testConsensusMixedGapAndPoorQualityOnly(self):
        self.nuc.count_nucleotides('N', 3)
        self.nuc.count_nucleotides('-', 2)
        consensus_max = self.nuc.get_consensus(MAX_CUTOFF)
        consensus_mix = self.nuc.get_consensus(0.1)

        expected_consensus_max = 'N'
        expected_consensus_mix = 'N'
        self.assertEqual(expected_consensus_max, consensus_max)
        self.assertEqual(expected_consensus_mix, consensus_mix)

    def testConsensusAllBelowCutoff(self):
        self.nuc.count_nucleotides('C', 101)
        self.nuc.count_nucleotides('T', 100)
        self.nuc.count_nucleotides('G', 99)
        consensus_max = self.nuc.get_consensus(MAX_CUTOFF)
        consensus_mix = self.nuc.get_consensus(0.5)

        expected_consensus_max = 'C'
        expected_consensus_mix = 'N'
        self.assertEqual(expected_consensus_max, consensus_max)
        self.assertEqual(expected_consensus_mix, consensus_mix)

    def testConsensusBetweenReads(self):
        """Lower-case n represents the gap between forward and reverse reads.

        Should not be counted in consensus totals"""
        self.nuc.count_nucleotides('C', 9)
        self.nuc.count_nucleotides('T', 1)
        self.nuc.count_nucleotides('n', 2)
        consensus_mix = self.nuc.get_consensus(0.1)

        expected_consensus = 'Y'
        self.assertEqual(expected_consensus, consensus_mix)

    def testConsensusMissingPositions(self):
        """ Positions that are never read are ignored in the consensus. """

        # No counts added

        consensus_max = self.nuc.get_consensus(MAX_CUTOFF)
        consensus_mix = self.nuc.get_consensus(0.1)

        expected_consensus = ''
        self.assertEqual(expected_consensus, consensus_max)
        self.assertEqual(expected_consensus, consensus_mix)

    def testOverlap(self):
        self.nuc.count_nucleotides('T', 4)
        other = SeedNucleotide()
        other.count_nucleotides('C', 5)
        expected_counts = {'T': 4}
        expected_v3_overlap = 5

        self.nuc.count_overlap(other)

        self.assertEqual(expected_counts, self.nuc.counts)
        self.assertEqual(expected_v3_overlap, self.nuc.v3_overlap)
