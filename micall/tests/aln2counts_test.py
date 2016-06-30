import csv
import StringIO
import sys
import unittest

from micall.core.aln2counts import SequenceReport, SeedNucleotide,\
    InsertionWriter, MAX_CUTOFF, SeedAmino
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
                     score=sys.maxint):
        self.overrides[(reference, query)] = (aligned_reference,
                                              aligned_query,
                                              score)


class SequenceReportTest(unittest.TestCase):
    def setUp(self):
        self.insertion_file = StringIO.StringIO()
        insert_writer = InsertionWriter(
            insert_file=self.insertion_file)
        projects = project_config.ProjectConfig()

        # Content of seed regions is irrelevant. For R-NO-COORD, there is
        # no coordinate reference, so we use the seed reference for display, but
        # only the length matters.
        projects.load(StringIO.StringIO("""\
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
    }
  },
  "regions": {
    "R1-seed": {
      "is_nucleotide": true,
      "reference": [
        "A"
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
        "A"
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
        "A"
      ]
    },
    "R3": {
      "is_nucleotide": false,
      "reference": [
        "KFQTPREH"
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
        self.report_file = StringIO.StringIO()

    def prepareReads(self, aligned_reads_text):
        full_text = "refname,qcut,rank,count,offset,seq\n" + aligned_reads_text
        dummy_file = StringIO.StringIO(full_text)
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
        self.report.write_consensus(self.report_file)

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
        self.report.write_consensus(self.report_file)

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
        self.report.write_consensus(self.report_file)

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
        self.report.write_consensus(self.report_file)

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
        self.report.write_consensus(self.report_file)

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
seed,region,q-cutoff,query.aa.pos,refseq.aa.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,2,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
"""

        self.report.write_amino_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_amino_counts(self.report_file)

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
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T
R1-seed,R1,15,1,1,9,0,0,0
R1-seed,R1,15,2,2,9,0,0,0
R1-seed,R1,15,3,3,9,0,0,0
R1-seed,R1,15,4,4,0,0,0,9
R1-seed,R1,15,5,5,0,0,0,9
R1-seed,R1,15,6,6,0,0,0,9
R1-seed,R1,15,,7,0,0,0,0
R1-seed,R1,15,,8,0,0,0,0
R1-seed,R1,15,,9,0,0,0,0
"""

        self.report.write_nuc_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_nuc_counts(self.report_file)

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
        self.report.write_amino_counts(self.report_file,
                                       coverage_summary=summary)

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
        self.report.write_amino_counts(self.report_file,
                                       coverage_summary=summary)

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
        self.report.write_amino_counts(self.report_file,
                                       coverage_summary=summary)

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

        # seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T
        expected_text = """\
R1-seed,R1,15,1,1,0,0,0,0
R1-seed,R1,15,2,2,0,0,0,0
R1-seed,R1,15,3,3,0,0,0,0
R1-seed,R1,15,4,4,0,0,0,1
R1-seed,R1,15,5,5,0,0,0,1
R1-seed,R1,15,6,6,0,0,0,9
R1-seed,R1,15,7,7,0,8,0,0
R1-seed,R1,15,8,8,0,0,8,0
R1-seed,R1,15,9,9,8,0,0,0
"""

        self.report.read(aligned_reads)
        self.report.write_nuc_counts(self.report_file)

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testPartialCodonNucleotideReport(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,9,0,AAATT
""")

        # seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T
        expected_text = """\
R1-seed,R1,15,1,1,9,0,0,0
R1-seed,R1,15,2,2,9,0,0,0
R1-seed,R1,15,3,3,9,0,0,0
R1-seed,R1,15,4,4,0,0,0,9
R1-seed,R1,15,5,5,0,0,0,9
R1-seed,R1,15,6,6,0,0,0,0
R1-seed,R1,15,,7,0,0,0,0
R1-seed,R1,15,,8,0,0,0,0
R1-seed,R1,15,,9,0,0,0,0
"""

        self.report.read(aligned_reads)
        self.report.write_nuc_counts(self.report_file)

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

        # seed,region,q-cutoff,query.aa.pos,refseq.aa.pos,
        #         A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
        expected_text = """\
R1-seed,R1,15,2,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,3,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,4,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0
"""

        self.report.read(aligned_reads)
        self.report.write_amino_counts(self.report_file)

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

        # seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T
        expected_text = """\
R2-seed,R2,15,1,1,9,0,0,0
R2-seed,R2,15,2,2,9,0,0,0
R2-seed,R2,15,3,3,9,0,0,0
R2-seed,R2,15,4,4,0,0,0,9
R2-seed,R2,15,5,5,0,0,0,9
R2-seed,R2,15,6,6,0,0,0,9
R2-seed,R2,15,,7,0,0,0,0
R2-seed,R2,15,,8,0,0,0,0
R2-seed,R2,15,,9,0,0,0,0
R2-seed,R2,15,7,10,0,9,0,0
R2-seed,R2,15,8,11,0,9,0,0
R2-seed,R2,15,9,12,0,9,0,0
R2-seed,R2,15,10,13,0,9,0,0
R2-seed,R2,15,11,14,0,0,9,0
R2-seed,R2,15,12,15,9,0,0,0
"""

        self.report.read(aligned_reads)
        self.report.write_nuc_counts(self.report_file)

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

        # seed,region,q-cutoff,query.aa.pos,refseq.aa.pos,
        #         A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
        expected_text = """\
R2-seed,R2,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,2,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,3,4,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0
R2-seed,R2,15,4,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0
"""

        self.report.read(aligned_reads)
        self.report.write_amino_counts(self.report_file)

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

        # seed,region,q-cutoff,query.aa.pos,refseq.aa.pos,
        #         A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
        expected_text = """\
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,7,0,0,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,2,2,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,3,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,0,0,0,0,0,0
"""

        self.report.read(aligned_reads)
        self.report.write_amino_counts(self.report_file)

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

        # seed,region,q-cutoff,query.aa.pos,refseq.aa.pos,
        #         A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
        expected_text = """\
R3-seed,R3,15,4,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0
R3-seed,R3,15,5,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R3-seed,R3,15,6,3,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0
R3-seed,R3,15,7,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0
R3-seed,R3,15,9,5,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0
R3-seed,R3,15,10,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0
R3-seed,R3,15,11,7,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R3-seed,R3,15,12,8,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0
"""
        expected_insertions = """\
seed,region,qcut,left,insert,count,before
R3-seed,R3,15,8,G,9,5
"""

        self.report.read(aligned_reads)
        self.report.write_amino_counts(self.report_file)
        self.report.write_insertions()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())
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
R3-seed,R3,15,5,G,9,5
"""

        self.report.read(aligned_reads)
        self.report.write_amino_counts(self.report_file)
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
R3-seed,R3,15,5,G,9,5
"""

        self.report.read(aligned_reads)
        self.report.write_amino_counts(self.report_file)
        self.report.write_insertions()

        self.assertMultiLineEqual(expected_insertions,
                                  self.insertion_file.getvalue())

    def testMultipleCoordinateInsertionReport(self):
        """ Two coordinate regions map the same seed region, the consensus
        has an insertion relative to only one of them.
        """
        self.report.projects.load(StringIO.StringIO("""\
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
        "A"
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
R3-seed,R3a,15,5,G,9,5
"""

        self.report.read(aligned_reads)
        self.report.write_insertions()

        self.assertMultiLineEqual(expected_insertions,
                                  self.insertion_file.getvalue())

    def testGapBetweenForwardAndReverse(self):
        """ Lower-case n represents a gap between forward and reverse reads.

        Region R2 has sequence KFGPR, so this read has a gap at the end of G
        and beginning of P. G is still unambiguous, but P is not.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R2-seed,15,0,5,0,AAATTTGGnnCCCGA
""")

        # seed,region,q-cutoff,query.aa.pos,refseq.aa.pos,
        #         A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
        expected_text = """\
R2-seed,R2,15,1,1,0,0,0,0,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,2,2,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,3,3,0,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,4,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,5,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,0,0,0,0,0,0
"""

        self.report.read(aligned_reads)
        self.report.write_amino_counts(self.report_file)

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
        self.report.write_failure(self.report_file)

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testFailedAlignmentWithHeadToTailMatch(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R3-seed,15,0,2,0,TTATCCTACTTATCCTACTTATCCAAA
""")

        expected_text = """\
R3-seed,R3,15,LSYLSYLSK,KFQTPREH
"""

        self.report.read(aligned_reads)
        self.report.write_failure(self.report_file)

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testGoodAlignmentWithTinyCoordinateReference(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,2,0,AAATTTCGATTATCCTACTTATCCTACTTATCCTACTTATCCTACTTATCCTACTTATCCTACTTATCCTACTTATCCTAC
""")

        expected_text = ""

        self.report.read(aligned_reads)
        self.report.write_failure(self.report_file)

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testGoodAlignmentWithGiantSeed(self):
        """ Short consensus with long seed and long coordinate reference.

        Even when the consensus maps to the end of the seed, it should still
        only require a low alignment score.
        """
        self.report.projects.load(StringIO.StringIO("""\
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
        "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
        "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
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
        expected_text = ""

        self.report.read(aligned_reads)
        self.report.write_failure(self.report_file)

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
R1-seed,R1,15,-LSY,KFR
"""

        self.report.read(aligned_reads)
        self.report.write_failure(self.report_file)

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testNoFailureReport(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,9,0,AAATTT
""")

        expected_text = ""

        self.report.read(aligned_reads)
        self.report.write_failure(self.report_file)

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testRegionWithoutCoordinateReferenceNucleotideReport(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R-NO-COORD,15,0,9,0,AAATTT
""")

        # seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T
        expected_text = """\
R-NO-COORD,R-NO-COORD,15,1,,9,0,0,0
R-NO-COORD,R-NO-COORD,15,2,,9,0,0,0
R-NO-COORD,R-NO-COORD,15,3,,9,0,0,0
R-NO-COORD,R-NO-COORD,15,4,,0,0,0,9
R-NO-COORD,R-NO-COORD,15,5,,0,0,0,9
R-NO-COORD,R-NO-COORD,15,6,,0,0,0,9
R-NO-COORD,R-NO-COORD,15,7,,0,0,0,0
R-NO-COORD,R-NO-COORD,15,8,,0,0,0,0
R-NO-COORD,R-NO-COORD,15,9,,0,0,0,0
"""

        self.report.read(aligned_reads)
        self.report.write_nuc_counts(self.report_file)

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testRegionWithoutCoordinateReferenceFailureReport(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R-NO-COORD,15,0,9,0,AAATTT
""")

        expected_text = ""

        self.report.read(aligned_reads)
        self.report.write_failure(self.report_file)

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testMultipleCoordinateAminoReport(self):
        """ Two coordinate regions map the same seed region, report both.
        """
        self.report.projects.load(StringIO.StringIO("""\
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
        "A"
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
R1-seed,R1a,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0
R1-seed,R1a,15,2,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R1-seed,R1a,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R1-seed,R1b,15,,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R1-seed,R1b,15,1,2,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0
R1-seed,R1b,15,2,3,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R1-seed,R1b,15,,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
"""

        self.report.read(aligned_reads)
        self.report.write_amino_counts(self.report_file)

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testVariantReport(self):
        """ Reference is KFR, variants are KFRx10, GFRx8, KFGx6, KFSx5
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,10,0,AAATTTCGA
R1-seed,15,1,8,0,GGGTTTCGA
R1-seed,15,2,6,0,AAATTTGGG
R1-seed,15,3,5,0,AAATTTTCA
""")

        expected_text = """\
seed,qcut,region,index,count,seq
R1-seed,15,R1,0,10,AAATTTCGA
R1-seed,15,R1,1,8,GGGTTTCGA
R1-seed,15,R1,2,6,AAATTTGGG
R1-seed,15,R1,3,5,AAATTTTCA
"""

        self.report.write_nuc_variants_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_nuc_variants(self.report_file)

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testVariantLimit(self):
        """ Only report top 3 variants
        """
        self.report.projects.load(StringIO.StringIO("""\
{
  "projects": {
    "R1": {
      "max_variants": 3,
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
        "A"
      ]
    },
    "R1": {
      "is_nucleotide": false,
      "reference": [
        "KFR"
      ]
    }
  }
}
"""))

        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,10,0,AAATTTCGA
R1-seed,15,1,8,0,GGGTTTCGA
R1-seed,15,2,6,0,AAATTTGGG
R1-seed,15,3,5,0,AAATTTTCA
""")

        # seed,qcut,region,index,count,seq
        expected_text = """\
R1-seed,15,R1,0,10,AAATTTCGA
R1-seed,15,R1,1,8,GGGTTTCGA
R1-seed,15,R1,2,6,AAATTTGGG
"""
        self.report.read(aligned_reads)
        self.report.write_nuc_variants(self.report_file)

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testVariantClippingEnd(self):
        """ Coordinate reference is KFR, so last codons should be clipped
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,10,0,AAATTTCGATCA
R1-seed,15,1,8,0,GGGTTTCGATCA
R1-seed,15,2,6,0,AAATTTGGGTCA
R1-seed,15,3,5,0,AAATTTTCA
""")

        # seed,qcut,region,index,count,seq
        expected_text = """\
R1-seed,15,R1,0,10,AAATTTCGA
R1-seed,15,R1,1,8,GGGTTTCGA
R1-seed,15,R1,2,6,AAATTTGGG
R1-seed,15,R1,3,5,AAATTTTCA
"""

        self.report.read(aligned_reads)
        self.report.write_nuc_variants(self.report_file)

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testVariantClippingStart(self):
        """ Coordinate reference is KFR, so first codons should be clipped
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,10,0,TCAAAATTTCGA
R1-seed,15,1,8,0,TCAGGGTTTCGA
R1-seed,15,2,6,0,TCAAAATTTGGG
R1-seed,15,3,5,0,TCAAAATTTTCA
""")

        # seed,qcut,region,index,count,seq
        expected_text = """\
R1-seed,15,R1,0,10,AAATTTCGA
R1-seed,15,R1,1,8,GGGTTTCGA
R1-seed,15,R1,2,6,AAATTTGGG
R1-seed,15,R1,3,5,AAATTTTCA
"""

        self.report.read(aligned_reads)
        self.report.write_nuc_variants(self.report_file)

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testVariantWhenAlignFails(self):
        """ Coordinate reference is KFR, read is SLS.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,10,0,TCACTCTCT
""")

        expected_text = ""

        self.report.read(aligned_reads)
        self.report.write_nuc_variants(self.report_file)

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testVariantWhenSomeAlignsFail(self):
        """ Coordinate reference and first read are KFR, second read is SLS,
        mapped to another part of the seed.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,10,0,AAATTTCGA
R1-seed,15,1,20,30,TCACTCTCT
""")

        expected_text = """\
R1-seed,15,R1,0,10,AAATTTCGA
"""

        self.report.read(aligned_reads)
        self.report.write_nuc_variants(self.report_file)

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testVariantShortRead(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,10,0,AAATTT
R1-seed,15,1,8,0,GGGTTT
""")

        # seed,qcut,region,index,count,seq
        expected_text = """\
R1-seed,15,R1,0,10,AAATTT
R1-seed,15,R1,1,8,GGGTTT
"""

        self.report.read(aligned_reads)
        self.report.write_nuc_variants(self.report_file)

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testVariantOffset(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,10,0,AAATTT
R1-seed,15,1,8,3,TTT
""")

        # seed,qcut,region,index,count,seq
        expected_text = """\
R1-seed,15,R1,0,10,AAATTT
R1-seed,15,R1,1,8,---TTT
"""

        self.report.read(aligned_reads)
        self.report.write_nuc_variants(self.report_file)

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testVariantsCombinedAfterClipping(self):
        """ If reads only differ outside the clipped region, their counts
        should be combined. max_variants is applied after combining the clipped
        reads.
        """
        self.report.projects.load(StringIO.StringIO("""\
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
        "A"
      ]
    },
    "R1": {
      "is_nucleotide": false,
      "reference": [
        "KFR"
      ]
    }
  }
}
"""))
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,10,0,TCAAAATTTCGA
R1-seed,15,1,9,0,TCAAAGTTTCGA
R1-seed,15,2,7,0,TCAAAGTTTGGG
R1-seed,15,3,5,0,CGAAAGTTTCGA
""")

        # seed,qcut,region,index,count,seq
        expected_text = """\
R1-seed,15,R1,0,14,AAGTTTCGA
R1-seed,15,R1,1,10,AAATTTCGA
"""
        self.report.read(aligned_reads)
        self.report.write_nuc_variants(self.report_file)

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testVariantInsertion(self):
        """ Coordinate reference is KFR, consensus is KFGR, and alignment decides
        on KF-R, KFGR. The G is an insertion inside the coordinate reference
        and should not be clipped.
        """
        self.report.add_override(reference='KFR',
                                 query='KFGR',
                                 aligned_query='KFGR',
                                 aligned_reference='KF-R')
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,10,0,AAATTTGGGCGA
""")

        # seed,qcut,region,index,count,seq
        expected_text = """\
R1-seed,15,R1,0,10,AAATTTGGGCGA
"""

        self.report.read(aligned_reads)
        self.report.write_nuc_variants(self.report_file)

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testVariantInsertionPastEnd(self):
        """ Coordinate reference is KFR, consensus is KFG, and alignment decides
        on KFR-, KF-G. The G is an insertion outside the coordinate reference
        and should be clipped.
        """
        self.report.add_override(reference='KFR',
                                 query='KFG',
                                 aligned_query='KF-G',
                                 aligned_reference='KFR-')
        # refname,qcut,rank,count,offset,seq
        aligned_reads = self.prepareReads("""\
R1-seed,15,0,10,0,AAATTTGGG
""")

        # seed,qcut,region,index,count,seq
        expected_text = """\
R1-seed,15,R1,0,10,AAATTT
"""

        self.report.read(aligned_reads)
        self.report.write_nuc_variants(self.report_file)

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testMultipleCoordinateVariantReport(self):
        """ Two coordinate regions map the same seed region, report both.
        """
        self.report.projects.load(StringIO.StringIO("""\
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
        "A"
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

        # seed,qcut,region,index,count,seq
        expected_text = """\
R1-seed,15,R1a,0,9,AAATTT
R1-seed,15,R1b,0,9,AAATTT
"""

        self.report.read(aligned_reads)
        self.report.write_nuc_variants(self.report_file)

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())


class InsertionWriterTest(unittest.TestCase):
    def setUp(self):
        self.insert_file = StringIO.StringIO()
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
R1-seed,R1,15,3,D,1,
"""

        self.writer.add_nuc_read(offset_sequence=self.nuc_seq_acdef, count=1)
        self.writer.write(inserts=[2], region='R1')

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())

    def testInsertDifferentReadingFrame(self):
        """ Add a partial codon at the start of the read to shift the reading
        frame.
        """
        expected_text = """\
seed,region,qcut,left,insert,count,before
R1-seed,R1,15,4,D,1,
"""

        self.writer.add_nuc_read(offset_sequence='A' + self.nuc_seq_acdef,
                                 count=1)
        self.writer.write(inserts=[3], region='R1', reading_frame=2)

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())

    def testInsertWithOffset(self):
        expected_text = """\
seed,region,qcut,left,insert,count,before
R1-seed,R1,15,3,D,1,
"""

        #                                            C  D  E  F
        self.writer.add_nuc_read(offset_sequence='---TGTGACGAGTTT', count=1)
        self.writer.write(inserts=[2], region='R1')

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())

    def testInsertWithDeletion(self):
        expected_text = """\
seed,region,qcut,left,insert,count,before
"""

        #                                         C  D     E  F
        self.writer.add_nuc_read(offset_sequence='TGTGAC---GAGTTT', count=1)
        self.writer.write(inserts=[1, 2], region='R1')

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())

    def testTwoInsertsWithOffset(self):
        expected_text = """\
seed,region,qcut,left,insert,count,before
R1-seed,R1,15,3,D,1,
R1-seed,R1,15,5,F,1,
"""

        #                                            C  D  E  F  G
        self.writer.add_nuc_read(offset_sequence='---TGTGACGAGTTTGGG', count=1)
        self.writer.write(inserts=[2, 4], region='R1')

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())

    def testInsertsWithVariants(self):
        expected_text = """\
seed,region,qcut,left,insert,count,before
R1-seed,R1,15,3,D,2,
"""

        self.writer.add_nuc_read(offset_sequence=self.nuc_seq_acdef, count=1)
        self.writer.add_nuc_read(offset_sequence=self.nuc_seq_afdef, count=1)
        self.writer.write(inserts=[2], region='R1')

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())

    def testDifferentInserts(self):
        expected_text = """\
seed,region,qcut,left,insert,count,before
R1-seed,R1,15,2,C,2,
R1-seed,R1,15,2,F,3,
"""

        self.writer.add_nuc_read(offset_sequence=self.nuc_seq_acdef, count=2)
        self.writer.add_nuc_read(offset_sequence=self.nuc_seq_afdef, count=3)
        self.writer.write(inserts=[1], region='R1')

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())

    def testMulticharacterInsert(self):
        expected_text = """\
seed,region,qcut,left,insert,count,before
R1-seed,R1,15,3,DE,1,
"""

        self.writer.add_nuc_read(offset_sequence=self.nuc_seq_acdef, count=1)
        self.writer.write(inserts=[2, 3], region='R1')

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())

    def testReadGapInInsert(self):
        nuc_seq = 'GCTCTnGACGAGTTT'

        expected_text = """\
seed,region,qcut,left,insert,count,before
"""

        self.writer.add_nuc_read(nuc_seq, count=1)
        self.writer.write(inserts=[1], region='R1')

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())

    def testUnsortedInserts(self):
        expected_text = """\
seed,region,qcut,left,insert,count,before
R1-seed,R1,15,3,DE,1,
"""

        self.writer.add_nuc_read(offset_sequence=self.nuc_seq_acdef, count=1)
        self.writer.write(inserts=(3, 2), region='R1')

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
        "Lower-case n represents a gap between the forward and reverse reads."

        nuc_seq = 'CTn'
        expected_consensus = 'L'

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
        "Positions that are never read are ignored in the consensus."

        # No counts added

        consensus_max = self.nuc.get_consensus(MAX_CUTOFF)
        consensus_mix = self.nuc.get_consensus(0.1)

        expected_consensus = ''
        self.assertEqual(expected_consensus, consensus_max)
        self.assertEqual(expected_consensus, consensus_mix)
