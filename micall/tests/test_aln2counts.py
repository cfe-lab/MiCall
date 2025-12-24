from io import StringIO
import unittest

import yaml

from micall.core.aln2counts import InsertionWriter, SeedAmino, \
    ReportAmino, ConsensusBuilder, ReportNucleotide, SeedNucleotide
from micall.tests.test_aln2counts_report import create_sequence_report, prepare_reads

LANDMARKS_YAML = """\
- seed_pattern: R1-.*
  coordinates: R1-seed
  landmarks:
    # Extra 3 positions for stop codon to get dropped.
    - {name: R1, start: 1, end: 12, colour: steelblue}
"""


# noinspection DuplicatedCode
class SequenceReportTest(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        self.report = create_sequence_report()
        self.report_file = StringIO()
        self.detail_report_file = StringIO()

    def testEmptyAminoReport(self):
        expected_text = ""

        self.report.write_amino_counts(self.report_file)

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testEmptyNucReport(self):
        expected_text = ""

        self.report.write_nuc_counts(self.report_file)

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testConsensusFromTwoReads(self):
        """ The second read is out voted by the first one.
        CCC -> P
        GGG -> G
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
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

    def testConsensusExactTie(self):
        """ There is an exact tie between sequences. """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R1-seed,15,0,5,0,AAATTT
R1-seed,15,0,5,0,CCCGGG
""")
        expected_text = """\
region,q-cutoff,consensus-percent-cutoff,offset,sequence
R1-seed,15,MAX,0,MMMKKK
R1-seed,15,0.100,0,MMMKKK
"""

        self.report.write_consensus_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_consensus()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testConsensusWithOffset(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
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

    def testConsensusFromPartialContig(self):
        """ Contigs with the -partial suffix report consensus. """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
1-R2-seed-partial,15,0,9,0,AAATTT
""")
        expected_text = """\
region,q-cutoff,consensus-percent-cutoff,offset,sequence
1-R2-seed-partial,15,MAX,0,AAATTT
1-R2-seed-partial,15,0.100,0,AAATTT
"""

        self.report.write_consensus_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_consensus()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testConsensusLowQualitySections(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R1-seed,15,0,9,3,NNNTTT
R1-seed,15,0,1,7,TTNGG
""")
        expected_text = """\
region,q-cutoff,consensus-percent-cutoff,offset,sequence
R1-seed,15,MAX,6,TTTxGG
R1-seed,15,0.100,6,TTTxGG
"""
        self.report.consensus_min_coverage = 1
        self.report.consensus_builder.consensus_min_coverage = 1

        self.report.write_consensus_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_consensus()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testConsensusLowQuality(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R1-seed,15,0,9,3,NNNNNN
R1-seed,15,0,1,7,NNNNN
""")
        expected_text = """\
region,q-cutoff,consensus-percent-cutoff,offset,sequence
"""
        self.report.consensus_min_coverage = 1
        self.report.consensus_builder.consensus_min_coverage = 1

        self.report.write_consensus_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_consensus()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testConsensusLowCoverageInMiddle(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R1-seed,15,0,9,0,AAATTTGGG
R1-seed,15,0,1,0,AAAT
R1-seed,15,0,1,6,GGG
""")
        expected_text = """\
region,q-cutoff,consensus-percent-cutoff,offset,sequence
R1-seed,15,MAX,0,AAATxxGGG
R1-seed,15,0.100,0,AAATxxGGG
"""
        self.report.consensus_min_coverage = 10
        self.report.consensus_builder.consensus_min_coverage = 10

        self.report.write_consensus_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_consensus()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testConsensusLowCoverageAtStart(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R1-seed,15,0,9,0,AAATTTGGG
R1-seed,15,0,1,4,TTGGG
""")
        expected_text = """\
region,q-cutoff,consensus-percent-cutoff,offset,sequence
R1-seed,15,MAX,4,TTGGG
R1-seed,15,0.100,4,TTGGG
"""
        self.report.consensus_min_coverage = 10
        self.report.consensus_builder.consensus_min_coverage = 10

        self.report.write_consensus_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_consensus()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testConsensusLowCoverageAtEnd(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R1-seed,15,0,9,0,AAATTTGGG
R1-seed,15,0,1,0,AAAT
""")
        expected_text = """\
region,q-cutoff,consensus-percent-cutoff,offset,sequence
R1-seed,15,MAX,0,AAAT
R1-seed,15,0.100,0,AAAT
"""
        self.report.consensus_min_coverage = 10
        self.report.consensus_builder.consensus_min_coverage = 10

        self.report.write_consensus_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_consensus()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testConsensusAllFromTwoReads(self):
        """ The second read is out voted by the first one.
        CCC -> P
        GGG -> G
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R1-seed,15,0,9,0,AAATTT
R1-seed,15,0,1,0,CCCGGG
""")
        expected_text = """\
seed,region,q-cutoff,consensus-percent-cutoff,seed-offset,region-offset,sequence
R1-seed,,15,MAX,0,,AAATTT
R1-seed,R1,15,MAX,0,0,AAATTT
"""

        self.report.write_consensus_all_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_consensus_all()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testConsensusAllExactTie(self):
        """ Exact ties still result in mixtures. """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R1-seed,15,0,5,0,AAATTT
R1-seed,15,0,5,0,CCCGGG
""")
        expected_text = """\
seed,region,q-cutoff,consensus-percent-cutoff,seed-offset,region-offset,sequence
R1-seed,,15,MAX,0,,MMMKKK
R1-seed,R1,15,MAX,0,0,MMMKKK
"""

        self.report.write_consensus_all_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_consensus_all()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testConsensusAllWithOffset(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R1-seed,15,0,9,3,AAATTT
R1-seed,15,0,1,7,TTGGG
""")
        expected_text = """\
seed,region,q-cutoff,consensus-percent-cutoff,seed-offset,region-offset,sequence
R1-seed,,15,MAX,3,,AAATTTGGG
R1-seed,R1,15,MAX,3,0,AAATTTGGG
"""

        self.report.write_consensus_all_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_consensus_all()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testConsensusAllLowQualitySections(self):
        """Low-quality bases still get reported as x."""
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R1-seed,15,0,9,3,NNNTTT
R1-seed,15,0,1,7,TTNGG
""")
        expected_text = """\
seed,region,q-cutoff,consensus-percent-cutoff,seed-offset,region-offset,sequence
R1-seed,,15,MAX,6,,TTTxGG
R1-seed,R1,15,MAX,6,3,TTTxGG
"""
        self.report.consensus_min_coverage = 1
        self.report.consensus_builder.consensus_min_coverage = 1

        self.report.write_consensus_all_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_consensus_all()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testConsensusAllLowQuality(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R1-seed,15,0,9,3,NNNNNN
R1-seed,15,0,1,7,NNNNN
""")
        expected_text = """\
seed,region,q-cutoff,consensus-percent-cutoff,seed-offset,region-offset,sequence
"""
        self.report.consensus_min_coverage = 1
        self.report.consensus_builder.consensus_min_coverage = 1

        self.report.write_consensus_all_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_consensus_all()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testConsensusAllLowCoverageInMiddle(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R1-seed,15,0,9,0,AAATTTGGG
R1-seed,15,0,1,0,AAAT
R1-seed,15,0,1,6,GGG
""")
        expected_text = """\
seed,region,q-cutoff,consensus-percent-cutoff,seed-offset,region-offset,sequence
R1-seed,,15,MAX,0,,AAATTTGGG
R1-seed,R1,15,MAX,0,0,AAATTTGGG
"""
        self.report.consensus_min_coverage = 10
        self.report.consensus_builder.consensus_min_coverage = 10

        self.report.write_consensus_all_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_consensus_all()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testConsensusAllLowCoverageAtStart(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R1-seed,15,0,9,0,AAATTTGGG
R1-seed,15,0,1,4,TTGGG
""")
        expected_text = """\
seed,region,q-cutoff,consensus-percent-cutoff,seed-offset,region-offset,sequence
R1-seed,,15,MAX,0,,AAATTTGGG
R1-seed,R1,15,MAX,0,0,AAATTTGGG
"""
        self.report.consensus_min_coverage = 10
        self.report.consensus_builder.consensus_min_coverage = 10

        self.report.write_consensus_all_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_consensus_all()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testConsensusAllLowCoverageAtEnd(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R1-seed,15,0,9,0,AAATTTGGG
R1-seed,15,0,1,0,AAAT
""")
        expected_text = """\
seed,region,q-cutoff,consensus-percent-cutoff,seed-offset,region-offset,sequence
R1-seed,,15,MAX,0,,AAATTTGGG
R1-seed,R1,15,MAX,0,0,AAATTTGGG
"""
        self.report.consensus_min_coverage = 10
        self.report.consensus_builder.consensus_min_coverage = 10

        self.report.write_consensus_all_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_consensus_all()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testConsensusAllMapToMultipleRegions(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R7-seed,15,0,9,0,AAATTTCAGACCCCACGAGAGCAT
R7-seed,15,0,1,0,AAATTGCAGACCCCACGAGAGCAT
""")
        expected_text = """\
seed,region,q-cutoff,consensus-percent-cutoff,seed-offset,region-offset,sequence
R7-seed,,15,MAX,0,,AAATTTCAGACCCCACGAGAGCAT
R7-seed,R7a,15,MAX,0,0,AAATTTCAG
R7-seed,R7b,15,MAX,15,0,CGAGAGCAT
"""
        self.report.consensus_min_coverage = 10
        self.report.consensus_builder.consensus_min_coverage = 10

        self.report.write_consensus_all_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_consensus_all()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testMultiplePrefixAminoReport(self):
        """ Assemble counts from three contigs to two references.

        Contig 1-R1 AAATTT -> KF
        Contig 2-R2 GGCCCG -> GP
        Contig 3-R1 TTTAGG -> FR

        Contig 1 and 3 should combine into R1 with KFR.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads1 = prepare_reads("1-R1-seed,15,0,5,0,AAATTT")
        aligned_reads2 = prepare_reads("2-R2-seed,15,0,4,0,GGCCCG")
        aligned_reads3 = prepare_reads("3-R1-seed,15,0,2,0,TTTAGG")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5
R1-seed,R1,15,,2,0,0,0,0,7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7
R1-seed,R1,15,4,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,2
R2-seed,R2,15,1,3,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4
R2-seed,R2,15,4,4,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4
"""

        expected_detail_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
1-R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5
1-R1-seed,R1,15,4,2,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5
2-R2-seed,R2,15,1,3,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4
2-R2-seed,R2,15,4,4,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4
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
        aligned_reads1 = prepare_reads("1-R1-seed,15,0,5,0,AAATTT")
        aligned_reads2 = prepare_reads("2-R1-seed,15,0,2,0,TT-AGG")
        aligned_reads3 = prepare_reads("3-R1-seed,15,0,3,0,AAATTTAGG\n"
                                       "3-R1-seed,15,0,1,0,AAA---AGG")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,,2,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,1,0,0,0,9
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,0,0,0,0,0,0,0,0,0,0,0,0,6
"""

        expected_detail_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
1-R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5
1-R1-seed,R1,15,4,2,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5
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

        self.assertEqual(expected_detail_text,
                         self.detail_report_file.getvalue())
        self.assertEqual(expected_text, self.report_file.getvalue())

    def testMultiplePrefixNucleotideReport(self):
        """ Assemble counts from three contigs to two references.

        Contig 1-R1 AAATTT -> KF
        Contig 2-R2 GGCCCG -> GP
        Contig 3-R1 TTTAGG -> FR

        Contig 1 and 3 should combine into R1 with KFR.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads1 = prepare_reads("1-R1-seed,15,0,5,0,AAATTT")
        aligned_reads2 = prepare_reads("2-R2-seed,15,0,4,0,GGCCCG")
        aligned_reads3 = prepare_reads("3-R1-seed,15,0,2,0,TTTAGG")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,genome.pos,\
A,C,G,T,N,del,ins,clip,v3_overlap,coverage,exact_coverage
R1-seed,R1,15,1,1,1,5,0,0,0,0,0,0,0,0,5,
R1-seed,R1,15,2,2,2,5,0,0,0,0,0,0,0,0,5,10
R1-seed,R1,15,3,3,3,5,0,0,0,0,0,0,0,0,5,10
R1-seed,R1,15,,4,4,0,0,0,7,0,0,0,0,0,7,
R1-seed,R1,15,,5,5,0,0,0,7,0,0,0,0,0,7,
R1-seed,R1,15,,6,6,0,0,0,7,0,0,0,0,0,7,
R1-seed,R1,15,4,7,7,2,0,0,0,0,0,0,0,0,2,10
R1-seed,R1,15,5,8,8,0,0,2,0,0,0,0,0,0,2,10
R1-seed,R1,15,6,9,9,0,0,2,0,0,0,0,0,0,2,
R2-seed,R2,15,1,7,7,0,0,4,0,0,0,0,0,0,4,
R2-seed,R2,15,2,8,8,0,0,4,0,0,0,0,0,0,4,
R2-seed,R2,15,3,9,9,0,4,0,0,0,0,0,0,0,4,
R2-seed,R2,15,4,10,10,0,4,0,0,0,0,0,0,0,4,
R2-seed,R2,15,5,11,11,0,4,0,0,0,0,0,0,0,4,
R2-seed,R2,15,6,12,12,0,0,4,0,0,0,0,0,0,4,
"""

        expected_detail_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,genome.pos,\
A,C,G,T,N,del,ins,clip,v3_overlap,coverage,exact_coverage
1-R1-seed,R1,15,1,1,1,5,0,0,0,0,0,0,0,0,5,
1-R1-seed,R1,15,2,2,2,5,0,0,0,0,0,0,0,0,5,
1-R1-seed,R1,15,3,3,3,5,0,0,0,0,0,0,0,0,5,
1-R1-seed,R1,15,4,4,4,0,0,0,5,0,0,0,0,0,5,
1-R1-seed,R1,15,5,5,5,0,0,0,5,0,0,0,0,0,5,
1-R1-seed,R1,15,6,6,6,0,0,0,5,0,0,0,0,0,5,
2-R2-seed,R2,15,1,7,7,0,0,4,0,0,0,0,0,0,4,
2-R2-seed,R2,15,2,8,8,0,0,4,0,0,0,0,0,0,4,
2-R2-seed,R2,15,3,9,9,0,4,0,0,0,0,0,0,0,4,
2-R2-seed,R2,15,4,10,10,0,4,0,0,0,0,0,0,0,4,
2-R2-seed,R2,15,5,11,11,0,4,0,0,0,0,0,0,0,4,
2-R2-seed,R2,15,6,12,12,0,0,4,0,0,0,0,0,0,4,
3-R1-seed,R1,15,1,4,4,0,0,0,2,0,0,0,0,0,2,
3-R1-seed,R1,15,2,5,5,0,0,0,2,0,0,0,0,0,2,
3-R1-seed,R1,15,3,6,6,0,0,0,2,0,0,0,0,0,2,
3-R1-seed,R1,15,4,7,7,2,0,0,0,0,0,0,0,0,2,
3-R1-seed,R1,15,5,8,8,0,0,2,0,0,0,0,0,0,2,
3-R1-seed,R1,15,6,9,9,0,0,2,0,0,0,0,0,0,2,
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

        assert self.detail_report_file.getvalue() == expected_detail_text
        assert self.report_file.getvalue() == expected_text

    def testNucleotideDetailReportOnlyPartials(self):
        """ The only contig is a partial BLAST match, not reported. """
        # refname,qcut,rank,count,offset,seq
        aligned_reads1 = prepare_reads("1-R1-seed-partial,15,0,5,0,AAATTT")
        aligned_reads2 = prepare_reads("2-R2-seed,15,0,4,0,GGCCCG")
        aligned_reads3 = prepare_reads("3-R1-seed,15,0,2,0,TTTAGG")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,genome.pos,\
A,C,G,T,N,del,ins,clip,v3_overlap,coverage,exact_coverage
R2-seed,R2,15,1,7,7,0,0,4,0,0,0,0,0,0,4,
R2-seed,R2,15,2,8,8,0,0,4,0,0,0,0,0,0,4,
R2-seed,R2,15,3,9,9,0,4,0,0,0,0,0,0,0,4,
R2-seed,R2,15,4,10,10,0,4,0,0,0,0,0,0,0,4,
R2-seed,R2,15,5,11,11,0,4,0,0,0,0,0,0,0,4,
R2-seed,R2,15,6,12,12,0,0,4,0,0,0,0,0,0,4,
R1-seed,R1,15,1,4,4,0,0,0,2,0,0,0,0,0,2,
R1-seed,R1,15,2,5,5,0,0,0,2,0,0,0,0,0,2,
R1-seed,R1,15,3,6,6,0,0,0,2,0,0,0,0,0,2,
R1-seed,R1,15,4,7,7,2,0,0,0,0,0,0,0,0,2,
R1-seed,R1,15,5,8,8,0,0,2,0,0,0,0,0,0,2,2
R1-seed,R1,15,6,9,9,0,0,2,0,0,0,0,0,0,2,2
"""

        expected_detail_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,genome.pos,\
A,C,G,T,N,del,ins,clip,v3_overlap,coverage,exact_coverage
2-R2-seed,R2,15,1,7,7,0,0,4,0,0,0,0,0,0,4,
2-R2-seed,R2,15,2,8,8,0,0,4,0,0,0,0,0,0,4,
2-R2-seed,R2,15,3,9,9,0,4,0,0,0,0,0,0,0,4,
2-R2-seed,R2,15,4,10,10,0,4,0,0,0,0,0,0,0,4,
2-R2-seed,R2,15,5,11,11,0,4,0,0,0,0,0,0,0,4,
2-R2-seed,R2,15,6,12,12,0,0,4,0,0,0,0,0,0,4,
3-R1-seed,R1,15,1,4,4,0,0,0,2,0,0,0,0,0,2,
3-R1-seed,R1,15,2,5,5,0,0,0,2,0,0,0,0,0,2,
3-R1-seed,R1,15,3,6,6,0,0,0,2,0,0,0,0,0,2,
3-R1-seed,R1,15,4,7,7,2,0,0,0,0,0,0,0,0,2,
3-R1-seed,R1,15,5,8,8,0,0,2,0,0,0,0,0,0,2,
3-R1-seed,R1,15,6,9,9,0,0,2,0,0,0,0,0,0,2,
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

        assert self.detail_report_file.getvalue() == expected_detail_text
        assert self.report_file.getvalue() == expected_text

    def testAminoReportForPartialContig(self):
        """ Contigs with the -partial suffix shouldn't be reported. """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
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

    def testSoftClippingNucleotideReport(self):
        """ Combine the soft clipping data with the read counts.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R1-seed,15,0,9,2,ATTTA
""")
        clipping = StringIO("""\
refname,pos,count
R1-seed,1,9
R1-seed,2,9
R1-seed,8,9
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,genome.pos,\
A,C,G,T,N,del,ins,clip,v3_overlap,coverage,exact_coverage
R1-seed,R1,15,,1,1,0,0,0,0,0,0,0,9,0,0,
R1-seed,R1,15,,2,2,0,0,0,0,0,0,0,9,0,0,
R1-seed,R1,15,3,3,3,9,0,0,0,0,0,0,0,0,9,
R1-seed,R1,15,4,4,4,0,0,0,9,0,0,0,0,0,9,
R1-seed,R1,15,5,5,5,0,0,0,9,0,0,0,0,0,9,
R1-seed,R1,15,6,6,6,0,0,0,9,0,0,0,0,0,9,
R1-seed,R1,15,7,7,7,9,0,0,0,0,0,0,0,0,9,
R1-seed,R1,15,,8,8,0,0,0,0,0,0,0,9,0,0,
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
        aligned_reads = prepare_reads("""\
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

        assert self.report_file.getvalue() == expected_text

    def testSoftClippingAminoReportMoreOffset(self):
        """ Combine the soft clipping data with the read counts.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R1-seed,15,0,9,5,TTTAGG
""")
        clipping = StringIO("""\
refname,pos,count
R1-seed,3,9
R1-seed,4,9
R1-seed,5,9
R1-seed,11,9
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R1-seed,R1,15,,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0
R1-seed,R1,15,6,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,9,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,9,0,9
"""

        self.report.read_clipping(clipping)
        self.report.write_amino_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_nuc_header(StringIO())
        self.report.write_nuc_counts()
        self.report.write_amino_counts()

        assert self.report_file.getvalue() == expected_text

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
        aligned_reads1 = prepare_reads("1-R1-seed,15,0,5,0,AAATTT")
        aligned_reads2 = prepare_reads("2-R2-seed,15,0,4,0,GGCCCG")
        aligned_reads3 = prepare_reads("3-R1-seed,15,0,2,0,TTTAGG")

        clipping = StringIO("""\
refname,pos,count
1-R1-seed,7,5
3-R1-seed,-1,2
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,5
R1-seed,R1,15,,2,0,0,0,0,7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7
R1-seed,R1,15,4,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,5,0,2
R2-seed,R2,15,1,3,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4
R2-seed,R2,15,4,4,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4
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
        aligned_reads = prepare_reads("""\
R1-seed,15,0,9,0,AAATTT
""")
        conseq_ins_csv = StringIO("""\
qname,fwd_rev,refname,pos,insert,qual
Example_read_1,F,R1-seed,3,AAC,AAA
Example_read_3,F,R2-seed,6,GTA,AAA
Example_read_2,F,R1-seed,3,AAC,AAA
Example_read_2,R,R1-seed,3,AAC,AAA
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,genome.pos,\
A,C,G,T,N,del,ins,clip,v3_overlap,coverage,exact_coverage
R1-seed,R1,15,1,1,1,9,0,0,0,0,0,0,0,0,9,
R1-seed,R1,15,2,2,2,9,0,0,0,0,0,0,0,0,9,18
R1-seed,R1,15,3,3,3,9,0,0,0,0,0,2,0,0,9,18
R1-seed,R1,15,4,4,4,0,0,0,9,0,0,0,0,0,9,18
R1-seed,R1,15,5,5,5,0,0,0,9,0,0,0,0,0,9,18
R1-seed,R1,15,6,6,6,0,0,0,9,0,0,0,0,0,9,
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
        aligned_reads = prepare_reads("""\
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
"""

        self.report.read_insertions(conseq_ins_csv)
        self.report.write_amino_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_nuc_header(StringIO())
        self.report.write_nuc_counts()  # calculates ins counts
        self.report.write_amino_counts()

        assert self.report_file.getvalue() == expected_text

    def testSubstitutionAtBoundary(self):
        """ In this sample, there are nine identical reads with six codons.
        ATG -> M
        GCA -> A
        AAC -> N
        TGG -> W
        ATC -> I
        AAT -> N
        GGG -> G
        The R4 coordinate reference is SING, so its first position will not map.
        However, the ING should map, so the first position should get treated as
        a substitution.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R4-seed,15,0,9,0,ATGGCAAACTGGATCAATGGG
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R4-seed,R4,15,10,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,9
R4-seed,R4,15,13,2,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R4-seed,R4,15,16,3,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R4-seed,R4,15,19,4,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
"""

        self.report.write_amino_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_amino_counts()

        assert self.report_file.getvalue() == expected_text

    def testCoverageSummary(self):
        """ R1 has coverage 9.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
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
        aligned_reads = prepare_reads("""\
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
        aligned_reads = prepare_reads("""\
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
        aligned_reads = prepare_reads("""\
R1-seed,15,0,1,3,TTT
R1-seed,15,0,8,5,TCGA
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,genome.pos,\
A,C,G,T,N,del,ins,clip,v3_overlap,coverage,exact_coverage
R1-seed,R1,15,4,4,4,0,0,0,1,0,0,0,0,0,1,
R1-seed,R1,15,5,5,5,0,0,0,1,0,0,0,0,0,1,
R1-seed,R1,15,6,6,6,0,0,0,9,0,0,0,0,0,9,
R1-seed,R1,15,7,7,7,0,8,0,0,0,0,0,0,0,8,
R1-seed,R1,15,8,8,8,0,0,8,0,0,0,0,0,0,8,
R1-seed,R1,15,9,9,9,8,0,0,0,0,0,0,0,0,8,
"""

        self.report.read(aligned_reads)
        self.report.write_nuc_header(self.report_file)
        self.report.write_nuc_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testPartialCodonNucleotideReport(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R1-seed,15,0,9,0,AAATT
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,genome.pos,\
A,C,G,T,N,del,ins,clip,v3_overlap,coverage,exact_coverage
R1-seed,R1,15,1,1,1,9,0,0,0,0,0,0,0,0,9,
R1-seed,R1,15,2,2,2,9,0,0,0,0,0,0,0,0,9,9
R1-seed,R1,15,3,3,3,9,0,0,0,0,0,0,0,0,9,18
R1-seed,R1,15,4,4,4,0,0,0,9,0,0,0,0,0,9,18
R1-seed,R1,15,5,5,5,0,0,0,9,0,0,0,0,0,9,9
"""

        self.report.read(aligned_reads)
        self.report.write_nuc_header(self.report_file)
        self.report.write_nuc_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testPartialStartCodonNucleotideReport(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R1-seed,15,0,9,0,TTAGG
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,genome.pos,\
A,C,G,T,N,del,ins,clip,v3_overlap,coverage,exact_coverage
R1-seed,R1,15,1,5,5,0,0,0,9,0,0,0,0,0,9,
R1-seed,R1,15,2,6,6,0,0,0,9,0,0,0,0,0,9,
R1-seed,R1,15,3,7,7,9,0,0,0,0,0,0,0,0,9,
R1-seed,R1,15,4,8,8,0,0,9,0,0,0,0,0,0,9,
R1-seed,R1,15,5,9,9,0,0,9,0,0,0,0,0,0,9,
"""

        self.report.read(aligned_reads)
        self.report.write_nuc_header(self.report_file)
        self.report.write_nuc_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testReadPairGapInMiddleOfAminoReport(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R3-seed,15,0,9,0,AAATTTCAGACCCCAnnnnnnnnnTACTAC
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R3-seed,R3,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R3-seed,R3,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R3-seed,R3,15,7,3,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R3-seed,R3,15,10,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,9
R3-seed,R3,15,13,5,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
"""

        self.report.read(aligned_reads)
        self.report.write_amino_header(self.report_file)
        self.report.write_amino_counts()

        self.assertEqual(expected_text, self.report_file.getvalue())

    def testLowQualityNucleotideReport(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R1-seed,15,0,9,0,AAATNT
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,genome.pos,\
A,C,G,T,N,del,ins,clip,v3_overlap,coverage,exact_coverage
R1-seed,R1,15,1,1,1,9,0,0,0,0,0,0,0,0,9,
R1-seed,R1,15,2,2,2,9,0,0,0,0,0,0,0,0,9,
R1-seed,R1,15,3,3,3,9,0,0,0,0,0,0,0,0,9,
R1-seed,R1,15,4,4,4,0,0,0,9,0,0,0,0,0,9,
R1-seed,R1,15,5,5,5,0,0,0,0,9,0,0,0,0,0,
R1-seed,R1,15,6,6,6,0,0,0,9,0,0,0,0,0,9,
"""

        self.report.read(aligned_reads)
        self.report.write_nuc_header(self.report_file)
        self.report.write_nuc_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testLowQualityAminoReport(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R1-seed,15,0,9,0,AAATNT
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,4,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0
"""

        self.report.read(aligned_reads)
        self.report.write_amino_header(self.report_file)
        self.report.write_amino_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testPartialDeletionAminoReport(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R1-seed,15,0,9,0,AAAT-T
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,4,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0
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
        aligned_reads = prepare_reads("""\
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

        self.assertEqual(expected_text, self.report_file.getvalue())

    def testShiftedReadingFrameNucleotideReport(self):
        """ The seed's reading frame doesn't match the coordinate reference's
        reading frame, so there is an extra nucleotide at the beginning of the
        reads.
        It will try padding the first codon to see which of the three possible
        reading frames gives the highest alignment score.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R1-seed,15,0,9,0,GAAATTTCGA
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,genome.pos,\
A,C,G,T,N,del,ins,clip,v3_overlap,coverage,exact_coverage
R1-seed,R1,15,2,1,1,9,0,0,0,0,0,0,0,0,9,
R1-seed,R1,15,3,2,2,9,0,0,0,0,0,0,0,0,9,
R1-seed,R1,15,4,3,3,9,0,0,0,0,0,0,0,0,9,
R1-seed,R1,15,5,4,4,0,0,0,9,0,0,0,0,0,9,
R1-seed,R1,15,6,5,5,0,0,0,9,0,0,0,0,0,9,
R1-seed,R1,15,7,6,6,0,0,0,9,0,0,0,0,0,9,
R1-seed,R1,15,8,7,7,0,9,0,0,0,0,0,0,0,9,
R1-seed,R1,15,9,8,8,0,0,9,0,0,0,0,0,0,9,
R1-seed,R1,15,10,9,9,9,0,0,0,0,0,0,0,0,9,
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
        aligned_reads = prepare_reads("""\
R1-seed,15,0,9,0,AAA---AGG
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,genome.pos,\
A,C,G,T,N,del,ins,clip,v3_overlap,coverage,exact_coverage
R1-seed,R1,15,1,1,1,9,0,0,0,0,0,0,0,0,9,
R1-seed,R1,15,2,2,2,9,0,0,0,0,0,0,0,0,9,
R1-seed,R1,15,3,3,3,9,0,0,0,0,0,0,0,0,9,
R1-seed,R1,15,4,4,4,0,0,0,0,0,9,0,0,0,9,
R1-seed,R1,15,5,5,5,0,0,0,0,0,9,0,0,0,9,
R1-seed,R1,15,6,6,6,0,0,0,0,0,9,0,0,0,9,
R1-seed,R1,15,7,7,7,9,0,0,0,0,0,0,0,0,9,
R1-seed,R1,15,8,8,8,0,0,9,0,0,0,0,0,0,9,
R1-seed,R1,15,9,9,9,0,0,9,0,0,0,0,0,0,9,
"""

        self.report.read(aligned_reads)
        self.report.write_nuc_header(self.report_file)
        self.report.write_nuc_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testDeletionBetweenSeedAndCoordinateNucleotideReport(self):
        """ Coordinate sequence is KFQTPREH, and this aligned read is KFQPREH.

        Must be a deletion in the seed reference with respect to the coordinate
        reference.
        """
        self.report.landmarks = yaml.safe_load("""\
- seed_pattern: R3-.*
  coordinates: R3-seed
  landmarks:
    - {name: R3, start: 1, end: 27}
    """)
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R3-seed,15,0,9,0,AAATTTCAGCCACGAGAGCAT
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,genome.pos,\
A,C,G,T,N,del,ins,clip,v3_overlap,coverage,exact_coverage
R3-seed,R3,15,1,1,1,9,0,0,0,0,0,0,0,0,9,
R3-seed,R3,15,2,2,2,9,0,0,0,0,0,0,0,0,9,
R3-seed,R3,15,3,3,3,9,0,0,0,0,0,0,0,0,9,
R3-seed,R3,15,4,4,4,0,0,0,9,0,0,0,0,0,9,
R3-seed,R3,15,5,5,5,0,0,0,9,0,0,0,0,0,9,
R3-seed,R3,15,6,6,6,0,0,0,9,0,0,0,0,0,9,
R3-seed,R3,15,7,7,7,0,9,0,0,0,0,0,0,0,9,
R3-seed,R3,15,8,8,8,9,0,0,0,0,0,0,0,0,9,
R3-seed,R3,15,9,9,9,0,0,9,0,0,0,0,0,0,9,
R3-seed,R3,15,,10,10,0,0,0,0,0,9,0,0,0,9,
R3-seed,R3,15,,11,11,0,0,0,0,0,9,0,0,0,9,
R3-seed,R3,15,,12,12,0,0,0,0,0,9,0,0,0,9,
R3-seed,R3,15,10,13,13,0,9,0,0,0,0,0,0,0,9,
R3-seed,R3,15,11,14,14,0,9,0,0,0,0,0,0,0,9,
R3-seed,R3,15,12,15,15,9,0,0,0,0,0,0,0,0,9,
R3-seed,R3,15,13,16,16,0,9,0,0,0,0,0,0,0,9,
R3-seed,R3,15,14,17,17,0,0,9,0,0,0,0,0,0,9,
R3-seed,R3,15,15,18,18,9,0,0,0,0,0,0,0,0,9,
R3-seed,R3,15,16,19,19,0,0,9,0,0,0,0,0,0,9,
R3-seed,R3,15,17,20,20,9,0,0,0,0,0,0,0,0,9,
R3-seed,R3,15,18,21,21,0,0,9,0,0,0,0,0,0,9,
R3-seed,R3,15,19,22,22,0,9,0,0,0,0,0,0,0,9,
R3-seed,R3,15,20,23,23,9,0,0,0,0,0,0,0,0,9,
R3-seed,R3,15,21,24,24,0,0,0,9,0,0,0,0,0,9,
"""

        self.report.read(aligned_reads)
        self.report.write_nuc_header(self.report_file)
        self.report.write_nuc_counts()

        self.assertEqual(expected_text, self.report_file.getvalue())

    def testDeletionBetweenSeedAndCoordinateAminoReport(self):
        """ Coordinate sequence is KFQTPREH, and this aligned read is KFQPREH.

        Must be a deletion in the seed reference with respect to the coordinate
        reference.
        """
        self.report.landmarks = yaml.safe_load("""\
- seed_pattern: R3-.*
  coordinates: R3-seed
  landmarks:
    - {name: R3, start: 1, end: 27}
    """)
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R3-seed,15,0,9,0,AAATTTCAGCCACGAGAGCAT
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R3-seed,R3,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R3-seed,R3,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R3-seed,R3,15,7,3,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R3-seed,R3,15,,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,9
R3-seed,R3,15,10,5,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R3-seed,R3,15,13,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,9
R3-seed,R3,15,16,7,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R3-seed,R3,15,19,8,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
"""

        self.report.read(aligned_reads)
        self.report.write_amino_header(self.report_file)
        self.report.write_amino_counts()

        self.assertEqual(expected_text, self.report_file.getvalue())

    def testDeletionBetweenSeedAndConsensusAminoReport(self):
        """ Coordinate and consensus are KFGPR, but seed is KFPR.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R5-seed,15,0,9,0,AAATTTGGCCCCCGACCTCAGGTCACTCTTTGG
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
R5-seed,R5,15,16,6,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R5-seed,R5,15,19,7,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R5-seed,R5,15,22,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,9
R5-seed,R5,15,25,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,9
R5-seed,R5,15,28,10,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
"""

        self.report.write_amino_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_amino_counts()

        self.assertEqual(expected_text, self.report_file.getvalue())

    def testDeletionWithMinorityVariant(self):
        """ Aligned reads are mostly K-R, but some are KFR.

        Must be a deletion in the sample with respect to the seed reference,
        but some variants in the sample do not have that deletion.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R1-seed,15,0,5,0,AAA---AGG
R1-seed,15,0,2,0,AAATTTAGG
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

        self.assertEqual(expected_text, self.report_file.getvalue())

    def testDeletionNotAlignedToCodons(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R1-seed,15,0,5,0,AAAA---GG
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

        self.assertEqual(expected_text, self.report_file.getvalue())

    def testInsertionBetweenSeedAndCoordinateAminoReport(self):
        """ Coordinate sequence is KFQTPREH, and this aligned read is HERKFQTGPREH.

        The G must be an insertion in the seed reference with respect to the
        coordinate reference.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R3-seed,15,0,9,0,CATGAGCGAAAATTTCAGACTGGGCCCCGAGAGCAT
""")
        self.report.landmarks = yaml.safe_load("""\
- seed_pattern: R3-.*
  coordinates: R3-seed
  landmarks:
    # Extra 3 positions for stop codons to get dropped, one codon overlaps.
    - {name: R3, start: 1, end: 27, colour: lightblue}
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
seed,mixture_cutoff,region,ref_region_pos,ref_genome_pos,query_pos,insertion
R3-seed,MAX,R3,12,12,21,GGG
R3-seed,0.100,R3,12,12,21,GGG
"""

        self.report.read(aligned_reads)
        self.report.write_insertions()
        self.report.write_nuc_header(StringIO())
        self.report.write_nuc_counts()  # calculates insertion counts
        self.report.write_amino_header(self.report_file)
        self.report.write_amino_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())
        self.assertEqual(expected_insertions,
                         self.report.insert_writer.insert_file.getvalue())

    def testInsertionBetweenSeedAndCoordinateNucleotideReport(self):
        """ Coordinate sequence is KFQTPREH, and this aligned read is HERKFQTGPREH.

        The G must be an insertion in the seed reference with respect to the
        coordinate reference.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R3-seed,15,0,9,0,CATGAGCGAAAATTTCAGACTGGGCCCCGAGAGCAT
""")
        self.report.landmarks = yaml.safe_load("""\
- seed_pattern: R3-.*
  coordinates: R3-seed
  landmarks:
    # Extra 3 positions for stop codons to get dropped, one codon overlaps.
    - {name: R3, start: 1, end: 27, colour: lightblue}
""")
        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,genome.pos,\
A,C,G,T,N,del,ins,clip,v3_overlap,coverage,exact_coverage
R3-seed,R3,15,10,1,1,9,0,0,0,0,0,0,0,0,9,
R3-seed,R3,15,11,2,2,9,0,0,0,0,0,0,0,0,9,
R3-seed,R3,15,12,3,3,9,0,0,0,0,0,0,0,0,9,
R3-seed,R3,15,13,4,4,0,0,0,9,0,0,0,0,0,9,
R3-seed,R3,15,14,5,5,0,0,0,9,0,0,0,0,0,9,
R3-seed,R3,15,15,6,6,0,0,0,9,0,0,0,0,0,9,
R3-seed,R3,15,16,7,7,0,9,0,0,0,0,0,0,0,9,
R3-seed,R3,15,17,8,8,9,0,0,0,0,0,0,0,0,9,
R3-seed,R3,15,18,9,9,0,0,9,0,0,0,0,0,0,9,
R3-seed,R3,15,19,10,10,9,0,0,0,0,0,0,0,0,9,
R3-seed,R3,15,20,11,11,0,9,0,0,0,0,0,0,0,9,
R3-seed,R3,15,21,12,12,0,0,0,9,0,0,9,0,0,9,
R3-seed,R3,15,25,13,13,0,9,0,0,0,0,0,0,0,9,
R3-seed,R3,15,26,14,14,0,9,0,0,0,0,0,0,0,9,
R3-seed,R3,15,27,15,15,0,9,0,0,0,0,0,0,0,9,
R3-seed,R3,15,28,16,16,0,9,0,0,0,0,0,0,0,9,
R3-seed,R3,15,29,17,17,0,0,9,0,0,0,0,0,0,9,
R3-seed,R3,15,30,18,18,9,0,0,0,0,0,0,0,0,9,
R3-seed,R3,15,31,19,19,0,0,9,0,0,0,0,0,0,9,
R3-seed,R3,15,32,20,20,9,0,0,0,0,0,0,0,0,9,
R3-seed,R3,15,33,21,21,0,0,9,0,0,0,0,0,0,9,
R3-seed,R3,15,34,22,22,0,9,0,0,0,0,0,0,0,9,
R3-seed,R3,15,35,23,23,9,0,0,0,0,0,0,0,0,9,
R3-seed,R3,15,36,24,24,0,0,0,9,0,0,0,0,0,9,
"""

        self.report.read(aligned_reads)
        self.report.write_insertions()
        self.report.write_nuc_header(self.report_file)
        self.report.write_nuc_counts()

        self.assertEqual(expected_text, self.report_file.getvalue())

    def testInsertionsSortedByCount(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R3-seed,15,0,9,0,CATGAGCGAAAATTTCAGACTGGGCCCCGAGAGCAT
R3-seed,15,0,8,0,CATGAGCGAAAATTTCAGACTAAACCCCGAGAGCAT
""")
        self.report.landmarks = yaml.safe_load("""\
- seed_pattern: R3-.*
  coordinates: R3-seed
  landmarks:
    # Extra 3 positions for stop codons to get dropped, one codon overlaps.
    - {name: R3, start: 1, end: 27, colour: lightblue}
""")
        expected_insertions = """\
seed,mixture_cutoff,region,ref_region_pos,ref_genome_pos,query_pos,insertion
R3-seed,MAX,R3,12,12,21,GGG
R3-seed,0.100,R3,12,12,21,RRR
"""

        self.report.read(aligned_reads)
        self.report.write_insertions()

        self.assertEqual(expected_insertions,
                         self.report.insert_writer.insert_file.getvalue())

    def testInsertionsSortedByLeft(self):
        """ Two insertions within a single consensus, sorted by position.

        Consensus is HERKFQTGPRKEHQFKL
        Reference is  ERKF-TGPRK-HQFKL (without the dashes)
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
        "CATGAGCGAAAATTTACTGGGCCCCGAAAACATCAGTTTAAACTC"
      ]
    },
    "R3": {
      "is_nucleotide": false,
      "reference": [
        "ERKFTGPRKHQFKL"
      ]
    }
  }
}
"""))
        self.report.landmarks = yaml.safe_load("""\
- seed_pattern: R3-seed
  coordinates: R3-seed
  landmarks:
    # Extra 3 nucleotides at end, because stop codons will get dropped.
    - {name: R3, start: 4, end: 48, frame: 0}
""")
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R3-seed,15,0,9,0,CATGAGCGAAAATTTCAGACTGGGCCCCGAAAAGAGCATCAGTTTAAACTC
""")
        expected_insertions = """\
seed,mixture_cutoff,region,ref_region_pos,ref_genome_pos,query_pos,insertion
R3-seed,MAX,R3,12,15,15,CAG
R3-seed,0.100,R3,12,15,15,CAG
R3-seed,MAX,R3,27,30,33,GAG
R3-seed,0.100,R3,27,30,33,GAG
"""

        self.report.read(aligned_reads)
        self.report.write_insertions()

        self.assertEqual(expected_insertions,
                         self.report.insert_writer.insert_file.getvalue())

    def testInsertionInDifferentReadingFrame(self):
        """ Delete part of the first codon to throw off the reading frame.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R3-seed,15,0,9,0,AATTTCAGACTGGGCCCCGAGAGCAT
""")
        self.report.landmarks = yaml.safe_load("""\
- seed_pattern: R3-.*
  coordinates: R3-seed
  landmarks:
    # Extra 3 positions for stop codons to get dropped, one codon overlaps.
    - {name: R3, start: 1, end: 27, colour: lightblue}
""")

        expected_insertions = """\
seed,mixture_cutoff,region,ref_region_pos,ref_genome_pos,query_pos,insertion
R3-seed,MAX,R3,12,12,11,GGG
R3-seed,0.100,R3,12,12,11,GGG
"""

        self.report.read(aligned_reads)
        self.report.write_amino_header(StringIO())
        self.report.write_amino_counts()
        self.report.write_insertions()

        self.assertEqual(expected_insertions,
                         self.report.insert_writer.insert_file.getvalue())

    def testInsertionInSomeReads(self):
        """ Not all reads have the insertion, some end before it.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R3-seed,15,0,9,0,AAATTTCAGACTGGGCCCCGAGAGCAT
R3-seed,15,1,5,0,AAATTTCAG
R3-seed,15,2,4,0,AAATTTCAGACTG
""")
        self.report.landmarks = yaml.safe_load("""\
- seed_pattern: R3-.*
  coordinates: R3-seed
  landmarks:
    # Extra 3 positions for stop codons to get dropped, one codon overlaps.
    - {name: R3, start: 1, end: 27, colour: lightblue}
""")

        expected_insertions = """\
seed,mixture_cutoff,region,ref_region_pos,ref_genome_pos,query_pos,insertion
R3-seed,MAX,R3,12,12,12,GGG
R3-seed,0.100,R3,12,12,12,Ggg
"""

        self.report.read(aligned_reads)
        self.report.write_amino_header(StringIO())
        self.report.write_amino_counts()
        self.report.write_insertions()

        self.assertEqual(expected_insertions,
                         self.report.insert_writer.insert_file.getvalue())

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
        self.report.landmarks = yaml.safe_load("""\
- seed_pattern: R3-seed
  coordinates: R3-seed
  landmarks:
    # Extra 3 nucleotides at end, because stop codons will get dropped.
    - {name: R3a, start: 1, end: 27, frame: 0}
    - {name: R3b, start: 1, end: 27, frame: 0}
""")
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R3-seed,15,0,9,0,AAATTTCAGACTGGGCCCCGAGAGCAT
""")

        expected_insertions = """\
seed,mixture_cutoff,region,ref_region_pos,ref_genome_pos,query_pos,insertion
R3-seed,MAX,R3a,12,12,12,GGG
R3-seed,0.100,R3a,12,12,12,GGG
"""

        self.report.read(aligned_reads)
        self.report.write_insertions()

        self.assertEqual(expected_insertions,
                         self.report.insert_writer.insert_file.getvalue())

    def testInsertionsRelativeToConsensus(self):
        """ Test that insertions relative to the consensus are handled correctly """
        aligned_reads = prepare_reads("""\
R1-seed,15,0,10,0,AAATTTAGG
""")
        conseq_ins_csv = StringIO("""\
qname,fwd_rev,refname,pos,insert,qual
Example_read_1,F,R1-seed,3,AAC,AAA
Example_read_2,F,R1-seed,3,AAC,AAA
""")

        expected_insertions = ("""\
seed,mixture_cutoff,region,ref_region_pos,ref_genome_pos,query_pos,insertion
R1-seed,0.100,R1,3,3,3,aac
""")

        self.report.read_insertions(conseq_ins_csv)
        self.report.write_amino_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_nuc_header(StringIO())
        self.report.write_nuc_counts()  # calculates ins counts
        self.report.write_amino_counts()
        self.report.insert_writer.write(self.report.inserts,
                                        self.report.detail_seed,
                                        self.report.reports,
                                        self.report.report_nucleotides,
                                        self.report.landmarks,
                                        self.report.consensus_builder)
        self.assertEqual(expected_insertions,
                         self.report.insert_writer.insert_file.getvalue())

    def testGapBetweenForwardAndReverse(self):
        """ Lower-case n represents a gap between forward and reverse reads.

        Region R2 has sequence KFGPR, so this read has a gap at the end of G
        and beginning of P. Partial codons at the ends of a read or next to the
        gap are ignored, even though G is still unambiguous.
        """
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R2-seed,15,0,5,0,AAATTTGGnnCCCGA
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R2-seed,R2,15,1,1,0,0,0,0,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5
R2-seed,R2,15,4,2,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5
R2-seed,R2,15,13,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,5
"""

        self.report.read(aligned_reads)
        self.report.write_amino_header(self.report_file)
        self.report.write_amino_counts()

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testFailedAlignmentAminoReport(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R1-seed,15,0,2,0,TTATCCTAC
""")

        expected_text = ""

        self.report.read(aligned_reads)
        self.report.write_amino_counts(self.report_file)

        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

    def testFailedAlignmentFailureReport(self):
        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
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
        aligned_reads = prepare_reads("""\
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
        aligned_reads = prepare_reads("""\
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
        aligned_reads = prepare_reads("""\
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
        aligned_reads = prepare_reads("""\
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
        aligned_reads = prepare_reads("""\
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
        aligned_reads = prepare_reads("""\
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
        aligned_reads = prepare_reads("""\
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
        self.report.landmarks = yaml.safe_load("""\
- seed_pattern: R1-seed
  coordinates: R1-seed
  landmarks:
    # Extra 3 nucleotides at end, because stop codons will get dropped.
    - {name: R1a, start: 4, end: 15, frame: 0}
    - {name: R1b, start: 1, end: 15, frame: 0}
""")

        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R1-seed,15,0,9,0,AAATTT
""")

        expected_text = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap,coverage
R1-seed,R1a,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1a,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1b,15,1,2,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1b,15,4,3,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
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
        self.report.landmarks = yaml.safe_load("""\
- seed_pattern: R1-seed
  coordinates: R1-seed
  landmarks:
    # Extra 3 nucleotides at end, because stop codons will get dropped.
    - {name: R1a, start: 4, end: 15, frame: 0}
    - {name: R1b, start: 1, end: 15, frame: 0}
""")

        # refname,qcut,rank,count,offset,seq
        aligned_reads = prepare_reads("""\
R1-seed,15,0,9,0,AAATTT
R1-seed,15,0,1,0,CCCGGG
""")

        expected_text = """\
seed,region,q-cutoff,consensus-percent-cutoff,offset,sequence
R1-seed,R1a,15,MAX,0,AAATTT
R1-seed,R1a,15,0.100,0,MMMKKK
R1-seed,R1b,15,MAX,3,AAATTT
R1-seed,R1b,15,0.100,3,MMMKKK
"""

        self.report.read(aligned_reads)
        self.report.write_consensus_regions_header(self.report_file)
        self.report.combine_reports()
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
        aligned_reads = prepare_reads("R1-seed,15,0,10,0,AAATTTAGG")
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
        aligned_reads = prepare_reads("R1-seed,15,0,10,0,AAA---AGG")
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
        aligned_reads = prepare_reads("R1-seed,15,0,10,0,AAAA---GG")
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
        aligned_reads = prepare_reads("R1-seed,15,0,10,0,AA---AAGG")
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
        aligned_reads = prepare_reads("R2-seed,15,0,10,0,AA------ACCGAGA")
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
        aligned_reads = prepare_reads("R1-seed,15,0,10,1,AA---AGG")
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
        aligned_reads = prepare_reads("R1-seed,15,0,10,2,AAA---AGG")
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
        aligned_reads = prepare_reads("R-seed,15,0,10,1,AAA---CAGTTTTTTTTC---AGCAT")
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
        aligned_reads = prepare_reads("R2-seed,15,0,10,0,AA-TCG--CCCGAGA")
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
        aligned_reads = prepare_reads("R2-seed,15,0,10,0,AA--TC----CGAGA")
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
        aligned_reads = prepare_reads("R3-seed,15,0,10,0,AA--TTCAGACCCC-CGAGAGCAT")
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
        aligned_reads = prepare_reads("R3-seed,15,0,10,0,AA--TTCAGACCCCA-GAGAGCAT")
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
        aligned_reads = prepare_reads("R3-seed,15,0,10,0,AA--TTCAGACCCC-CGA---CAT")
        self.report.remap_conseqs = {'R3-seed': 'AAATTTCAGACCCCACGAGAGCAT'}
        expected_reads = [dict(refname='R3-seed',
                               qcut='15',
                               rank='0',
                               count='10',
                               offset='0',
                               seq='AA--TTCAGACCCC-CGA---CAT')]

        reads = list(self.report.align_deletions(aligned_reads))

        self.assertEqual(expected_reads, reads)

    def testCombinedCoordinateConcordance(self):
        aligned_reads = prepare_reads("R1A-seed,15,0,10,0,AAATTTAGGTAG")
        expected_file = """\
reference,region,pct_concordance,pct_covered
R1A-seed,R1A,100.0,100.0
R1A-seed,R1A_second,0.0,0.0
"""

        self.report.read(aligned_reads)
        concordance_file = StringIO()
        self.report.write_concordance_header(concordance_file)
        self.report.write_nuc_header(StringIO())
        self.report.write_nuc_counts()
        self.report.combine_reports()
        self.report.write_coordinate_concordance(self.report.concordance_writer, use_combined_reports=True)
        self.assertMultiLineEqual(expected_file, concordance_file.getvalue())

    def testReportCoordinateConcordance(self):
        aligned_reads = prepare_reads("R1A-seed,15,0,10,0,AAATTTAGGTAG")
        expected_file = """\
reference,region,pct_concordance,pct_covered
R1A-seed,R1A,100.0,100.0
R1A-seed,R1A_second,0.0,0.0
"""

        self.report.read(aligned_reads)
        concordance_file = StringIO()
        self.report.write_concordance_header(concordance_file)
        self.report.write_nuc_header(StringIO())
        self.report.write_nuc_counts()
        self.report.write_coordinate_concordance(self.report.concordance_writer)
        self.assertMultiLineEqual(expected_file, concordance_file.getvalue())

    def testCoordinateConcordanceCoverage(self):
        aligned_reads = prepare_reads("R1A-seed,15,0,10,0,AAATTT")
        expected_file = """\
reference,region,pct_concordance,pct_covered
R1A-seed,R1A,100.0,50.0
R1A-seed,R1A_second,0.0,0.0
"""

        self.report.read(aligned_reads)
        concordance_file = StringIO()
        self.report.write_concordance_header(concordance_file)
        self.report.write_nuc_header(StringIO())
        self.report.write_nuc_counts()
        self.report.write_coordinate_concordance(self.report.concordance_writer)
        self.assertMultiLineEqual(expected_file, concordance_file.getvalue())

    def testCoordinateConcordanceMismatch(self):
        aligned_reads = prepare_reads("R1A-seed,15,0,10,0,AAATTTGGGTAG")
        # 1 different nuc here:                                 ^
        expected_file = """\
reference,region,pct_concordance,pct_covered
R1A-seed,R1A,91.66666666666667,100.0
R1A-seed,R1A_second,0.0,0.0
"""

        self.report.read(aligned_reads)
        concordance_file = StringIO()
        self.report.write_concordance_header(concordance_file)
        self.report.write_nuc_header(StringIO())
        self.report.write_nuc_counts()
        self.report.write_coordinate_concordance(self.report.concordance_writer)
        self.assertMultiLineEqual(expected_file, concordance_file.getvalue())

    def testCoordinateConcordanceMismatchCoverage(self):
        aligned_reads = prepare_reads("R1A-seed,15,0,10,0,AAATTTGGG")
        # 1 different nuc here:                                 ^
        expected_file = """\
reference,region,pct_concordance,pct_covered
R1A-seed,R1A,88.88888888888889,75.0
R1A-seed,R1A_second,0.0,0.0
"""

        self.report.read(aligned_reads)
        concordance_file = StringIO()
        self.report.write_concordance_header(concordance_file)
        self.report.write_nuc_header(StringIO())
        self.report.write_nuc_counts()
        self.report.write_coordinate_concordance(self.report.concordance_writer)
        self.assertMultiLineEqual(expected_file, concordance_file.getvalue())

    def testCoordinateConcordanceDeletion(self):
        aligned_reads = prepare_reads("R1A-seed,15,0,10,0,AAATTT---TAG")
        expected_file = """\
reference,region,pct_concordance,pct_covered
R1A-seed,R1A,100.0,75.0
R1A-seed,R1A_second,0.0,0.0
"""

        self.report.read(aligned_reads)
        concordance_file = StringIO()
        self.report.write_concordance_header(concordance_file)
        self.report.write_nuc_header(StringIO())
        self.report.write_nuc_counts()
        self.report.write_coordinate_concordance(self.report.concordance_writer)
        self.assertMultiLineEqual(expected_file, concordance_file.getvalue())

    def testDetailedCombinedCoordinateConcordance(self):
        aligned_reads = prepare_reads("R1A-seed,15,0,10,12,CCGAGACCTCAGGTCACTCTTTGGTAG")
        expected_file = """\
reference,region,pct_concordance,pct_covered,position
R1A-seed,R1A_second,100.0,100.0,10
R1A-seed,R1A_second,100.0,100.0,11
R1A-seed,R1A_second,100.0,100.0,12
R1A-seed,R1A_second,100.0,100.0,13
R1A-seed,R1A_second,100.0,100.0,14
R1A-seed,R1A_second,100.0,100.0,15
R1A-seed,R1A_second,100.0,100.0,16
R1A-seed,R1A_second,100.0,100.0,17
"""

        self.report.read(aligned_reads)
        concordance_file = StringIO()
        self.report.write_concordance_header(StringIO())
        self.report.write_concordance_header(concordance_file, is_detailed=True)
        self.report.write_nuc_header(StringIO())
        self.report.write_nuc_counts()
        self.report.combine_reports()
        self.report.write_coordinate_concordance(self.report.concordance_writer,
                                                 self.report.detailed_concordance_writer,
                                                 use_combined_reports=True)
        self.assertMultiLineEqual(expected_file, concordance_file.getvalue())

    def testDetailedCoordinateConcordance(self):
        aligned_reads = prepare_reads("R1A-seed,15,0,10,12,CCGAGACCTCAGGTCACTCTTTGGTAG")
        expected_file = """\
reference,region,pct_concordance,pct_covered,position
R1A-seed,R1A_second,100.0,100.0,10
R1A-seed,R1A_second,100.0,100.0,11
R1A-seed,R1A_second,100.0,100.0,12
R1A-seed,R1A_second,100.0,100.0,13
R1A-seed,R1A_second,100.0,100.0,14
R1A-seed,R1A_second,100.0,100.0,15
R1A-seed,R1A_second,100.0,100.0,16
R1A-seed,R1A_second,100.0,100.0,17
"""

        self.report.read(aligned_reads)
        concordance_file = StringIO()
        self.report.write_concordance_header(StringIO())
        self.report.write_concordance_header(concordance_file, is_detailed=True)
        self.report.write_nuc_header(StringIO())
        self.report.write_nuc_counts()
        self.report.write_coordinate_concordance(self.report.concordance_writer,
                                                 self.report.detailed_concordance_writer)
        self.assertMultiLineEqual(expected_file, concordance_file.getvalue())

    def testDetailedCoordinateConcordanceCoverage(self):
        aligned_reads = prepare_reads("R1A-seed,15,0,10,12,CCGAGACCTCAGGTCACTCTTTGG")
        expected_file = """\
reference,region,pct_concordance,pct_covered,position
R1A-seed,R1A_second,100.0,100.0,10
R1A-seed,R1A_second,100.0,100.0,11
R1A-seed,R1A_second,100.0,100.0,12
R1A-seed,R1A_second,100.0,100.0,13
R1A-seed,R1A_second,100.0,100.0,14
R1A-seed,R1A_second,95.0,95.0,15
R1A-seed,R1A_second,90.0,90.0,16
R1A-seed,R1A_second,85.0,85.0,17
"""

        self.report.read(aligned_reads)
        concordance_file = StringIO()
        self.report.write_concordance_header(StringIO())
        self.report.write_concordance_header(concordance_file, is_detailed=True)
        self.report.write_nuc_header(StringIO())
        self.report.write_nuc_counts()
        self.report.write_coordinate_concordance(self.report.concordance_writer,
                                                 self.report.detailed_concordance_writer)
        self.assertMultiLineEqual(expected_file, concordance_file.getvalue())

    def testDetailedCoordinateConcordanceMismatch(self):
        aligned_reads = prepare_reads("R1A-seed,15,0,10,12,CCGAGCCCTCTGGTCACTCTGTGGTAG")
        # mismatch:                                             ^    ^         ^
        expected_file = """\
reference,region,pct_concordance,pct_covered,position
R1A-seed,R1A_second,90.0,100.0,10
R1A-seed,R1A_second,85.0,100.0,11
R1A-seed,R1A_second,85.0,100.0,12
R1A-seed,R1A_second,85.0,100.0,13
R1A-seed,R1A_second,85.0,100.0,14
R1A-seed,R1A_second,85.0,100.0,15
R1A-seed,R1A_second,90.0,100.0,16
R1A-seed,R1A_second,90.0,100.0,17
"""

        self.report.read(aligned_reads)
        concordance_file = StringIO()
        self.report.write_concordance_header(StringIO())
        self.report.write_concordance_header(concordance_file, is_detailed=True)
        self.report.write_nuc_header(StringIO())
        self.report.write_nuc_counts()
        self.report.write_coordinate_concordance(self.report.concordance_writer,
                                                 self.report.detailed_concordance_writer)
        self.assertMultiLineEqual(expected_file, concordance_file.getvalue())

    def testDetailedCoordinateConcordanceDeletion(self):
        aligned_reads = prepare_reads("R1A-seed,15,0,10,12,CCGAGACCTCAGGTCACTCTT---TAG")
        expected_file = """\
reference,region,pct_concordance,pct_covered,position
R1A-seed,R1A_second,100.0,100.0,10
R1A-seed,R1A_second,100.0,100.0,11
R1A-seed,R1A_second,95.0,95.0,12
R1A-seed,R1A_second,90.0,90.0,13
R1A-seed,R1A_second,85.0,85.0,14
R1A-seed,R1A_second,85.0,85.0,15
R1A-seed,R1A_second,85.0,85.0,16
R1A-seed,R1A_second,85.0,85.0,17
"""

        self.report.read(aligned_reads)
        concordance_file = StringIO()
        self.report.write_concordance_header(StringIO())
        self.report.write_concordance_header(concordance_file, is_detailed=True)
        self.report.write_nuc_header(StringIO())
        self.report.write_nuc_counts()
        self.report.write_coordinate_concordance(self.report.concordance_writer,
                                                 self.report.detailed_concordance_writer)
        self.assertMultiLineEqual(expected_file, concordance_file.getvalue())


class InsertionWriterTest(unittest.TestCase):
    def setUp(self):
        self.insert_file = StringIO()
        self.writer = InsertionWriter(self.insert_file)
        self.writer.start_group(seed='R1-seed', qcut=15)
        self.nuc_seq_acdef = 'GCTTGTGACGAGTTT'
        self.nuc_seq_afdef = 'GCTTTTGACGAGTTT'

    def testNoInserts(self):
        expected_text = """\
seed,mixture_cutoff,region,ref_region_pos,ref_genome_pos,query_pos,insertion
"""

        self.writer.add_nuc_read(offset_sequence=self.nuc_seq_acdef, count=1)
        self.writer.write(insertions={},
                          seed_name='',
                          report_aminos_all=[],
                          report_nucleotides_all=[],
                          landmarks=None,
                          consensus_builder=ConsensusBuilder([0.1, 'MAX'], 0))

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())

    def testInsert(self):
        expected_text = """\
seed,mixture_cutoff,region,ref_region_pos,ref_genome_pos,query_pos,insertion
R1-seed,MAX,R1,6,6,6,GAC
R1-seed,0.100,R1,6,6,6,GAC
"""
        expected_counts = {('R1-seed', 'R1'): {6: 1}}

        report_aminos = {'R1': [ReportAmino(SeedAmino(0), 1),
                                ReportAmino(SeedAmino(3), 2),
                                ReportAmino(SeedAmino(9), 3),
                                ReportAmino(SeedAmino(12), 4)]}

        self.writer.add_nuc_read(offset_sequence=self.nuc_seq_acdef, count=1)
        self.writer.write(insertions={'R1': [6]},
                          seed_name='R1-seed',
                          report_aminos_all=report_aminos,
                          report_nucleotides_all={'R1': []},
                          landmarks=yaml.safe_load(LANDMARKS_YAML),
                          consensus_builder=ConsensusBuilder(['MAX', 0.1], 0))

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())
        self.assertEqual(expected_counts, self.writer.insert_pos_counts)

    def testInsertNuc(self):
        expected_text = """\
seed,mixture_cutoff,region,ref_region_pos,ref_genome_pos,query_pos,insertion
R1-seed,MAX,R1,2,2,2,T
R1-seed,0.100,R1,2,2,2,T
"""

        report_nucleotides = {'R1': [ReportNucleotide(1, SeedNucleotide()),
                                     ReportNucleotide(2, SeedNucleotide()),
                                     ReportNucleotide(3, SeedNucleotide()),
                                     ReportNucleotide(4, SeedNucleotide())]}
        report_nucleotides['R1'][0].seed_nucleotide.consensus_index = 0
        report_nucleotides['R1'][1].seed_nucleotide.consensus_index = 1
        report_nucleotides['R1'][2].seed_nucleotide.consensus_index = 3
        report_nucleotides['R1'][3].seed_nucleotide.consensus_index = 4

        self.writer.add_nuc_read(offset_sequence=self.nuc_seq_acdef, count=1)
        self.writer.write(insertions={'R1': [2]},
                          seed_name='R1-seed',
                          report_aminos_all={'R1': []},
                          report_nucleotides_all=report_nucleotides,
                          landmarks=yaml.safe_load(LANDMARKS_YAML),
                          consensus_builder=ConsensusBuilder(['MAX', 0.1], 0))

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())

    def testInsertDifferentReadingFrame(self):
        """ Add a partial codon at the start of the read to shift the reading
        frame.
        """
        expected_text = """\
seed,mixture_cutoff,region,ref_region_pos,ref_genome_pos,query_pos,insertion
R1-seed,MAX,R1,6,6,7,GAC
R1-seed,0.100,R1,6,6,7,GAC
"""
        report_aminos = {'R1': [ReportAmino(SeedAmino(1), 1),
                                ReportAmino(SeedAmino(4), 2),
                                ReportAmino(SeedAmino(10), 3),
                                ReportAmino(SeedAmino(13), 4)]}

        self.writer.add_nuc_read(offset_sequence='A' + self.nuc_seq_acdef,
                                 count=1)
        self.writer.write(insertions={'R1': [7]},
                          seed_name='R1-seed',
                          report_aminos_all=report_aminos,
                          report_nucleotides_all={'R1': []},
                          landmarks=yaml.safe_load(LANDMARKS_YAML),
                          consensus_builder=ConsensusBuilder(['MAX', 0.1], 0))

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())

    def testInsertWithOffset(self):
        expected_text = """\
seed,mixture_cutoff,region,ref_region_pos,ref_genome_pos,query_pos,insertion
R1-seed,MAX,R1,6,6,6,GAC
R1-seed,0.100,R1,6,6,6,GAC
"""
        report_aminos = {'R1': [ReportAmino(SeedAmino(0), 1),
                                ReportAmino(SeedAmino(3), 2),
                                ReportAmino(SeedAmino(9), 3),
                                ReportAmino(SeedAmino(12), 4)]}

        #                                            C  D  E  F
        self.writer.add_nuc_read(offset_sequence='---TGTGACGAGTTT', count=1)
        self.writer.write(insertions={'R1': [6]},
                          seed_name='R1-seed',
                          report_aminos_all=report_aminos,
                          report_nucleotides_all={'R1': []},
                          landmarks=yaml.safe_load(LANDMARKS_YAML),
                          consensus_builder=ConsensusBuilder(['MAX', 0.1], 0))

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())

    def testInsertWithDeletion(self):
        expected_text = """\
seed,mixture_cutoff,region,ref_region_pos,ref_genome_pos,query_pos,insertion
"""

        report_aminos = {'R1': [ReportAmino(SeedAmino(0), 1),
                                ReportAmino(SeedAmino(3), 2),
                                ReportAmino(SeedAmino(9), 3),
                                ReportAmino(SeedAmino(12), 4)]}

        #                                         C  D     E  F
        self.writer.add_nuc_read(offset_sequence='TGTGAC---GAGTTT', count=1)
        self.writer.write(insertions={'R1': [3, 6]},
                          seed_name='R1-seed',
                          report_aminos_all=report_aminos,
                          report_nucleotides_all={'R1': []},
                          landmarks=yaml.safe_load(LANDMARKS_YAML),
                          consensus_builder=ConsensusBuilder(['MAX', 0.1], 0))

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())

    def testTwoInsertsWithOffset(self):
        expected_text = """\
seed,mixture_cutoff,region,ref_region_pos,ref_genome_pos,query_pos,insertion
R1-seed,MAX,R1,6,6,6,GAC
R1-seed,0.100,R1,6,6,6,GAC
R1-seed,MAX,R1,9,9,12,TTT
R1-seed,0.100,R1,9,9,12,TTT
"""
        report_aminos = {'R1': [ReportAmino(SeedAmino(0), 1),
                                ReportAmino(SeedAmino(3), 2),
                                ReportAmino(SeedAmino(9), 3),
                                ReportAmino(SeedAmino(15), 4)]}

        #                                            C  D  E  F  G
        self.writer.add_nuc_read(offset_sequence='---TGTGACGAGTTTGGG', count=1)
        self.writer.write(insertions={'R1': [6, 12]},
                          seed_name='R1-seed',
                          report_aminos_all=report_aminos,
                          report_nucleotides_all={'R1': []},
                          landmarks=yaml.safe_load(LANDMARKS_YAML),
                          consensus_builder=ConsensusBuilder(['MAX', 0.1], 0))

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())

    def testInsertsWithVariants(self):
        expected_text = """\
seed,mixture_cutoff,region,ref_region_pos,ref_genome_pos,query_pos,insertion
R1-seed,MAX,R1,6,6,6,GAC
R1-seed,0.100,R1,6,6,6,GAC
"""
        report_aminos = {'R1': [ReportAmino(SeedAmino(0), 1),
                                ReportAmino(SeedAmino(3), 2),
                                ReportAmino(SeedAmino(9), 3),
                                ReportAmino(SeedAmino(12), 4)]}

        self.writer.add_nuc_read(offset_sequence=self.nuc_seq_acdef, count=1)
        self.writer.add_nuc_read(offset_sequence=self.nuc_seq_afdef, count=1)
        self.writer.write(insertions={'R1': [6]},
                          seed_name='R1-seed',
                          report_aminos_all=report_aminos,
                          report_nucleotides_all={'R1': []},
                          landmarks=yaml.safe_load(LANDMARKS_YAML),
                          consensus_builder=ConsensusBuilder(['MAX', 0.1], 0))

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())

    def testDifferentInserts(self):
        expected_text = """\
seed,mixture_cutoff,region,ref_region_pos,ref_genome_pos,query_pos,insertion
R1-seed,MAX,R1,3,3,3,TTT
R1-seed,0.100,R1,3,3,3,TKT
"""
        report_aminos = {'R1': [ReportAmino(SeedAmino(0), 1),
                                ReportAmino(SeedAmino(6), 2),
                                ReportAmino(SeedAmino(9), 3),
                                ReportAmino(SeedAmino(12), 4)]}

        self.writer.add_nuc_read(offset_sequence=self.nuc_seq_acdef, count=2)
        self.writer.add_nuc_read(offset_sequence=self.nuc_seq_afdef, count=3)
        self.writer.write(insertions={'R1': [3]},
                          seed_name='R1-seed',
                          report_aminos_all=report_aminos,
                          report_nucleotides_all={'R1': []},
                          landmarks=yaml.safe_load(LANDMARKS_YAML),
                          consensus_builder=ConsensusBuilder(['MAX', 0.1], 0))

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())

    def testMulticharacterInsert(self):
        expected_text = """\
seed,mixture_cutoff,region,ref_region_pos,ref_genome_pos,query_pos,insertion
R1-seed,MAX,R1,6,6,6,GACGAG
R1-seed,0.100,R1,6,6,6,GACGAG
"""
        report_aminos = {'R1': [ReportAmino(SeedAmino(0), 1),
                                ReportAmino(SeedAmino(3), 2),
                                ReportAmino(SeedAmino(12), 3)]}

        self.writer.add_nuc_read(offset_sequence=self.nuc_seq_acdef, count=1)
        self.writer.write(insertions={'R1': [6, 9]},
                          seed_name='R1-seed',
                          report_aminos_all=report_aminos,
                          report_nucleotides_all={'R1': []},
                          landmarks=yaml.safe_load(LANDMARKS_YAML),
                          consensus_builder=ConsensusBuilder(['MAX', 0.1], 0))

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())

    def testReadGapInInsert(self):
        nuc_seq = 'GCTCTnGACGAGTTT'

        expected_text = """\
seed,mixture_cutoff,region,ref_region_pos,ref_genome_pos,query_pos,insertion
"""
        report_aminos = {'R1': [ReportAmino(SeedAmino(0), 1),
                                ReportAmino(SeedAmino(6), 2),
                                ReportAmino(SeedAmino(9), 3)]}

        self.writer.add_nuc_read(nuc_seq, count=1)
        self.writer.write(insertions={'R1': [3]},
                          seed_name='R1-seed',
                          report_aminos_all=report_aminos,
                          report_nucleotides_all={'R1': []},
                          landmarks=yaml.safe_load(LANDMARKS_YAML),
                          consensus_builder=ConsensusBuilder(['MAX', 0.1], 0))

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())

    def testUnsortedInserts(self):
        expected_text = """\
seed,mixture_cutoff,region,ref_region_pos,ref_genome_pos,query_pos,insertion
R1-seed,MAX,R1,6,6,6,GACGAG
R1-seed,0.100,R1,6,6,6,GACGAG
"""
        report_aminos = {'R1': [ReportAmino(SeedAmino(0), 1),
                                ReportAmino(SeedAmino(3), 2),
                                ReportAmino(SeedAmino(12), 3)]}

        self.writer.add_nuc_read(offset_sequence=self.nuc_seq_acdef, count=1)
        self.writer.write(insertions={'R1': [9, 6]},
                          seed_name='R1-seed',
                          report_aminos_all=report_aminos,
                          report_nucleotides_all={'R1': []},
                          landmarks=yaml.safe_load(LANDMARKS_YAML),
                          consensus_builder=ConsensusBuilder(['MAX', 0.1], 0))

        self.assertMultiLineEqual(expected_text, self.insert_file.getvalue())
