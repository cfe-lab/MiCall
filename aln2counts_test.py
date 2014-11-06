import StringIO
import unittest

import aln2counts
import project_config

class StubbedSequenceReport(aln2counts.SequenceReport):
    def __init__(self, 
                 insert_writer, 
                 projects, 
                 conseq_mixture_cutoffs):
        aln2counts.SequenceReport.__init__(self,
                                           insert_writer,
                                           projects,
                                           conseq_mixture_cutoffs)
        self.overrides = {}
        
    def _pair_align(self, reference, query):
        override = self.overrides.get((reference, query))
        return (override
                if override is not None
                else aln2counts.SequenceReport._pair_align(self, reference, query))
    
    def add_override(self, reference, query, aligned_query, aligned_reference):
        self.overrides[(reference, query)] = (aligned_query, aligned_reference)

class SequenceReportTest(unittest.TestCase):
    def setUp(self):
        self.insertion_file = StringIO.StringIO()
        insert_writer = aln2counts.InsertionWriter(
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
          "seed_region": "R1-seed"
        }
      ]
    },
    "R2": {
      "max_variants": 10,
      "regions": [
        {
          "coordinate_region": "R2",
          "seed_region": "R2-seed"
        }
      ]
    },
    "R3": {
      "max_variants": 10,
      "regions": [
        {
          "coordinate_region": "R3",
          "seed_region": "R3-seed"
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
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = """\
E1234_S1,R1-seed,15,0,9,0,AAATTT
""".splitlines(True)
        expected_text = """\
sample,region,q-cutoff,s-number,consensus-percent-cutoff,sequence
E1234,R1-seed,15,S1,MAX,AAATTT
E1234,R1-seed,15,S1,0.100,AAATTT
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
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = """\
E1234_S1,R1-seed,15,0,9,0,AAATTT
E1234_S1,R1-seed,15,0,1,0,CCCGGG
""".splitlines(True)
        expected_consensus = "KF"
         
        self.report.read(aligned_reads)
        consensus = self.report.consensus
         
        self.assertSequenceEqual(expected_consensus, consensus)
     
    def testSingleReadAminoReport(self):
        """ In this sample, there is a single read with two codons.
        AAA -> K
        TTT -> F
        The coordinate reference has three codons, so the third position is
        empty.
        """
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = """\
E1234_S1,R1-seed,15,0,9,0,AAATTT
""".splitlines(True)
         
        expected_text = """\
sample,seed,region,q-cutoff,query.aa.pos,refseq.aa.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
E1234_S1,R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0
E1234_S1,R1-seed,R1,15,2,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
E1234_S1,R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
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
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = """\
E1234_S1,R1-seed,15,0,9,0,AAATTT
""".splitlines(True)
          
        expected_text = """\
sample,seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T
E1234_S1,R1-seed,R1,15,1,1,9,0,0,0
E1234_S1,R1-seed,R1,15,2,2,9,0,0,0
E1234_S1,R1-seed,R1,15,3,3,9,0,0,0
E1234_S1,R1-seed,R1,15,4,4,0,0,0,9
E1234_S1,R1-seed,R1,15,5,5,0,0,0,9
E1234_S1,R1-seed,R1,15,6,6,0,0,0,9
E1234_S1,R1-seed,R1,15,,7,0,0,0,0
E1234_S1,R1-seed,R1,15,,8,0,0,0,0
E1234_S1,R1-seed,R1,15,,9,0,0,0,0
"""
        
        self.report.write_nuc_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_nuc_counts(self.report_file)
        
        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())
     
    def testOffsetNucleotideReport(self):
        """ The first row provides alignment so the partial codon at the start
        of the second row will map to the reference.
        """
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = """\
E1234_S1,R1-seed,15,0,1,3,TTT
E1234_S1,R1-seed,15,0,8,5,TCGA
""".splitlines(True)
          
        #sample,seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T
        expected_text = """\
E1234_S1,R1-seed,R1,15,1,1,0,0,0,0
E1234_S1,R1-seed,R1,15,2,2,0,0,0,0
E1234_S1,R1-seed,R1,15,3,3,0,0,0,0
E1234_S1,R1-seed,R1,15,4,4,0,0,0,1
E1234_S1,R1-seed,R1,15,5,5,0,0,0,1
E1234_S1,R1-seed,R1,15,6,6,0,0,0,9
E1234_S1,R1-seed,R1,15,7,7,0,8,0,0
E1234_S1,R1-seed,R1,15,8,8,0,0,8,0
E1234_S1,R1-seed,R1,15,9,9,8,0,0,0
"""
        
        self.report.read(aligned_reads)
        self.report.write_nuc_counts(self.report_file)
        
        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())
    
    def testPartialCodonNucleotideReport(self):
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = """\
E1234_S1,R1-seed,15,0,9,0,AAATT
""".splitlines(True)
        
        #sample,seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T
        expected_text = """\
E1234_S1,R1-seed,R1,15,1,1,9,0,0,0
E1234_S1,R1-seed,R1,15,2,2,9,0,0,0
E1234_S1,R1-seed,R1,15,3,3,9,0,0,0
E1234_S1,R1-seed,R1,15,4,4,0,0,0,9
E1234_S1,R1-seed,R1,15,5,5,0,0,0,9
E1234_S1,R1-seed,R1,15,6,6,0,0,0,0
E1234_S1,R1-seed,R1,15,,7,0,0,0,0
E1234_S1,R1-seed,R1,15,,8,0,0,0,0
E1234_S1,R1-seed,R1,15,,9,0,0,0,0
"""
          
        self.report.read(aligned_reads)
        self.report.write_nuc_counts(self.report_file)
         
        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())
    
    def testDeletionBetweenSeedAndCoordinateNucleotideReport(self):
        """ Coordinate sequence is KFGPR, and this aligned read is KFPR.
         
        Must be a deletion in the seed reference with respect to the coordinate
        reference.
        """
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = """\
E1234_S1,R2-seed,15,0,9,0,AAATTTCCCCGA
""".splitlines(True)
           
        #sample,seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T
        expected_text = """\
E1234_S1,R2-seed,R2,15,1,1,9,0,0,0
E1234_S1,R2-seed,R2,15,2,2,9,0,0,0
E1234_S1,R2-seed,R2,15,3,3,9,0,0,0
E1234_S1,R2-seed,R2,15,4,4,0,0,0,9
E1234_S1,R2-seed,R2,15,5,5,0,0,0,9
E1234_S1,R2-seed,R2,15,6,6,0,0,0,9
E1234_S1,R2-seed,R2,15,,7,0,0,0,0
E1234_S1,R2-seed,R2,15,,8,0,0,0,0
E1234_S1,R2-seed,R2,15,,9,0,0,0,0
E1234_S1,R2-seed,R2,15,7,10,0,9,0,0
E1234_S1,R2-seed,R2,15,8,11,0,9,0,0
E1234_S1,R2-seed,R2,15,9,12,0,9,0,0
E1234_S1,R2-seed,R2,15,10,13,0,9,0,0
E1234_S1,R2-seed,R2,15,11,14,0,0,9,0
E1234_S1,R2-seed,R2,15,12,15,9,0,0,0
"""
           
        self.report.read(aligned_reads)
        self.report.write_nuc_counts(self.report_file)
        
        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())
    
    def testDeletionBetweenSeedAndCoordinateAminoReport(self):
        """ Coordinate sequence is KFGPR, and this aligned read is KFPR.
         
        Must be a deletion in the seed reference with respect to the coordinate
        reference.
        """
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = """\
E1234_S1,R2-seed,15,0,9,0,AAATTTCCCCGA
""".splitlines(True)
           
        #sample,seed,region,q-cutoff,query.aa.pos,refseq.aa.pos,
        #                  A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
        expected_text = """\
E1234_S1,R2-seed,R2,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0
E1234_S1,R2-seed,R2,15,2,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
E1234_S1,R2-seed,R2,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
E1234_S1,R2-seed,R2,15,3,4,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0
E1234_S1,R2-seed,R2,15,4,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0
"""
        
        self.report.read(aligned_reads)
        self.report.write_amino_counts(self.report_file)
        
        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())
    
    def testDeletionWithMinorityVariant(self):
        """ Aligned reads are mostly K-R, but some are KFR.
        
        Must be a deletion in the sample with respect to the seed reference,
        but some variants in the sample do not have that deletion.
        """
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = """\
E1234_S1,R1-seed,15,0,5,0,AAA---CGA
E1234_S1,R1-seed,15,0,2,0,AAATTTCGA
""".splitlines(True)
           
        #sample,seed,region,q-cutoff,query.aa.pos,refseq.aa.pos,
        #                  A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
        expected_text = """\
E1234_S1,R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,7,0,0,0,0,0,0,0,0,0,0,0,0
E1234_S1,R1-seed,R1,15,2,2,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
E1234_S1,R1-seed,R1,15,3,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,0,0,0,0,0,0
"""
        
        self.report.read(aligned_reads)
        self.report.write_amino_counts(self.report_file)
        
        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())
    
    def testInsertionBetweenSeedAndCoordinateAminoReport(self):
        """ Coordinate sequence is KFQTPREH, and this aligned read is KFQTGPREH.
         
        The G must be an insertion in the seed reference with respect to the
        coordinate reference.
        """
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = """\
E1234_S1,R3-seed,15,0,9,0,AAATTTCAGACTGGGCCCCGAGAGCAT
""".splitlines(True)
           
        #sample,seed,region,q-cutoff,query.aa.pos,refseq.aa.pos,
        #                  A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
        expected_text = """\
E1234_S1,R3-seed,R3,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0
E1234_S1,R3-seed,R3,15,2,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
E1234_S1,R3-seed,R3,15,3,3,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0
E1234_S1,R3-seed,R3,15,4,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0
E1234_S1,R3-seed,R3,15,6,5,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0
E1234_S1,R3-seed,R3,15,7,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0
E1234_S1,R3-seed,R3,15,8,7,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
E1234_S1,R3-seed,R3,15,9,8,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0
"""
        expected_insertions = """\
sample,seed,region,qcut,left,insert,count
E1234_S1,R3-seed,R3,15,5,G,9
"""

        self.report.read(aligned_reads)
        self.report.write_amino_counts(self.report_file)
        self.report.write_insertions()
        
        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())
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
          "seed_region": "R3-seed"
        },
        {
          "coordinate_region": "R3b",
          "seed_region": "R3-seed"
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
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = """\
E1234_S1,R3-seed,15,0,9,0,AAATTTCAGACTGGGCCCCGAGAGCAT
""".splitlines(True)
        
        expected_insertions = """\
sample,seed,region,qcut,left,insert,count
E1234_S1,R3-seed,R3a,15,5,G,9
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
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = """\
E1234_S1,R2-seed,15,0,5,0,AAATTTGGnnCCCGA
""".splitlines(True)
           
        #sample,seed,region,q-cutoff,query.aa.pos,refseq.aa.pos,
        #                  A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
        expected_text = """\
E1234_S1,R2-seed,R2,15,1,1,0,0,0,0,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0
E1234_S1,R2-seed,R2,15,2,2,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
E1234_S1,R2-seed,R2,15,3,3,0,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
E1234_S1,R2-seed,R2,15,4,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
E1234_S1,R2-seed,R2,15,5,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,0,0,0,0,0,0
"""
        
        self.report.read(aligned_reads)
        self.report.write_amino_counts(self.report_file)
        
        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())
    
    def testFailedAlignmentAminoReport(self):
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = """\
E1234_S1,R1-seed,15,0,2,0,TTATCCTAC
""".splitlines(True)
           
        expected_text = ""
        
        self.report.read(aligned_reads)
        self.report.write_amino_counts(self.report_file)
        
        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())
    
    def testFailedAlignmentFailureReport(self):
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = """\
E1234_S1,R1-seed,15,0,2,0,TTATCCTAC
""".splitlines(True)
        
        expected_text = """\
sample,seed,region,qcut,queryseq,refseq
E1234_S1,R1-seed,R1,15,LSY,KFR
"""
        
        self.report.write_failure_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_failure(self.report_file)
        
        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())
    
    def testFailedAlignmentWithOffset(self):
        """ Be careful that an offset from the seed reference doesn't match
        the dashes in the failed alignment.
        """
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = """\
E1234_S1,R1-seed,15,0,2,3,TTATCCTAC
""".splitlines(True)
        
        expected_text = """\
sample,seed,region,qcut,queryseq,refseq
E1234_S1,R1-seed,R1,15,-LSY,KFR
"""
        
        self.report.write_failure_header(self.report_file)
        self.report.read(aligned_reads)
        self.report.write_failure(self.report_file)
        
        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())
       
    def testNoFailureReport(self):
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = """\
E1234_S1,R1-seed,15,0,9,0,AAATTT
""".splitlines(True)
         
        expected_text = ""
        
        self.report.read(aligned_reads)
        self.report.write_failure(self.report_file)
        
        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())
     
    def testRegionWithoutCoordinateReferenceNucleotideReport(self):
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = """\
E1234_S1,R-NO-COORD,15,0,9,0,AAATTT
""".splitlines(True)
          
        #sample,seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T
        expected_text = """\
E1234_S1,R-NO-COORD,R-NO-COORD,15,1,,9,0,0,0
E1234_S1,R-NO-COORD,R-NO-COORD,15,2,,9,0,0,0
E1234_S1,R-NO-COORD,R-NO-COORD,15,3,,9,0,0,0
E1234_S1,R-NO-COORD,R-NO-COORD,15,4,,0,0,0,9
E1234_S1,R-NO-COORD,R-NO-COORD,15,5,,0,0,0,9
E1234_S1,R-NO-COORD,R-NO-COORD,15,6,,0,0,0,9
E1234_S1,R-NO-COORD,R-NO-COORD,15,7,,0,0,0,0
E1234_S1,R-NO-COORD,R-NO-COORD,15,8,,0,0,0,0
E1234_S1,R-NO-COORD,R-NO-COORD,15,9,,0,0,0,0
"""
        
        self.report.read(aligned_reads)
        self.report.write_nuc_counts(self.report_file)
        
        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())
     
    def testRegionWithoutCoordinateReferenceFailureReport(self):
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = """\
E1234_S1,R-NO-COORD,15,0,9,0,AAATTT
""".splitlines(True)
          
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
          "seed_region": "R1-seed"
        },
        {
          "coordinate_region": "R1b",
          "seed_region": "R1-seed"
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
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = """\
E1234_S1,R1-seed,15,0,9,0,AAATTT
""".splitlines(True)
        
        expected_text = """\
E1234_S1,R1-seed,R1a,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0
E1234_S1,R1-seed,R1a,15,2,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
E1234_S1,R1-seed,R1a,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
E1234_S1,R1-seed,R1b,15,,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
E1234_S1,R1-seed,R1b,15,1,2,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0
E1234_S1,R1-seed,R1b,15,2,3,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
E1234_S1,R1-seed,R1b,15,,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
"""        
        
        self.report.read(aligned_reads)
        self.report.write_amino_counts(self.report_file)
        
        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())
    
    def testVariantReport(self):
        """ Reference is KFR, variants are KFRx10, GFRx8, KFGx6, KFSx5
        """
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = """\
E1234_S1,R1-seed,15,0,10,0,AAATTTCGA
E1234_S1,R1-seed,15,1,8,0,GGGTTTCGA
E1234_S1,R1-seed,15,2,6,0,AAATTTGGG
E1234_S1,R1-seed,15,3,5,0,AAATTTTCA
""".splitlines(True)
        
        expected_text = """\
sample,seed,qcut,region,index,count,seq
E1234_S1,R1-seed,15,R1,0,10,AAATTTCGA
E1234_S1,R1-seed,15,R1,1,8,GGGTTTCGA
E1234_S1,R1-seed,15,R1,2,6,AAATTTGGG
E1234_S1,R1-seed,15,R1,3,5,AAATTTTCA
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
          "seed_region": "R1-seed"
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
        
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = """\
E1234_S1,R1-seed,15,0,10,0,AAATTTCGA
E1234_S1,R1-seed,15,1,8,0,GGGTTTCGA
E1234_S1,R1-seed,15,2,6,0,AAATTTGGG
E1234_S1,R1-seed,15,3,5,0,AAATTTTCA
""".splitlines(True)
        
        # sample,seed,qcut,region,index,count,seq
        expected_text = """\
E1234_S1,R1-seed,15,R1,0,10,AAATTTCGA
E1234_S1,R1-seed,15,R1,1,8,GGGTTTCGA
E1234_S1,R1-seed,15,R1,2,6,AAATTTGGG
"""
        self.report.read(aligned_reads)
        self.report.write_nuc_variants(self.report_file)
        
        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())
    
    def testVariantClippingEnd(self):
        """ Coordinate reference is KFR, so last codons should be clipped
        """
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = """\
E1234_S1,R1-seed,15,0,10,0,AAATTTCGATCA
E1234_S1,R1-seed,15,1,8,0,GGGTTTCGATCA
E1234_S1,R1-seed,15,2,6,0,AAATTTGGGTCA
E1234_S1,R1-seed,15,3,5,0,AAATTTTCA
""".splitlines(True)
        
        # sample,seed,qcut,region,index,count,seq
        expected_text = """\
E1234_S1,R1-seed,15,R1,0,10,AAATTTCGA
E1234_S1,R1-seed,15,R1,1,8,GGGTTTCGA
E1234_S1,R1-seed,15,R1,2,6,AAATTTGGG
E1234_S1,R1-seed,15,R1,3,5,AAATTTTCA
"""
        
        self.report.read(aligned_reads)
        self.report.write_nuc_variants(self.report_file)
        
        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())
    
    def testVariantClippingStart(self):
        """ Coordinate reference is KFR, so first codons should be clipped
        """
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = """\
E1234_S1,R1-seed,15,0,10,0,TCAAAATTTCGA
E1234_S1,R1-seed,15,1,8,0,TCAGGGTTTCGA
E1234_S1,R1-seed,15,2,6,0,TCAAAATTTGGG
E1234_S1,R1-seed,15,3,5,0,TCAAAATTTTCA
""".splitlines(True)
        
        # sample,seed,qcut,region,index,count,seq
        expected_text = """\
E1234_S1,R1-seed,15,R1,0,10,AAATTTCGA
E1234_S1,R1-seed,15,R1,1,8,GGGTTTCGA
E1234_S1,R1-seed,15,R1,2,6,AAATTTGGG
E1234_S1,R1-seed,15,R1,3,5,AAATTTTCA
"""
        
        self.report.read(aligned_reads)
        self.report.write_nuc_variants(self.report_file)
        
        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())
    
    def testVariantWhenAlignFails(self):
        """ Coordinate reference is KFR, read is SGG.
        """
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = """\
E1234_S1,R1-seed,15,0,10,0,TCAGGGGGG
""".splitlines(True)
          
        expected_text = ""
        
        self.report.read(aligned_reads)
        self.report.write_nuc_variants(self.report_file)
        
        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())
    
    def testVariantShortRead(self):
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = """\
E1234_S1,R1-seed,15,0,10,0,AAATTT
E1234_S1,R1-seed,15,1,8,0,GGGTTT
""".splitlines(True)
        
        # sample,seed,qcut,region,index,count,seq
        expected_text = """\
E1234_S1,R1-seed,15,R1,0,10,AAATTT
E1234_S1,R1-seed,15,R1,1,8,GGGTTT
"""
        
        self.report.read(aligned_reads)
        self.report.write_nuc_variants(self.report_file)
        
        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())
    
    def testVariantOffset(self):
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = """\
E1234_S1,R1-seed,15,0,10,0,AAATTT
E1234_S1,R1-seed,15,1,8,3,TTT
""".splitlines(True)
          
        # sample,seed,qcut,region,index,count,seq
        expected_text = """\
E1234_S1,R1-seed,15,R1,0,10,AAATTT
E1234_S1,R1-seed,15,R1,1,8,---TTT
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
          "seed_region": "R1-seed"
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
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = """\
E1234_S1,R1-seed,15,0,10,0,TCAAAATTTCGA
E1234_S1,R1-seed,15,1,9,0,TCAGGGTTTCGA
E1234_S1,R1-seed,15,2,7,0,TCAGGGTTTGGG
E1234_S1,R1-seed,15,3,5,0,CGAGGGTTTCGA
""".splitlines(True)
          
        # sample,seed,qcut,region,index,count,seq
        expected_text = """\
E1234_S1,R1-seed,15,R1,0,14,GGGTTTCGA
E1234_S1,R1-seed,15,R1,1,10,AAATTTCGA
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
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = """\
E1234_S1,R1-seed,15,0,10,0,AAATTTGGGCGA
""".splitlines(True)
          
        # sample,seed,qcut,region,index,count,seq
        expected_text = """\
E1234_S1,R1-seed,15,R1,0,10,AAATTTGGGCGA
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
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = """\
E1234_S1,R1-seed,15,0,10,0,AAATTTGGG
""".splitlines(True)
          
        # sample,seed,qcut,region,index,count,seq
        expected_text = """\
E1234_S1,R1-seed,15,R1,0,10,AAATTT
"""
        
        self.report.read(aligned_reads)
        self.report.write_nuc_variants(self.report_file)
        
        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())
    
    def testVariantCanProcessReadsTwice(self):
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = iter("""\
E1234_S1,R1-seed,15,0,10,0,AAATTT
E1234_S1,R1-seed,15,1,8,3,TTT
""".splitlines(True))
          
        # sample,seed,qcut,region,index,count,seq
        expected_text = """\
E1234_S1,R1-seed,15,R1,0,10,AAATTT
E1234_S1,R1-seed,15,R1,1,8,---TTT
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
          "seed_region": "R1-seed"
        },
        {
          "coordinate_region": "R1b",
          "seed_region": "R1-seed"
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
        #sample,refname,qcut,rank,count,offset,seq
        aligned_reads = """\
E1234_S1,R1-seed,15,0,9,0,AAATTT
""".splitlines(True)
        
        # sample,seed,qcut,region,index,count,seq
        expected_text = """\
E1234_S1,R1-seed,15,R1a,0,9,AAATTT
E1234_S1,R1-seed,15,R1b,0,9,AAATTT
"""        
        
        self.report.read(aligned_reads)
        self.report.write_nuc_variants(self.report_file)
        
        self.assertMultiLineEqual(expected_text, self.report_file.getvalue())

class InsertionWriterTest(unittest.TestCase):
    def setUp(self):
        self.writer = aln2counts.InsertionWriter(insert_file=StringIO.StringIO())
        self.writer.start_group(sample_name='E1234_S1', seed='R1-seed', qcut=15)
        
    def testNoInserts(self):
        expected_text = """\
sample,seed,region,qcut,left,insert,count
"""
        
        self.writer.add_read(offset_sequence='ACDEF', count=1)
        self.writer.write(inserts=[], region='R1')
        
        self.assertMultiLineEqual(expected_text, self.writer.insert_file.getvalue())
        
    def testInsert(self):
        expected_text = """\
sample,seed,region,qcut,left,insert,count
E1234_S1,R1-seed,R1,15,3,D,1
"""
        
        self.writer.add_read(offset_sequence='ACDEF', count=1)
        self.writer.write(inserts=[2], region='R1')
        
        self.assertMultiLineEqual(expected_text, self.writer.insert_file.getvalue())
        
    def testInsertWithOffset(self):
        expected_text = """\
sample,seed,region,qcut,left,insert,count
E1234_S1,R1-seed,R1,15,3,D,1
"""
        
        self.writer.add_read(offset_sequence='-CDEFG', count=1)
        self.writer.write(inserts=[2], region='R1')
        
        self.assertMultiLineEqual(expected_text, self.writer.insert_file.getvalue())
        
    def testTwoInsertsWithOffset(self):
        expected_text = """\
sample,seed,region,qcut,left,insert,count
E1234_S1,R1-seed,R1,15,3,D,1
E1234_S1,R1-seed,R1,15,5,F,1
"""
        
        self.writer.add_read(offset_sequence='-CDEFG', count=1)
        self.writer.write(inserts=[2, 4], region='R1')
        
        self.assertMultiLineEqual(expected_text, self.writer.insert_file.getvalue())

    def testInsertsWithVariants(self):
        expected_text = """\
sample,seed,region,qcut,left,insert,count
E1234_S1,R1-seed,R1,15,3,D,2
"""
        
        self.writer.add_read(offset_sequence='ACDEF', count=1)
        self.writer.add_read(offset_sequence='AFDEF', count=1)
        self.writer.write(inserts=[2], region='R1')
        
        self.assertMultiLineEqual(expected_text, self.writer.insert_file.getvalue())

    def testDifferentInserts(self):
        expected_text = """\
sample,seed,region,qcut,left,insert,count
E1234_S1,R1-seed,R1,15,3,D,2
E1234_S1,R1-seed,R1,15,3,F,3
"""
        
        self.writer.add_read(offset_sequence='ACDEF', count=2)
        self.writer.add_read(offset_sequence='ACFEF', count=3)
        self.writer.write(inserts=[2], region='R1')
        
        self.assertMultiLineEqual(expected_text, self.writer.insert_file.getvalue())

    def testMulticharacterInsert(self):
        expected_text = """\
sample,seed,region,qcut,left,insert,count
E1234_S1,R1-seed,R1,15,3,DE,1
"""
        
        self.writer.add_read(offset_sequence='ACDEF', count=1)
        self.writer.write(inserts=[2,3], region='R1')
        
        self.assertMultiLineEqual(expected_text, self.writer.insert_file.getvalue())

    def testUnsortedInserts(self):
        expected_text = """\
sample,seed,region,qcut,left,insert,count
E1234_S1,R1-seed,R1,15,3,DE,1
"""
        
        self.writer.add_read(offset_sequence='ACDEF', count=1)
        self.writer.write(inserts=(3, 2), region='R1')
        
        self.assertMultiLineEqual(expected_text, self.writer.insert_file.getvalue())

class TranslateTest(unittest.TestCase):
    def testSingleCodon(self):
        nucs = 'TTT'
        expected_aminos = 'F'

        aminos = aln2counts.translate(nucs)
        
        self.assertEqual(expected_aminos, aminos)
        
    def testPartialCodon(self):
        nucs = 'TTTC'
        expected_aminos = 'F'

        aminos = aln2counts.translate(nucs)
        
        self.assertEqual(expected_aminos, aminos)
        
    def testTwoCodons(self):
        nucs = 'TTTCCT'
        expected_aminos = 'FP'

        aminos = aln2counts.translate(nucs)
        
        self.assertEqual(expected_aminos, aminos)
        
    def testOffset(self):
        nucs = "TTTCCT"
        offset = 3
        expected_aminos = "-FP"
        
        aminos = aln2counts.translate(nucs, offset)
        
        self.assertEqual(expected_aminos, aminos)
        
    def testSingleDashAmbiguous(self):
        nucs = '-TT'
        expected_aminos = '?'

        aminos = aln2counts.translate(nucs)
        
        self.assertEqual(expected_aminos, aminos)
        
    def testSingleDashUnambiguous(self):
        nucs = 'CG-' # CGA, CGC, CGG, CGT all map to R
        expected_aminos = 'R'

        aminos = aln2counts.translate(nucs)
        
        self.assertEqual(expected_aminos, aminos)
        
    def testTwoDashes(self):
        nucs = '--T'
        expected_aminos = '?'

        aminos = aln2counts.translate(nucs)
        
        self.assertEqual(expected_aminos, aminos)
        
    def testThreeDashes(self):
        nucs = '---'
        expected_aminos = '-'

        aminos = aln2counts.translate(nucs)
        
        self.assertEqual(expected_aminos, aminos)
        
    def testAmbiguousBasesThatAreSynonyms(self):
        nucs = 'TTY' # TTC or TTT: both map to F
        expected_aminos = 'F'

        aminos = aln2counts.translate(nucs)
        
        self.assertEqual(expected_aminos, aminos)
        
    def testTwoAmbiguousBasesThatAreSynonyms(self):
        nucs = 'MGR' # CGA, CGG, AGA, or AGG: all map to R
        expected_aminos = '?'

        aminos = aln2counts.translate(nucs)
        
        self.assertEqual(expected_aminos, aminos)

class SeedAminoTest(unittest.TestCase):
    def setUp(self):
        self.amino = aln2counts.SeedAmino(None)
    
    def testSingleRead(self):
        """ Read a single codon, and report on counts.
        Columns are:       A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
        """
        nuc_seq = 'AAA' # -> K
        expected_counts = '0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0'
        
        self.amino.count_nucleotides(nuc_seq, 8)
        counts = self.amino.get_report()
        
        self.assertSequenceEqual(expected_counts, counts)
    
    def testDifferentCodon(self):
        """ Read two different codons, and report on counts.
        Columns are:       A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
        """
        nuc_seq1 = 'AAA' # -> K
        nuc_seq2 = 'GGG' # -> G
        expected_counts = '0,0,0,0,0,5,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0'
        
        self.amino.count_nucleotides(nuc_seq1, 8)
        self.amino.count_nucleotides(nuc_seq2, 5)
        counts = self.amino.get_report()
        
        self.assertSequenceEqual(expected_counts, counts)
    
    def testSameAminoAcid(self):
        """ Read same codon twice, and report on counts.
        Columns are:       A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
        """
        nuc_seq1 = 'AAA' # -> K
        nuc_seq2 = 'AAG' # -> K
        expected_counts = '0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0'
        
        self.amino.count_nucleotides(nuc_seq1, 4)
        self.amino.count_nucleotides(nuc_seq2, 5)
        counts = self.amino.get_report()
        
        self.assertSequenceEqual(expected_counts, counts)
    
    def testNucleotides(self):
        nuc_seq1 = 'AAA' # -> K
        nuc_seq2 = 'AAG' # -> K
        expected_nuc_counts = '4,0,5,0'
        
        self.amino.count_nucleotides(nuc_seq1, 4)
        self.amino.count_nucleotides(nuc_seq2, 5)
        counts = self.amino.nucleotides[2].get_report()
        
        self.assertSequenceEqual(expected_nuc_counts, counts)
    
    def testConsensus(self):
        nuc_seq1 = 'AAA' # -> K
        nuc_seq2 = 'GGG' # -> G
        expected_consensus = 'G'
        
        self.amino.count_nucleotides(nuc_seq1, 4)
        self.amino.count_nucleotides(nuc_seq2, 5)
        consensus = self.amino.get_consensus()
        
        self.assertSequenceEqual(expected_consensus, consensus)
    
    def testConsensusMixture(self):
        nuc_seq1 = 'AAA' # -> K
        nuc_seq2 = 'GGG' # -> G
        nuc_seq3 = 'TTT' # -> F
        allowed_consensus_values = ('G', 'K')
        
        self.amino.count_nucleotides(nuc_seq1, 4)
        self.amino.count_nucleotides(nuc_seq2, 4)
        self.amino.count_nucleotides(nuc_seq3, 3)
        consensus = self.amino.get_consensus()
        
        self.assertIn(consensus, allowed_consensus_values)
    
    def testConsensusWithNoReads(self):
        consensus = self.amino.get_consensus()
        
        self.assertEqual(consensus, '-')
        
    def testMissingData(self):
        "Lower-case n represents a gap between the forward and reverse reads."
        
        nuc_seq = 'CTn'
        expected_consensus = 'L'
        
        self.amino.count_nucleotides(nuc_seq, 1)
        consensus = self.amino.get_consensus()
        
        self.assertEqual(expected_consensus, consensus)

class SeedNucleotideTest(unittest.TestCase):
    def setUp(self):
        self.nuc = aln2counts.SeedNucleotide()
    
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
        consensus_max = self.nuc.get_consensus(aln2counts.MAX_CUTOFF)
        consensus_mix = self.nuc.get_consensus(0.1)

        expected_consensus = 'C'
        self.assertEqual(expected_consensus, consensus_max)
        self.assertEqual(expected_consensus, consensus_mix)

    def testConsensusMixed(self):
        self.nuc.count_nucleotides('C', 2)
        self.nuc.count_nucleotides('T', 1)
        consensus_max = self.nuc.get_consensus(aln2counts.MAX_CUTOFF)
        consensus_mix = self.nuc.get_consensus(0.1)

        expected_consensus_max = 'C'
        expected_consensus_mix = 'Y'
        self.assertEqual(expected_consensus_max, consensus_max)
        self.assertEqual(expected_consensus_mix, consensus_mix)

    def testConsensusMixedThree(self):
        self.nuc.count_nucleotides('C', 2)
        self.nuc.count_nucleotides('T', 1)
        self.nuc.count_nucleotides('G', 1)
        consensus_max = self.nuc.get_consensus(aln2counts.MAX_CUTOFF)
        consensus_mix = self.nuc.get_consensus(0.1)

        expected_consensus_max = 'C'
        expected_consensus_mix = 'B' # B is a mix of T, G, and C
        self.assertEqual(expected_consensus_max, consensus_max)
        self.assertEqual(expected_consensus_mix, consensus_mix)

    def testConsensusMixedAll(self):
        self.nuc.count_nucleotides('C', 2)
        self.nuc.count_nucleotides('T', 1)
        self.nuc.count_nucleotides('G', 1)
        self.nuc.count_nucleotides('A', 1)
        consensus_max = self.nuc.get_consensus(aln2counts.MAX_CUTOFF)
        consensus_mix = self.nuc.get_consensus(0.1)

        expected_consensus_max = 'C'
        expected_consensus_mix = 'N' # All four are reported as N
        self.assertEqual(expected_consensus_max, consensus_max)
        self.assertEqual(expected_consensus_mix, consensus_mix)

    def testConsensusMixedMax(self):
        self.nuc.count_nucleotides('C', 2)
        self.nuc.count_nucleotides('T', 2)
        self.nuc.count_nucleotides('G', 1)
        consensus_max = self.nuc.get_consensus(aln2counts.MAX_CUTOFF)
        consensus_mix = self.nuc.get_consensus(0.1)

        expected_consensus_max = 'Y' # C and T tie for max, mix is Y
        expected_consensus_mix = 'B' # C, T, and G mix is B
        self.assertEqual(expected_consensus_max, consensus_max)
        self.assertEqual(expected_consensus_mix, consensus_mix)

    def testConsensusCutoff(self):
        self.nuc.count_nucleotides('C', 2)
        self.nuc.count_nucleotides('T', 1)
        consensus_mix = self.nuc.get_consensus(0.5)

        expected_consensus = 'C' # T was below the cutoff
        self.assertEqual(expected_consensus, consensus_mix)

    def testConsensusCutoffAtBoundary(self):
        self.nuc.count_nucleotides('C', 9000)
        self.nuc.count_nucleotides('T', 1000)
        consensus_mix = self.nuc.get_consensus(0.1)

        expected_consensus = 'Y' # T was at the cutoff
        self.assertEqual(expected_consensus, consensus_mix)

    def testConsensusCutoffBelowBoundary(self):
        self.nuc.count_nucleotides('C', 9001)
        self.nuc.count_nucleotides('T', 999)
        consensus_mix = self.nuc.get_consensus(0.5)

        expected_consensus = 'C' # T was below the cutoff
        self.assertEqual(expected_consensus, consensus_mix)

    def testConsensusMixedWithPoorQuality(self):
        self.nuc.count_nucleotides('N', 2)
        self.nuc.count_nucleotides('T', 1)
        consensus_max = self.nuc.get_consensus(aln2counts.MAX_CUTOFF)
        consensus_mix = self.nuc.get_consensus(0.1)

        expected_consensus_max = 'T' # N always overruled
        expected_consensus_mix = 'T'
        self.assertEqual(expected_consensus_max, consensus_max)
        self.assertEqual(expected_consensus_mix, consensus_mix)

    def testConsensusMixedWithGap(self):
        self.nuc.count_nucleotides('-', 2)
        self.nuc.count_nucleotides('T', 1)
        consensus_max = self.nuc.get_consensus(aln2counts.MAX_CUTOFF)
        consensus_mix = self.nuc.get_consensus(0.1)

        expected_consensus_max = 'T' # dash always overruled
        expected_consensus_mix = 'T'
        self.assertEqual(expected_consensus_max, consensus_max)
        self.assertEqual(expected_consensus_mix, consensus_mix)

    def testConsensusMixedWithGapAndPoorQuality(self):
        self.nuc.count_nucleotides('N', 3)
        self.nuc.count_nucleotides('-', 2)
        self.nuc.count_nucleotides('T', 1)
        consensus_max = self.nuc.get_consensus(aln2counts.MAX_CUTOFF)
        consensus_mix = self.nuc.get_consensus(0.1)

        expected_consensus_max = 'T'
        expected_consensus_mix = 'T'
        self.assertEqual(expected_consensus_max, consensus_max)
        self.assertEqual(expected_consensus_mix, consensus_mix)

    def testConsensusPoorQualityOnly(self):
        self.nuc.count_nucleotides('N', 1)
        consensus_max = self.nuc.get_consensus(aln2counts.MAX_CUTOFF)
        consensus_mix = self.nuc.get_consensus(0.1)

        expected_consensus_max = 'N'
        expected_consensus_mix = 'N'
        self.assertEqual(expected_consensus_max, consensus_max)
        self.assertEqual(expected_consensus_mix, consensus_mix)

    def testConsensusMixedGapAndPoorQualityOnly(self):
        self.nuc.count_nucleotides('N', 3)
        self.nuc.count_nucleotides('-', 2)
        consensus_max = self.nuc.get_consensus(aln2counts.MAX_CUTOFF)
        consensus_mix = self.nuc.get_consensus(0.1)

        expected_consensus_max = 'N'
        expected_consensus_mix = 'N'
        self.assertEqual(expected_consensus_max, consensus_max)
        self.assertEqual(expected_consensus_mix, consensus_mix)

    def testConsensusAllBelowCutoff(self):
        self.nuc.count_nucleotides('C', 101)
        self.nuc.count_nucleotides('T', 100)
        self.nuc.count_nucleotides('G', 99)
        consensus_max = self.nuc.get_consensus(aln2counts.MAX_CUTOFF)
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
        
        #No counts added
        
        consensus_max = self.nuc.get_consensus(aln2counts.MAX_CUTOFF)
        consensus_mix = self.nuc.get_consensus(0.1)

        expected_consensus = ''
        self.assertEqual(expected_consensus, consensus_max)
        self.assertEqual(expected_consensus, consensus_mix)
