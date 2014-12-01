require 'test/unit'
require 'sam_g2p'

class CigarTest < Test::Unit::TestCase
  def test_trivial
    cigar = '9M'
    seq     = 'AAACAACCA'
    quality = 'BBBBBBBBB'
    expected_seq = seq
    expected_quality = quality
    
    clipped_seq, clipped_quality = apply_cigar_and_clip(cigar, seq, quality)
    
    assert_equal(expected_seq, clipped_seq)
    assert_equal(expected_quality, clipped_quality)
  end

  def test_deletion
    cigar = '6M3D3M'
    seq              = 'AAACAACCA'
    quality          = 'BBBDDDEEE'
    expected_seq     = 'AAACAA---CCA'
    expected_quality = 'BBBDDD   EEE'
    
    clipped_seq, clipped_quality = apply_cigar_and_clip(cigar, seq, quality)
    
    assert_equal(expected_seq, clipped_seq)
    assert_equal(expected_quality, clipped_quality)
  end

  def test_soft_clip
    cigar = '3S6M'
    seq              = 'AAACAACCA'
    quality          = 'BBBDDDEEE'
    expected_seq     =    'CAACCA'
    expected_quality =    'DDDEEE'
    
    clipped_seq, clipped_quality = apply_cigar_and_clip(cigar, seq, quality)
    
    assert_equal(expected_seq, clipped_seq)
    assert_equal(expected_quality, clipped_quality)
  end

  def test_insertion
    cigar = '3M3I6M'
    seq              = 'AAACAACCACCC'
    quality          = 'BBBDDDEEEFFF'
    expected_seq     = 'AAACAACCACCC'
    expected_quality = 'BBBDDDEEEFFF'
    
    clipped_seq, clipped_quality = apply_cigar_and_clip(cigar, seq, quality)
    
    assert_equal(expected_seq, clipped_seq)
    assert_equal(expected_quality, clipped_quality)
  end

  def test_insertion_low_quality
    cigar = '3M3I6M'
    seq              = 'AAACAACCACCC'
    quality          = 'BBBD*DEEEFFF'
    expected_seq     = 'AAACCACCC'
    expected_quality = 'BBBEEEFFF'
    
    clipped_seq, clipped_quality = apply_cigar_and_clip(cigar, seq, quality)
    
    assert_equal(expected_seq, clipped_seq)
    assert_equal(expected_quality, clipped_quality)
  end
  
  def test_large_token
    cigar = '12M'
    seq     = 'AAACAACCACCC'
    quality = 'BBBBBBBBBBBB'
    expected_seq = seq
    expected_quality = quality
    
    clipped_seq, clipped_quality = apply_cigar_and_clip(cigar, seq, quality)
    
    assert_equal(expected_seq, clipped_seq)
    assert_equal(expected_quality, clipped_quality)
  end
  
  def test_padding
    cigar = '12M'
    seq              = 'AAACAACCACCC'
    quality          = 'BBBDDDEEEFFF'
    pos = 3
    expected_seq     = '---AAACAACCACCC'
    expected_quality = '!!!BBBDDDEEEFFF'
    
    clipped_seq, clipped_quality = apply_cigar_and_clip(cigar, seq, quality, pos)
    
    assert_equal(expected_seq, clipped_seq)
    assert_equal(expected_quality, clipped_quality)
  end
  
  def test_clipping
    cigar = '12M'
    seq              = 'AAACAACCACCC'
    quality          = 'BBBDDDEEEFFF'
    pos = 0
    clip_from = 3
    clip_to = 8
    expected_seq     = 'CAACCA'
    expected_quality = 'DDDEEE'
    
    clipped_seq, clipped_quality = apply_cigar_and_clip(
      cigar,
      seq,
      quality,
      pos,
      clip_from,
      clip_to)
    
    assert_equal(expected_seq, clipped_seq)
    assert_equal(expected_quality, clipped_quality)
  end

  def test_invalid_cigar
    cigar = '3M...6M'
    seq     = 'AAACAACCACCC'
    quality = 'BBBDDDEEEFFF'
    
    begin
      clipped_seq, clipped_quality = apply_cigar_and_clip(cigar, seq, quality)
      
      error_message = nil
    rescue ArgumentError => e
      error_message = e.message
    end
    
    assert_equal("Invalid CIGAR string: '3M...6M'.", error_message)
  end

  def test_unsupported_cigar_token
    cigar = '3M3X6M'
    seq     = 'AAACAACCACCC'
    quality = 'BBBDDDEEEFFF'
    
    begin
      clipped_seq, clipped_quality = apply_cigar_and_clip(cigar, seq, quality)
      
      error_message = nil
    rescue ArgumentError => e
      error_message = e.message
    end
    
    assert_equal("Unsupported CIGAR token: '3X'.", error_message)
  end

  def test_short_cigar
    cigar = '8M'
    seq     = 'AAACAACCA'
    quality = 'BBBDDDEEE'
    
    begin
      clipped_seq, clipped_quality = apply_cigar_and_clip(cigar, seq, quality)
      
      error_message = nil
    rescue ArgumentError => e
      error_message = e.message
    end
    
    assert_equal(
      "CIGAR string '8M' is too short for sequence 'AAACAACCA'.",
      error_message)
  end

  def test_long_cigar
    cigar = '10M'
    seq     = 'AAACAACCA'
    quality = 'BBBDDDEEE'
    
    begin
      clipped_seq, clipped_quality = apply_cigar_and_clip(cigar, seq, quality)
      
      error_message = nil
    rescue ArgumentError => e
      error_message = e.message
    end
    
    assert_equal(
      "CIGAR string '10M' is too long for sequence 'AAACAACCA'.",
      error_message)
  end
end

class MergeTest < Test::Unit::TestCase
  def test_trivial
    seq1  = "AAACAA"
    seq2  = seq1
    qual1 = "BBBDDD"
    qual2 = qual1
    expected_seq = seq1
    
    merged_seq = merge_pairs(seq1, seq2, qual1, qual2)
    
    assert_equal(expected_seq, merged_seq)
  end
  
  def test_overlap
    seq1         = "AAACAA"
    seq2         = "---CAACCA"
    qual1        = "BBBDDD"
    qual2        = "!!!DDDEEE"
    expected_seq = "AAACAACCA"
    
    merged_seq = merge_pairs(seq1, seq2, qual1, qual2)
    
    assert_equal(expected_seq, merged_seq)
  end
  
  # Test a bunch of combinations when the two reads have the same base.
  # A quality of '+' is below the quality cutoff, and a ',' is just above it.
  def test_matching_overlap_with_quality_combinations
    seq1         = "AAAAA"
    seq2         = "AAAAA"
    qual1        = "DD++,"
    qual2        = "D+D+,"
    expected_seq = "AAANA"
    
    merged_seq = merge_pairs(seq1, seq2, qual1, qual2)
    
    assert_equal(expected_seq, merged_seq)
  end
  
  # Test a bunch of combinations when the two reads differ.
  # A quality score of 'H' isn't high enough to outvote 'D', but 'I' is.
  def test_conflict_with_quality_combinations
    seq1         = "GGGGGGGG"
    seq2         = "TTTTTTTT"
    qual1        = "D++#DDDI"
    qual2        = "+D++DHID"
    expected_seq = "GTNNNNTG"
    
    merged_seq = merge_pairs(seq1, seq2, qual1, qual2)
    
    assert_equal(expected_seq, merged_seq)
  end
  
  def test_gap
    seq1         = "AAA"
    seq2         = "------CCA"
    qual1        = "BBBDDD"
    qual2        = "!!!!!!EEE"
    expected_seq = "AAAnnnCCA"
    
    merged_seq = merge_pairs(seq1, seq2, qual1, qual2)
    
    assert_equal(expected_seq, merged_seq)
  end
  
  def test_deletion_not_gap
    seq1         = "AAA"
    seq2         = "AAA---CCA"
    qual1        = "BBBDDD"
    qual2        = "BBB   EEE"
    expected_seq = "AAANNNCCA"
    
    merged_seq = merge_pairs(seq1, seq2, qual1, qual2)
    
    assert_equal(expected_seq, merged_seq)
  end
  
  def test_offset
    seq1         = "---AAACAACCA"
    seq2         = "------CAACCACCC"
    qual1        = "!!!BBBDDDEEE"
    qual2        = "!!!!!!DDDEEEFFF"
    expected_seq = "---AAACAACCACCC"
    
    merged_seq = merge_pairs(seq1, seq2, qual1, qual2)
    
    assert_equal(expected_seq, merged_seq)
  end
  
  # Pass in the reverse read first
  def test_reverse_first
    seq1         = "---CAACCA"
    seq2         = "AAACAA"
    qual1        = "!!!DDDEEE"
    qual2        = "BBBDDD"
    expected_seq = "AAACAACCA"
    
    merged_seq = merge_pairs(seq1, seq2, qual1, qual2)
    
    assert_equal(expected_seq, merged_seq)
  end
end
