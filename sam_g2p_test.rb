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
  
  def test_clip_insertion
    cigar = '6M3I6M'
    seq              = 'AAACAAGGGCCACCC'
    quality          = 'BBBDDDHHHEEEFFF'
    pos = 0
    clip_from = 3
    clip_to = 8
    expected_seq     = 'CAAGGGCCA'
    expected_quality = 'DDDHHHEEE'
    
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
  
  def test_clip_insertion_low_quality
    cigar = '6M3I6M'
    seq              = 'AAACAAGGGCCACCC'
    quality          = 'BBBDDDHH*EEEFFF'
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
  
  def test_insertion_before_clip
    cigar = '3M3I9M'
    seq              = 'AAAGGGCAACCACCC'
    quality          = 'BBBHHHDDDEEEFFF'
    pos = 0
    clip_from = 6
    clip_to = 8
    expected_seq     = 'CCA'
    expected_quality = 'EEE'
    
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
  
  def test_clipping_everything
    cigar = '12M'
    seq              = 'AAACAACCACCC'
    quality          = 'BBBDDDEEEFFF'
    pos = 0
    clip_from = 100
    clip_to = 108
    expected_seq     = ''
    expected_quality = ''
    
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

class RegionTrackerTest < Test::Unit::TestCase
  def setup
    @tracked_region = "R1"
    @untracked_region = "R2"
    @tracker = RegionTracker.new(@tracked_region)
  end
  
  def test_get_range
    sample_name = "S1"
    seed = "R-seed"
    region = @tracked_region
    query_pos1 = 10
    query_pos2 = 11
    query_pos3 = 12
    expected_min = 10
    expected_max = 12
    
    @tracker.add_nuc(sample_name, seed, region, query_pos1)
    @tracker.add_nuc(sample_name, seed, region, query_pos2)
    @tracker.add_nuc(sample_name, seed, region, query_pos3)
    min, max = @tracker.get_range(seed)
    
    assert_equal(expected_min, min, "minimum position")
    assert_equal(expected_max, max, "maximum position")
  end
  
  def test_untracked_region
    sample_name = "S1"
    seed = "R-seed"
    region = @untracked_region
    query_pos1 = 10
    query_pos2 = 11
    query_pos3 = 12
    expected_min = nil
    expected_max = nil
    
    @tracker.add_nuc(sample_name, seed, region, query_pos1)
    @tracker.add_nuc(sample_name, seed, region, query_pos2)
    @tracker.add_nuc(sample_name, seed, region, query_pos3)
    min, max = @tracker.get_range(seed)
    
    assert_equal(expected_min, min, "minimum position")
    assert_equal(expected_max, max, "maximum position")
  end
  
  def test_changed_sample
    sample_name = "S1"
    changed_sample_name = "S2"
    expected_error = "Two sample names found: 'S1' and 'S2'."
    seed = "R-seed"
    region = @tracked_region
    query_pos1 = 10
    query_pos2 = 11
    query_pos3 = 12

    @tracker.add_nuc(sample_name, seed, region, query_pos1)
    @tracker.add_nuc(sample_name, seed, region, query_pos2)
    begin
      @tracker.add_nuc(changed_sample_name, seed, region, query_pos3)
      
      message = nil
    rescue ArgumentError => e
      message = e.message
    end
    
    assert_equal(expected_error, message)
  end
  
  def test_multiple_seeds
    sample_name = "S1"
    seeda = "Ra-seed"
    seedb = "Rb-seed"
    region = @tracked_region
    query_pos1a = 10
    query_pos2a = 11
    query_pos1b = 100
    query_pos2b = 101
    expected_mina = 10
    expected_maxa = 11
    expected_minb = 100
    expected_maxb = 101
    
    @tracker.add_nuc(sample_name, seeda, region, query_pos1a)
    @tracker.add_nuc(sample_name, seeda, region, query_pos2a)
    @tracker.add_nuc(sample_name, seedb, region, query_pos1b)
    @tracker.add_nuc(sample_name, seedb, region, query_pos2b)
    mina, maxa = @tracker.get_range(seeda)
    minb, maxb = @tracker.get_range(seedb)
    
    assert_equal(expected_mina, mina, "minimum position for seed A")
    assert_equal(expected_maxa, maxa, "maximum position for seed A")
    assert_equal(expected_minb, minb, "minimum position for seed B")
    assert_equal(expected_maxb, maxb, "maximum position for seed B")
  end
end