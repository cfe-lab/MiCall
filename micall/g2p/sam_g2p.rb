=begin
Iterate through a CSV (produced by remap.py) where each row contains
 the contents of a SAM file.
Sequences have already been restricted to HIV env V3 if the reads were
 mapped to the 'V3LOOP' region.
Do not attempt to clean sequence based on CIGAR string.
Apply Conan's implementation of the geno2pheno (g2p) coreceptor
 tropism prediction algorithm to each entry.  Output tuples of header,
 g2p score, FPR, and g2p-aligned protein sequence.

Call:
ruby sam_to_g2p.rb <INPUT SAM csv> <INPUT NUC csv> <OUTPUT scored csv>

Dependencies:
BioRuby
pssm_lib.rb
g2p.matrix
CSV [Ruby module]
=end

require 'faster_csv'
require './pssm_lib'
require 'bio'
require 'ostruct'

QMIN = 20   # minimum base quality within insertions
QCUT = 10   # minimum base quality to not be censored
QDELTA = 5

def parse_options()
  options = OpenStruct.new
  options.remap_csv, options.nuc_csv, options.g2p_csv = ARGV
  raise "Usage: #{__FILE__} remap_csv nuc_csv g2p_csv" unless options.g2p_csv
  
  options
end

# Applies a cigar string to recreate a read, then clips the read.
#
# * <tt>:cigar</tt> - a string in the CIGAR format, describing the relationship
#   between the read sequence and the consensus sequence
# * <tt>:seq</tt> - the sequence that was read
# * <tt>:qual</tt> - quality codes for each base in the read
# * <tt>:pos</tt> - first position of the read, given in consensus coordinates
# * <tt>:clip_from</tt> - first position to include after clipping, given in
#   consensus coordinates
# * <tt>:clip_to</tt> - last position to include after clipping, given in
#   consensus coordinates, nil means no clipping at the end
#
# Returns the new sequence and the new quality string. If none of the read was
# within the clipped range, then both strings will be blank.
#
#   apply_cigar_and_clip("12M", "AAACAACCACCC", "AAAAAAAAAAAA", 0, 3, 9)
#   # => "CAACCA", "AAAAAA"
def apply_cigar_and_clip(cigar, seq, qual, pos=0, clip_from=0, clip_to=nil)
  clip_to = clip_to || -1
  newseq = '-' * pos
  newqual = '!' * pos
  is_valid =        /^((\d+)([MIDNSHPX=]))*$/.match(cigar)
  tokens = cigar.scan(/(\d+)([MIDNSHPX=])/)
  raise ArgumentError.new("Invalid CIGAR string: '#{cigar}'.") unless is_valid
  left = 0
  tokens.each do |token|
      length, operation = token
      length = length.to_i
      right = left + length - 1
      if operation == 'M'
          newseq += seq[left..right]
          newqual += qual[left..right]
          left += length
      elsif operation == 'D'
          newseq += '-'*length
          newqual += ' '*length
      elsif operation == 'I'
          its_quals = qual[left..right].unpack('C*') # array of ASCII codes
          if its_quals.map{|asc| asc-33 >= QMIN}.all?
              # accept only high quality insertions
              newseq += seq[left..right]
              newqual += qual[left..right]
              clip_from += length if left < clip_from
              clip_to += length if clip_to >= 0
          end
          left += length
      elsif operation == 'S'
          left += length
      else
          raise ArgumentError.new("Unsupported CIGAR token: '#{token}'.")
      end
      if left > seq.length
        message = "CIGAR string '#{cigar}' is too long for sequence '#{seq}'."
        raise ArgumentError.new(message)
      end
  end
  if left < seq.length
    message = "CIGAR string '#{cigar}' is too short for sequence '#{seq}'."
    raise ArgumentError.new(message)
  end
  newseq = newseq[clip_from..clip_to] || ''
  newqual = newqual[clip_from..clip_to] || ''
  return newseq, newqual
end

def merge_pairs(seq1, seq2, qual1, qual2)
    mseq = ''
    if seq2.length < seq1.length # make sure seq2 is longer
      seq1, seq2 = seq2, seq1 
      qual1, qual2 = qual2, qual1
    end
    i = 0
    qual1_ints = qual1.unpack('C*')
    qual2_ints = qual2.unpack('C*')
    qual1_ints.map! {|x| x - 33}
    qual2_ints.map! {|x| x - 33}
    
    seq2.each_byte { |b2|
        c2 = b2.chr
        q2 = qual2_ints[i]
        if i < seq1.length
            c1 = seq1[i,1]
            q1 = qual1_ints[i]
            if c1 == '-' and c2 == '-'
                mseq += '-'
            elsif c1 == c2
                if q1 > QCUT or q2 > QCUT
                    mseq += c1
                else
                    mseq += 'N'
                end
            else
                if (q2 - q1).abs >= QDELTA
                    if q1 > [q2, QCUT].max
                        mseq += c1
                    elsif q2 > QCUT
                        mseq += c2
                    else
                        mseq += 'N'
                    end
                else
                    mseq += 'N'
                end
            end
        else
            # past the end of read 1
            if c2 == '-' and q2 == 0
                mseq += 'n' # interval between reads
            else
                if q2 > QCUT
                    mseq += c2
                else
                    mseq += 'N'
                end
            end
        end
        i += 1 # increment at end of loop
    }
    return mseq
end

# Keep track of the clipping region for each seed by looking at which query
# positions mapped to the coordinate region.
class RegionTracker
  # Initialize an instance.
  #
  # * <tt>:tracked_region</tt> - the only region to track. All others are
  # ignored.
  def initialize(tracked_region)
    @tracked_region = tracked_region
    @ranges = {}
  end
  
  # Add a nucleotide position to the tracker.
  #
  # * <tt>:seed</tt> - name of the seed region
  # * <tt>:region</tt> - name of the coordinate region
  # * <tt>:query_pos</tt> - query position in the consensus coordinates
  def add_nuc(seed, region, query_pos)
    if region != @tracked_region
      return
    end
    
    range = @ranges[seed]
    if range.nil?
      @ranges[seed] = [query_pos, query_pos]
    else
      if range[1] < query_pos
        range[1] = query_pos
      elsif query_pos < range[0]
        range[0] = query_pos
      end
    end
  end
  
  # Get the minimum and maximum query positions that were seen for a seed.
  # * <tt>:seed</tt> - name of the seed region
  # Returns an array of two integers.
  def get_range(seed)
    return @ranges[seed]
  end
end

def main()
  options = parse_options()
  ## convert file contents into hash to collect identical variants
  
  pairs = Hash.new  # cache read for pairing
  merged = Hash.new # tabulate merged sequence variants
  tracker = RegionTracker.new("V3LOOP")
  
  # Look up clipping region for each seed
  FasterCSV.foreach(
    options.nuc_csv,
    :headers => true,
    :return_headers => false) do |row|

    # Conan's fix: if deletion, then this value is empty string that returns 0 on calling "to_i"
    next if(row['query.nuc.pos'].nil?)
    
    tracker.add_nuc(
      row['seed'],
      row['region'],
      row['query.nuc.pos'].to_i - 1)
  end
  
  FasterCSV.foreach(
      options.remap_csv,
      :headers => true,
      :return_headers => false) do |row|
      
      clip_from, clip_to = tracker.get_range(row['rname'])
      
      if clip_from.nil?
          # uninteresting region
          next
      end
      seq2, qual2 = apply_cigar_and_clip(
        row['cigar'],
        row['seq'],
        row['qual'],
        row['pos'].to_i-1,
        clip_from,
        clip_to)
  
      pair = pairs[row['qname']]
      if pair == nil
        # cache this read
        pairs[row['qname']] = {'seq'=>seq2, 'qual'=>qual2}
      else
        # merge reads
        pairs.delete(row['qname'])
        seq1 = pair['seq'] # retrieve from cache
        qual1 = pair['qual']
      
        mseq = merge_pairs(seq1, seq2, qual1, qual2)
      
        count = merged.fetch(mseq, 0)
        merged[mseq] = count + 1
      end
  end  
  
  sorted = merged.sort_by {|k,v| v}.reverse  # sort by count in descending order
  
  ## apply g2p algorithm to merged reads
  f = File.open(options.g2p_csv, mode='w') # write-only
  f.write("rank,count,g2p,fpr,aligned,error\n")  # CSV header
  
  rank = 0
  sorted.each do |s, count|
      rank += 1
      prefix = "#{rank},#{count}"
      seq = s.delete '-'
      if (seq.upcase.count 'N') > (0.5*seq.length)
          # if more than 50% of sequence is garbage
          f.write("#{prefix},,,,low quality\n")
          next
      end
      
      if s.length % 3 != 0
          f.write("#{prefix},,,,notdiv3\n")
          next
      end
    
      if seq.length == 0
        f.write("#{prefix},,,,zerolength\n")
        next
      end
    
      dna = Bio::Sequence::NA.new(seq)
      prot = dna.translate
  
      # sanity check 1 - bounded by cysteines
      if !prot.match(/^C/) || !prot.match(/C$/)
          f.write("#{prefix},,,#{prot},cysteines\n")
          next
      end
  
      # sanity check 2 - no ambiguous codons
      if prot.count('X') > 0
          f.write("#{prefix},,,#{prot},ambiguous\n")
          next
      end
  
      # sanity check 3 - no stop codons
      if prot.count('*') > 0
          f.write("#{prefix},,,#{prot},stop codons\n")
          next
      end
  
      # sanity check 4 - V3 length in range 32-40 inclusive
      if prot.length < 32 || prot.length > 40
          f.write("#{prefix},,,#{prot},length\n")
          next
      end
  
      score, aligned = run_g2p(seq, $std_v3_g2p, load_matrix('g2p'))
      begin
          aligned2 = aligned.flatten * ""
      rescue
          # sequence failed to align
          f.write("#{prefix},#{score},,,failed to align\n")
          next
      end
      fpr = g2p_to_fpr(score)
      f.write("#{prefix},#{score},#{fpr},#{aligned2},\n")
  end
  
  f.close()
end

if __FILE__ == $0
  main()
end
