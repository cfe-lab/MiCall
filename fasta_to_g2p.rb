=begin
Iterate through a CSV (produced by sam2csf) where each row contains
region2, qcut2, rank, count, offset, seq
Restrict to sequences of HIV env V3 (region == V3LOOP).
Apply Conan's implementation of the geno2pheno (g2p) coreceptor
tropism prediction algorithm to each entry.  Output tuples of header,
g2p score, FPR, and g2p-aligned protein sequence.

Call:
ruby fasta_to_g2p.rb <INPUT fasta csv> <OUTPUT scored csv>

Dependencies:
BioRuby
pssm_lib.rb
g2p.matrix
CSV [Ruby module]
=end

require 'csv'
require './pssm_lib'
require 'bio'

f = File.open(ARGV[1], mode='w') # write-only

CSV.foreach(ARGV[0]) do |row|
  region, qcut, rank, count, offset, seq = row
  if region != 'V3LOOP'
    # skip header and other regions
    next
  end
  
  # skip partial sequences
  offset = offset.to_i
  if offset > 0
    f.write("#{qcut},#{rank},#{count},,,,incomplete\n")
    next
  end
  
  # remove all gaps
  seq2 = seq.delete '-'
  if seq2.length % 3 != 0
    # sequence no longer in frame, ignore
    f.write("#{qcut},#{rank},#{count},,,,frameshift\n")
    next
  end
  
  dna = Bio::Sequence.auto(seq2)
  prot = dna.translate
  
  # sanity check 1 - bounded by cysteines
  if !prot.match(/^C/) || !prot.match(/C$/)
    f.write("#{qcut},#{rank},#{count},,,,cysteines\n")
    next
  end
  
  # sanity check 2 - no ambiguous codons
  if prot.count('X') > 0
    f.write("#{qcut},#{rank},#{count},,,,ambiguous\n")
    next
  end
  
  # sanity check 3 - no stop codons
  if prot.count('*') > 0
    f.write("#{qcut},#{rank},#{count},,,,stop codons\n")
    next
  end
  
  # sanity check 4 - V3 length in range 32-40 inclusive
  if prot.length < 32 || prot.length > 40
    f.write("#{qcut},#{rank},#{count},,,,length\n")
    next
  end
  
  score, aligned = run_g2p(seq, $std_v3_g2p, load_matrix('g2p'))
  aligned2 = aligned.flatten * ""
  fpr = g2p_to_fpr(score)
  f.write("#{qcut},#{rank},#{count},#{score},#{fpr},#{aligned2},\n")
end

f.close()
