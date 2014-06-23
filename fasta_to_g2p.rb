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
pssm_lib.rb
g2p.matrix
CSV [Ruby module]
=end

require 'CSV'
require './pssm_lib'

f = File.open(ARGV[1], mode='w') # write-only

CSV.foreach(ARGV[0]) do |row|
  region, qcut, rank, count, offset, seq = row
  if region != 'V3LOOP'
    next
  end
  offset = offset.to_i
  seq2 = '-'*offset + seq
  score, aligned = run_g2p(seq2, $std_v3_g2p, load_matrix('g2p'))
  aligned2 = aligned.flatten * ""
  fpr = g2p_to_fpr(score)
  f.write("#{qcut},#{rank},#{count},#{score},#{fpr},#{aligned2}\n")
end

f.close()
