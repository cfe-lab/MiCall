=begin
Iterate through a CSV comprising tuples of header and nucleotide
sequences of HIV env V3.  Apply Conan's implementation of the 
geno2pheno (g2p) coreceptor tropism prediction algorithm to 
each entry.  Output tuples of header, g2p score and FPR.
=end

require 'CSV'
require './pssm_lib'

CSV.foreach(ARGV[0])

score, aligned = run_g2p(ARGV[0], $std_v3_g2p, load_matrix('g2p'))  
fpr = g2p_to_fpr(score)
puts score, fpr, aligned.join('')
