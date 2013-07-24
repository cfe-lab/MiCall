require '/Users/apoon/Python/pssm_lib'
score, aligned = run_g2p(ARGV[0], $std_v3_g2p, load_matrix('g2p'))  
fpr = g2p_to_fpr(score)
puts score, fpr, aligned.join('')
