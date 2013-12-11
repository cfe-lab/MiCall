import miseq_modules,sys

csf_file = sys.argv[1]
g2p_alignment_cutoff = float(sys.argv[2])

print "miseq_modules.g2p_scoring({}, {})".format(csf_file, g2p_alignment_cutoff)
miseq_modules.g2p_scoring(csf_file, g2p_alignment_cutoff)
