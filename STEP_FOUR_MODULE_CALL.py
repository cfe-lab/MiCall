import miseq_modules,sys

file = sys.argv[1]
mode = sys.argv[2]
consensus_mixture_cutoffs = map(float,(sys.argv[3].split(",")))
final_alignment_ref_path = sys.argv[4]

print "miseq_modules.csf2counts({},{},{},{})".format(file, mode, consensus_mixture_cutoffs, final_alignment_ref_path)
miseq_modules.csf2counts(file,mode,consensus_mixture_cutoffs,final_alignment_ref_path)
