"""
pairAlign v1.2

pairwise alignment support functions for g2pscorer.py

revision log
- May 10, 2011 - was bouncing alignments with score of 0 - in fact there is a chance that a correct
				alignment can have a score of exactly 0 when terminal gap penalty is active.
				- removed this, errors will be raised downstream by empty result

- May 17, 2011 - replace empirical HIV 25% matrix with gonnet matrix from g2p
				- also use g2p gap penalties
	--> v1.2
	
- May 31, 2011 - found bug in pairwise_align(), aligned_ref and aligned_query were
					swapped in parsing alignment output
"""



import HyPhy
from seqUtils import *
import re

gap_prefix = re.compile('^[-]+')
gap_suffix = re.compile('[-]+$')

# =====================================================
# align protein sequences


scoreMatrixGonnet = """\
{{2.4,-0.6,-0.3,-0.3,0.5,0.0,-0.2,0.5,-0.8,-0.8,-1.2,-0.4,-0.7,-2.3,0.3,1.1,0.6,-3.6,-2.2,0.1,-5.0,-8.0},\
{-0.6,4.7,0.3,-0.3,-2.2,0.4,1.5,-1.0,0.6,-2.4,-2.2,2.7,-1.7,-3.2,-0.9,-0.2,-0.2,-1.6,-1.8,-2.0,-5.0,-8.0},\
{-0.3,0.3,3.8,2.2,-1.8,0.9,0.7,0.4,1.2,-2.8,-3.0,0.8,-2.2,-3.1,-0.9,0.9,0.5,-3.6,-1.4,-2.2,-5.0,-8.0},\
{-0.3,-0.3,2.2,4.7,-3.2,2.7,0.9,0.1,0.4,-3.8,-4.0,0.5,-3.0,-4.5,-0.7,0.5,0.0,-5.2,-2.8,-2.9,-5.0,-8.0},\
{0.5,-2.2,-1.8,-3.2,11.5,-3.0,-2.4,-2.0,-1.3,-1.1,-1.5,-2.8,-0.9,-0.8,-3.1,0.1,-0.5,-1.0,-0.5,0.0,-5.0,-8.0},\
{0.0,0.4,0.9,2.7,-3.0,3.6,1.7,-0.8,0.4,-2.7,-2.8,1.2,-2.0,-3.9,-0.5,0.2,-0.1,-4.3,-2.7,-1.9,-5.0,-8.0},\
{-0.2,1.5,0.7,0.9,-2.4,1.7,2.7,-1.0,1.2,-1.9,-1.6,1.5,-1.0,-2.6,-0.2,0.2,0.0,-2.7,-1.7,-1.5,-5.0,-8.0},\
{0.5,-1.0,0.4,0.1,-2.0,-0.8,-1.0,6.6,-1.4,-4.5,-4.4,-1.1,-3.5,-5.2,-1.6,0.4,-1.1,-4.0,-4.0,-3.3,-5.0,-8.0},\
{-0.8,0.6,1.2,0.4,-1.3,0.4,1.2,-1.4,6.0,-2.2,-1.9,0.6,-1.3,-0.1,-1.1,-0.2,-0.3,-0.8,2.2,-2.0,-5.0,-8.0},\
{-0.8,-2.4,-2.8,-3.8,-1.1,-2.7,-1.9,-4.5,-2.2,4.0,2.8,-2.1,2.5,1.0,-2.6,-1.8,-0.6,-1.8,-0.7,3.1,-5.0,-8.0},\
{-1.2,-2.2,-3.0,-4.0,-1.5,-2.8,-1.6,-4.4,-1.9,2.8,4.0,-2.1,2.8,2.0,-2.3,-2.1,-1.3,-0.7,0.0,1.8,-5.0,-8.0},\
{-0.4,2.7,0.8,0.5,-2.8,1.2,1.5,-1.1,0.6,-2.1,-2.1,3.2,-1.4,-3.3,-0.6,0.1,0.1,-3.5,-2.1,-1.7,-5.0,-8.0},\
{-0.7,-1.7,-2.2,-3.0,-0.9,-2.0,-1.0,-3.5,-1.3,2.5,2.8,-1.4,4.3,1.6,-2.4,-1.4,-0.6,-1.0,-0.2,1.6,-5.0,-8.0},\
{-2.3,-3.2,-3.1,-4.5,-0.8,-3.9,-2.6,-5.2,-0.1,1.0,2.0,-3.3,1.6,7.0,-3.8,-2.8,-2.2,3.6,5.1,0.1,-5.0,-8.0},\
{0.3,-0.9,-0.9,-0.7,-3.1,-0.5,-0.2,-1.6,-1.1,-2.6,-2.3,-0.6,-2.4,-3.8,7.6,0.4,0.1,-5.0,-3.1,-1.8,-5.0,-8.0},\
{1.1,-0.2,0.9,0.5,0.1,0.2,0.2,0.4,-0.2,-1.8,-2.1,0.1,-1.4,-2.8,0.4,2.2,1.5,-3.3,-1.9,-1.0,-5.0,-8.0},\
{0.6,-0.2,0.5,0.0,-0.5,-0.1,0.0,-1.1,-0.3,-0.6,-1.3,0.1,-0.6,-2.2,0.1,1.5,2.5,-3.5,-1.9,0.0,-5.0,-8.0},\
{-3.6,-1.6,-3.6,-5.2,-1.0,-4.3,-2.7,-4.0,-0.8,-1.8,-0.7,-3.5,-1.0,3.6,-5.0,-3.3,-3.5,14.2,4.1,-2.6,-5.0,-8.0},\
{-2.2,-1.8,-1.4,-2.8,-0.5,-2.7,-1.7,-4.0,2.2,-0.7,0.0,-2.1,-0.2,5.1,-3.1,-1.9,-1.9,4.1,7.8,-1.1,-5.0,-8.0},\
{0.1,-2.0,-2.2,-2.9,0.0,-1.9,-1.5,-3.3,-2.0,3.1,1.8,-1.7,1.6,0.1,-1.8,-1.0,0.0,-2.6,-1.1,3.4,-5.0,-8.0},\
{-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,1.0,-5.0},\
{-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,1.0}};
"""

scoreMatrixHIV25 = """\
{{7,-7,-7,-4,-10,-11,-4,-3,-10,-6,-9,-9,-7,-13,-3,-2,1,-16,-15,0,-5,-5,-10,-17},\
{-7,7,-5,-11,-8,-2,-7,-2,0,-6,-6,2,-3,-12,-4,-2,-2,-5,-9,-10,-7,-3,-10,-17},\
{-7,-5,8,2,-9,-6,-6,-7,0,-6,-12,0,-10,-12,-9,1,0,-17,-3,-10,6,-6,-10,-17},\
{-4,-11,2,8,-14,-10,0,-2,-3,-11,-15,-7,-13,-15,-13,-5,-6,-16,-6,-5,7,0,-10,-17},\
{-10,-8,-9,-14,11,-16,-15,-5,-7,-11,-9,-13,-14,0,-12,-1,-6,-2,0,-8,-10,-16,-10,-17},\
{-11,-2,-6,-10,-16,8,-2,-10,0,-12,-4,0,-8,-12,-1,-9,-8,-14,-9,-13,-7,6,-10,-17},\
{-4,-7,-6,0,-15,-2,7,-1,-9,-12,-15,-1,-10,-17,-13,-11,-8,-15,-12,-5,0,6,-10,-17},\
{-3,-2,-7,-2,-5,-10,-1,7,-10,-11,-14,-6,-12,-9,-11,-1,-7,-5,-14,-5,-4,-3,-10,-17},\
{-10,0,0,-3,-7,0,-9,-10,10,-10,-4,-5,-10,-6,-3,-6,-6,-11,2,-14,-1,-2,-10,-17},\
{-6,-6,-6,-11,-11,-12,-12,-11,-10,7,0,-7,0,-2,-10,-4,0,-14,-9,2,-7,-12,-10,-17},\
{-9,-6,-12,-15,-9,-4,-15,-14,-4,0,6,-10,0,0,-3,-5,-8,-6,-8,-4,-13,-6,-10,-17},\
{-9,2,0,-7,-13,0,-1,-6,-5,-7,-10,7,-4,-14,-9,-5,-1,-12,-13,-9,-1,-1,-10,-17},\
{-7,-3,-10,-13,-14,-8,-10,-12,-10,0,0,-4,10,-7,-11,-9,-1,-11,-15,0,-11,-9,-10,-17},\
{-13,-12,-12,-15,0,-12,-17,-9,-6,-2,0,-14,-7,10,-11,-5,-10,-5,1,-5,-13,-14,-10,-17},\
{-3,-4,-9,-13,-12,-1,-13,-11,-3,-10,-3,-9,-11,-11,8,-1,-3,-13,-11,-12,-10,-3,-10,-17},\
{-2,-2,1,-5,-1,-9,-11,-1,-6,-4,-5,-5,-9,-5,-1,8,0,-12,-6,-9,0,-10,-10,-17},\
{1,-2,0,-6,-6,-8,-8,-7,-6,0,-8,-1,-1,-10,-3,0,7,-16,-10,-4,-2,-8,-10,-17},\
{-16,-5,-17,-16,-2,-14,-15,-5,-11,-14,-6,-12,-11,-5,-13,-12,-16,10,-4,-16,-16,-14,-10,-17},\
{-15,-9,-3,-6,0,-9,-12,-14,2,-9,-8,-13,-15,1,-11,-6,-10,-4,10,-12,-4,-10,-10,-17},\
{0,-10,-10,-5,-8,-13,-5,-5,-14,2,-4,-9,0,-5,-12,-9,-4,-16,-12,7,-7,-7,-10,-17},\
{-5,-7,6,7,-10,-7,0,-4,-1,-7,-13,-1,-11,-13,-10,0,-2,-16,-4,-7,7,-2,-10,-17},\
{-5,-3,-6,0,-16,6,6,-3,-2,-12,-6,-1,-9,-14,-3,-10,-8,-14,-10,-7,-2,6,-10,-17},\
{-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,1,-17},\
{-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,1}};
"""


scoreMatrixNuc = """\
{{ 5,-4,-4,-4},\
 {-4, 5,-4,-4},\
 {-4,-4, 5,-4},\
 {-4,-4,-4, 5}};
"""

def init_hyphy_nuc_align (hy):
	# create associative array with alignment options
	#  second argument is set to False for object persistence
	hy.ExecuteBF("alignOptions = {};", False)
	hy.ExecuteBF("alignOptions [\"SEQ_ALIGN_CHARACTER_MAP\"]=\"ACGT\";", False)
	hy.ExecuteBF("alignOptions [\"SEQ_ALIGN_SCORE_MATRIX\"] = "+scoreMatrixNuc, False)
	hy.ExecuteBF("alignOptions [\"SEQ_ALIGN_GAP_OPEN\"] = 10;", False)
	hy.ExecuteBF("alignOptions [\"SEQ_ALIGN_GAP_OPEN2\"] = 5;", False)
	hy.ExecuteBF("alignOptions [\"SEQ_ALIGN_GAP_EXTEND\"] = 1;", False)
	hy.ExecuteBF("alignOptions [\"SEQ_ALIGN_GAP_EXTEND2\"] = 1;", False)
	hy.ExecuteBF("alignOptions [\"SEQ_ALIGN_AFFINE\"] = 1;", False)
	# do not penalize prefix or suffix indels
	hy.ExecuteBF("alignOptions [\"SEQ_ALIGN_NO_TP\"] = 1;", False)



def init_hyphy_protein_align (hyphy):
	hyphy.ExecuteBF("alignOptions = {};", False)
	hyphy.ExecuteBF("alignOptions [\"SEQ_ALIGN_CHARACTER_MAP\"]=\"ARNDCEQGHILKMFPSTWYVX*\";", False)
	hyphy.ExecuteBF("alignOptions [\"SEQ_ALIGN_SCORE_MATRIX\"] = "+scoreMatrixGonnet, False)
	hyphy.ExecuteBF("alignOptions [\"SEQ_ALIGN_GAP_OPEN\"] = 20;", False)
	hyphy.ExecuteBF("alignOptions [\"SEQ_ALIGN_GAP_OPEN2\"] = 5;", False)
	hyphy.ExecuteBF("alignOptions [\"SEQ_ALIGN_GAP_EXTEND\"] = 2;", False)
	hyphy.ExecuteBF("alignOptions [\"SEQ_ALIGN_GAP_EXTEND2\"] = 1;", False)
	hyphy.ExecuteBF("alignOptions [\"SEQ_ALIGN_AFFINE\"] = 1;", False)
	hyphy.ExecuteBF("alignOptions [\"SEQ_ALIGN_NO_TP\"] = 0;", False)




def pairwise_align (query_seq, ref_seq, hyphy):
	"""
	Returns a tuple containing aligned query and reference sequences using
	Smith-Wasserman algorithm.
	alignOptions is a persistent HyPhy object initialized in init_hyphy_nuc_align.
	"""
	#hy.ExecuteBF("inSeq = {{\""+ref_seq+"\",\""+query_seq+"\"}};", False)
	#hy.ExecuteBF("AlignSequences(result,inSeq,alignOptions);", False)
	
	dump = hyphy.ExecuteBF ('inStr={{"'+ref_seq+'","'+query_seq+'"}};', False);
	dump = hyphy.ExecuteBF ('AlignSequences(aligned, inStr, alignOptions);', False);
	aligned = hyphy.ExecuteBF ('return aligned;', False);
	exec "d = " + aligned.sData
	
	align_score = int(d['0']['0'])
	aligned_ref = d['0']['1']
	aligned_query = d['0']['2']
	
	return (aligned_query, aligned_ref, align_score)
	

def get_boundaries (str):
	# return a tuple giving indices of subsequence without gap prefix and suffix
	res = [0,len(str)]
	left = gap_prefix.findall(str)
	right = gap_suffix.findall(str)
	if left:
		res[0] = len(left[0])
	if right:
		res[1] = len(str) - len(right[0])
		
	return res


def get_overlap (str1, str2):
	bounds1 = get_boundaries(str1)
	bounds2 = get_boundaries(str2)
	left_bound = max(bounds1[0], bounds2[0])
	right_bound = min(bounds1[1], bounds2[1])
	return (left_bound, right_bound)



# ===========================================

def align (hyphy, sd, refseq):
	# screen sequence for likeness with V3 - determine direction and reading frame
	init_hyphy_protein_align(hyphy)
	
	
	# loop through headers
	for h in sd.iterkeys():
		# we trust that the sequence is in-frame, but we don't necessarily
		# trust the alignment
		aligned_query, aligned_ref, new_seq = codon_align (hyphy, sd[h]['rawseq'], refseq)
		
		sd[h].update({'aligned_ref':aligned_ref, 'aligned_query':aligned_query, 'clipped_nuc':new_seq, 'clipped_aa':query_v3, 'offset':best_offset, 'is_rc':best_is_rc, 'align_score':align_score})
	
	return sd



def codon_align (hyphy, seq, refseq):
	"""
	Align a nucleotide sequence [seq] to a protein reference sequence [refseq]
	by finding the best reading frame of the query sequence, performing the 
	protein alignment and then using the alignment to pad the nucleotide sequence 
	so that it stays in reading frame.
	"""
	
	# strip out gaps from query
	seq = seq.strip('-').replace('---','')
	
	refprot = translate_nuc(refseq, 0)
	
	best_alignment = {'0':{'0':-999999}}
	best_offset = 0
	best_is_rc = False
	
	init_hyphy_protein_align (hyphy)
	
	# attempt all 3 reading frames in forward and reverse directions
	for offset in range(3):
		for is_rc in [False, True]:
			if is_rc:
				aaseq = translate_nuc(reverse_and_complement(seq), offset, True)
			else:
				aaseq = translate_nuc(seq, offset, True)
			
			dump = hyphy.ExecuteBF ('inStr={{"'+refprot+'","'+aaseq+'"}};', False)
			dump = hyphy.ExecuteBF ('AlignSequences(aligned, inStr, alignOptions);', False)
			aligned = hyphy.ExecuteBF ('return aligned;', False)
			exec "d = " + aligned.sData
			
			if d['0']['0'] > best_alignment['0']['0']:
				best_alignment = d
				best_is_rc = is_rc
				best_offset = offset
	
	# 0 = score, 1 = reference, 2 = query
	align_score = int(best_alignment['0']['0'])
	try:
		aligned_ref = best_alignment['0']['1']
	except:
		print best_alignment
		print aaseq
		raise
	
	aligned_query = best_alignment['0']['2']
	
	# clip to overlapping region
	gp_ref = gap_prefix.findall(aligned_ref)
	gp_query = gap_prefix.findall(aligned_query)
	if len(gp_ref) > 0 and len(gp_query) > 0:
		left = max(len(gp_ref[0]), len(gp_query[0]))
	elif len(gp_ref) > 0:
		left = len(gp_ref[0])
	elif len(gp_query) > 0:
		left = len(gp_query[0])
	else:
		left = 0
		
	gs_ref = gap_suffix.findall(aligned_ref)
	gs_query = gap_suffix.findall(aligned_query)
	right = len(aligned_ref)
	if len(gs_ref) > 0 and len(gs_query) > 0:
		right -= max(len(gs_ref[0]), len(gs_query[0]))
	elif len(gs_ref) > 0:
		right -= len(gs_ref[0])
	elif len(gs_query) > 0:
		right -= len(gs_query[0])
	else:
		pass
	
	
	clipped_ref = aligned_ref[ left : right ]
	clipped_query = aligned_query[ left : right ]
	
	
	# apply AA alignment pattern back to the original codon sequence
	# i.e., pad with deletions relative to reference
	codon_query = '-'*best_offset + seq
	codon_query = codon_query[(3*left):(3*right)]
	codon_ref = refseq[(3*left):(3*right)]
	
	new_seq = ''
	old_pos = 0
	for pos in range(len(clipped_query)):
		codon = codon_query[(3*old_pos):(3*(old_pos+1))]
		if clipped_query[pos] == '-' and codon != '---':
			new_seq += '---'
		else:
			new_seq += codon
			old_pos += 1
	
	aligned_query = new_seq
	
	
	new_seq = ''
	old_pos = 0
	for pos in range(len(clipped_ref)):
		codon = codon_ref[(3*old_pos):(3*(old_pos+1))]
		if clipped_query[pos] == '-' and codon != '---':
			new_seq += '---'
		else:
			new_seq += codon
			old_pos += 1
	
	aligned_ref = new_seq
	
	# re-translate leaving in the '?'s
	#aligned_query = translate_nuc(new_seq, 0)
	
	return (aligned_query, aligned_ref, align_score)



# ====================================================================
	
def what_is (hyphy, varname):
	# a debugging utility for retrieving objects from HyPhy instance
	dump = hyphy.ExecuteBF('fprintf(stdout, '+varname+');', False)
	hyout = hyphy.GetStdout();
	return hyout.sData

	
