import os,sys
from glob import glob
from proportion_x4 import prop_x4

summary_file = open(root + '/v3prot.summary', 'w')
v3_files = glob(root + '/*.v3prot')
for i, file in enumerate(v3_files):
prefix, gene, qCutoff = file.split('.')[:3]
total_x4_count, total_count = prop_x4(file, fpr_cutoff, mincount)
proportion_x4 = total_x4_count / total_count
summary_file.write("{}".format(proportion_x4))
summary_file.close()
timestamp('Generated v3prot.summary file\n', my_rank, nprocs)
