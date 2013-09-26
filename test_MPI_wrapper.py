import os
import sys

load_mpi = "module load openmpi/gnu"
script = "/usr/local/share/miseq/scripts/generate_remapped_consensus_and_align_fastas.py"

script = "/usr/local/share/miseq/scripts/generate_remapped_consensus_and_align_fastas_2.py"
script_argument = "/data/miseq/130820_M01841_0018_000000000-A4V8D/"

#script_argument = "/data/miseq/0_testing/"
os.system("{}; mpirun -machinefile mfile python {} {}".format(load_mpi, script, script_argument))
