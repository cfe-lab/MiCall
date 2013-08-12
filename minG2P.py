"""
just a wrapper for Conan's Ruby script in Python
"""

import sys, os, math
import subprocess

def conan_g2p (aaseq):
	p = subprocess.Popen(['/usr/local/bin/ruby', 'g2p.rb', aaseq], stdout=subprocess.PIPE)
	stdout, stderr = p.communicate()
	try:
		g2p, fpr, aligned = stdout.split('\n')[:3]
	except ValueError:
		return (None, None, None)
	
	return (g2p, fpr, aligned)

