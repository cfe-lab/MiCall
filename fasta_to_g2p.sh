#!/bin/bash -e

# The #! only allows one argument, so I couldn't get the -rubygems to work.
source $HOME/.rvm/scripts/rvm
ruby -rubygems fasta_to_g2p.rb $*
