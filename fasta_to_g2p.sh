#!/bin/bash -e

# The #! only allows one argument, so I couldn't get the -rubygems to work.
if [ -f $HOME/.rvm/scripts/rvm ]; then
    source $HOME/.rvm/scripts/rvm
fi
ruby -rubygems fasta_to_g2p.rb $*
