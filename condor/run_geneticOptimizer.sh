#!/bin/bash

echo "Starting genetic optimizer"

cd /afs/cern.ch/work/j/jniedzie/private/disapp_tracks/disappTracks/
. setenv.sh

echo "Im in `pwd`"
./geneticOptimizer

# > /afs/cern.ch/work/j/jniedzie/private/disapp_tracks/disappTracks/results/geneticOutput_$1_$2.txt
