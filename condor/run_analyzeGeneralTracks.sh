#!/bin/bash

echo "Starting analyzeClusters"

cd /afs/cern.ch/work/j/jniedzie/private/disapp_tracks/disappTracks/
. setenv.sh

./analyzeGeneralTracks $1 $2
