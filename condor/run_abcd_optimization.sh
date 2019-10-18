#!/bin/bash

echo "Starting helix tagger"

cd /afs/cern.ch/work/j/jniedzie/private/disapp_tracks/disappTracks/
. setenv.sh

category="3-layers"
sample=$1

echo "Im in `pwd`"
./prepareABCDplots results/abcd_optimization_$category_sample$sample $sample $category
