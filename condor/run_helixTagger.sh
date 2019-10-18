#!/bin/bash

echo "Starting helix tagger"

cd /afs/cern.ch/work/j/jniedzie/private/disapp_tracks/disappTracks/
. setenv.sh

echo "Im in `pwd`"
./helixTagger afterHelixTagging/chunk$1 $1 1
