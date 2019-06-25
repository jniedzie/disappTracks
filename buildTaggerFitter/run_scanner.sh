#!/bin/bash

echo "Starting tagger fitter"

export lcgenv=/cvmfs/sft.cern.ch/lcg/releases/lcgenv/0.2-bbd0d/x86_64-centos7-gcc48-dbg/lcgenv

source /cvmfs/sft.cern.ch/lcg/contrib/gcc/4.8/x86_64-centos7-gcc48-opt/setup.sh
source /cvmfs/sft.cern.ch/lcg/releases/pcre/8.34-2c9d9/x86_64-centos7-gcc49-opt/pcre-env.sh
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.16.00/x86_64-centos7-gcc48-opt/bin/thisroot.sh
echo "Which root:`which root`"
echo "Lib path:`echo $LD_LIBRARY_PATH`"

cd /afs/cern.ch/work/j/jniedzie/private/disapp_tracks/disappTracks/buildTaggerFitter

echo "Im in `pwd`"
./taggerScanner
