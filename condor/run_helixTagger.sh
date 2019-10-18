#!/bin/bash

echo "Starting helix tagger"

#export lcgenv=/cvmfs/sft.cern.ch/lcg/releases/lcgenv/0.2-bbd0d/x86_64-centos7-gcc48-dbg/lcgenv

# setup GCC
#source /cvmfs/sft.cern.ch/lcg/contrib/gcc/4.8/x86_64-centos7-gcc48-opt/setup.sh
. /cvmfs/sft.cern.ch/lcg/contrib/gcc/8/x86_64-centos7/setup.sh

# setup PCRE
source /cvmfs/sft.cern.ch/lcg/releases/pcre/8.34-2c9d9/x86_64-centos7-gcc49-opt/pcre-env.sh

# setup ROOT
#source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.16.00/x86_64-centos7-gcc48-opt/bin/thisroot.sh
. /cvmfs/sft.cern.ch/lcg/releases/LCG_96/ROOT/6.18.00/x86_64-centos7-gcc8-opt/bin/thisroot.sh

cd /afs/cern.ch/work/j/jniedzie/private/disapp_tracks/disappTracks/

echo "Im in `pwd`"
./helixTagger $1 $2
