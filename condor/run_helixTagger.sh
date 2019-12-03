#!/bin/bash

echo "Starting helix tagger"

cd /afs/cern.ch/work/j/jniedzie/private/disapp_tracks/disappTracks/
. setenv.sh

configPath=""
outputPath=""

if [ $2 -eq 0 ];
then
  configPath="configs/helixTagger.md"
  outputPath="after_L1/all/afterHelixTagging_default"
elif [ $2 -eq 1 ];
then
  configPath="configs/helixTaggerBkg.md"
  outputPath="after_L0/afterHelixTagging_default"
elif [ $2 -eq 2 ];
then
  configPath="configs/helixTaggerHighMerging.md"
  outputPath="after_L1/all/afterHelixTagging_highMerging"
elif [ $2 -eq 3 ];
then
  configPath="configs/helixTaggerHighMergingBkg.md"
  outputPath="after_L0/afterHelixTagging_highMerging"
elif [ $2 -eq 4 ];
then
  configPath="configs/helixTaggerLowSeedChi.md"
  outputPath="after_L1/all/afterHelixTagging_lowSeedChi"
elif [ $2 -eq 5 ];
then
  configPath="configs/helixTaggerLowSeedChiBkg.md"
  outputPath="after_L0/afterHelixTagging_lowSeedChi"
elif [ $2 -eq 6 ];
then
  configPath="configs/helixTaggerLowTrackChi.md"
  outputPath="after_L1/all/afterHelixTagging_lowTrackChi"
elif [ $2 -eq 7 ];
then
  configPath="configs/helixTaggerLowTrackChiBkg.md"
  outputPath="after_L0/afterHelixTagging_lowTrackChi"
elif [ $2 -eq 8 ];
then
  configPath="configs/helixTaggerNoMissing.md"
  outputPath="after_L1/all/afterHelixTagging_noMissing"
elif [ $2 -eq 9 ];
then
  configPath="configs/helixTaggerNoMissingBkg.md"
  outputPath="after_L0/afterHelixTagging_noMissing"
elif [ $2 -eq 10 ];
then
  configPath="configs/helixTaggerNoPionHits.md"
  outputPath="after_L1/all/afterHelixTagging_noPionHits"
else
  echo "Unkown option!!"
fi

echo "Im in `pwd`"
./helixTagger $outputPath/chunk$1 $1 1 $configPath
