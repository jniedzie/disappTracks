#!/bin/bash

echo "Starting helix tagger"

cd /afs/cern.ch/work/j/jniedzie/private/disapp_tracks/disappTracks/
. setenv.sh

configPath="configs/helixTagger.md"
outputPath="afterHelixTagging"

if [ $2 -eq 0 ];
then
  configPath="configs/helixTagger.md"
  outputPath="afterHelixTagging"
elif [ $2 -eq 1 ];
then
  configPath="configs/helixTaggerWithEndcaps.md"
  outputPath="afterHelixTagging_withEndcaps"
elif [ $2 -eq 2 ];
then
  configPath="configs/helixTaggerAsymmetric.md"
  outputPath="afterHelixTagging_asymmetric"
elif [ $2 -eq 3 ];
then
  configPath="configs/helixTaggerHighMerging.md"
  outputPath="afterHelixTagging_highMerging"
elif [ $2 -eq 4 ];
then
  configPath="configs/helixTaggerLowSeedChi.md"
  outputPath="afterHelixTagging_lowSeedChi"
elif [ $2 -eq 5 ];
then
  configPath="configs/helixTaggerLowTrackChi.md"
  outputPath="afterHelixTagging_lowTrackChi"
elif [ $2 -eq 6 ];
then
  configPath="configs/helixTaggerNoMissing.md"
  outputPath="afterHelixTagging_noMissing"
elif [ $2 -eq 7 ];
then
  configPath="configs/helixTaggerNoNoise.md"
  outputPath="afterHelixTagging_noNoise"
elif [ $2 -eq 8 ];
then
  configPath="configs/helixTaggerNoNoiseBkg.md"
  outputPath="afterHelixTagging_noNoise"
elif [ $2 -eq 9 ];
then
  configPath="configs/helixTaggerLowSeedChiBkg.md"
  outputPath="afterHelixTagging_lowSeedChi"
elif [ $2 -eq 10 ];
then
  configPath="configs/helixTaggerLowTrackChiBkg.md"
  outputPath="afterHelixTagging_lowTrackChi"
fi

echo "Im in `pwd`"
./helixTagger $outputPath/chunk$1 $1 1 $configPath
