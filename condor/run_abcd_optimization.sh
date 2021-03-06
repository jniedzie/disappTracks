#!/bin/bash

# arguments: sampleIndex category nDedxBins nMetBins

echo "Starting ABCD optimization"

cd /afs/cern.ch/work/j/jniedzie/private/disapp_tracks/disappTracks/
. setenv.sh

category="3-layers"
if [ $2 -eq 0 ]
then
  category="3-layers"
elif [ $2 -eq 1 ]
then
  category="4-layers"
elif [ $2 -eq 2 ]
then
  category="5-6-layers"
fi

sample=$1
nDedxBins=$3
nMetBins=$4

echo "Im in `pwd`"
output_path=results/abcd_optimization_${nDedxBins}x${nMetBins}_${category}_sample${sample}.txt

echo "Output path:$output_path"
./prepareABCDplots $output_path $sample $category $nDedxBins $nMetBins
