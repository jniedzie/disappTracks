#!/bin/bash

# add additional modules:
echo "\nAdding additional CMSSW modules\n"

#git cms-addpkg Configuration/StandardSequences
#git cms-addpkg Configuration/DataProcessing
#git cms-addpkg SimG4Core/CustomPhysics
cd src

echo "\nAdding chargino card and other necessary files\n"
mkdir -p Configuration/GenProduction/python/ThirteenTeV/AMSB_chargino/
cd Configuration/GenProduction/python/ThirteenTeV/AMSB_chargino/
wget https://raw.githubusercontent.com/jniedzie/disappTracks/master/eventsGeneration/AMSB_chargino300GeV_ctau10cm_NoFilter_13TeV.py

cd ../..
wget https://raw.githubusercontent.com/jniedzie/disappTracks/master/eventsGeneration/RandomSeed_cfi.py
cd ../../..

echo "\nBuilding CMSSW\n"
scram b -j8
