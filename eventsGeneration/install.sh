#!/bin/bash

# add additional modules:
echo "\nAdding additional CMSSW modules\n"

git cms-addpkg Configuration/StandardSequences
git cms-addpkg Configuration/DataProcessing
git cms-addpkg SimG4Core/CustomPhysics
cd src

echo "\nAdding chargino card and other necessary files\n"
mkdir -p Configuration/GenProduction/python/ThirteenTeV/AMSB_chargino/
mkdir -p DisappTrks/SignalMC/data/

cd Configuration/GenProduction/python/ThirteenTeV/AMSB_chargino/
wget https://raw.githubusercontent.com/jniedzie/disappTracks/master/eventsGeneration/AMSB_chargino300GeV_ctau10cm_NoFilter_13TeV.py

cd ../..
wget https://raw.githubusercontent.com/jniedzie/disappTracks/master/eventsGeneration/RandomSeed_cfi.py
cd ../../..

cd DisappTrks/SignalMC/data/
wget https://github.com/jniedzie/disappTracks/raw/master/eventsGeneration/geant4cards.tar
tar xvf geant4cards.tar
rm geant4cards.tar
cd ../../..

mkdir -p error
mkdir -p output
mkdir -p log
mkdir -p scripts
mkdir -p generatedEvents

echo "\nBuilding CMSSW\n"
scram b -j8
