//
//  Helpers.hpp
//
//  Created by Jeremi Niedziela on 13/06/2018.
//

#ifndef Helpers_h
#define Helpers_h

#include "Config.hpp"

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TColor.h>

#include <vector>
#include <iostream>
#include <map>
#include <numeric>

using namespace std;

// Plotting style
const double fillOpacity = 0.1;
const int fillStyleBack = 1000;
const int fillStyleSignal = 3003;

const vector<int> signalMarkers = {
  20, // wino m=300 cτ=10
  21, // wino m=300 cτ=3
  22, // wino m=300 cτ=30
  23, // wino m=500 cτ=10
};

const vector<vector<int>> backColors = {
  {230, 25 , 75 },  // Z->μμ + jets
  {60 , 180, 75 },  // tops
  {0  , 130, 200},  // VV
  {245, 130, 48 },  // W->μν + jets
  {145, 30 , 180},  // QCD
};
//  {255, 225, 25},{70, 240, 240},{240, 50, 230},{250, 190, 190},{0, 128, 128},{230, 190, 255}

const vector<vector<int>> signalColors = {
  {170, 110, 40 },  // wino m=300 cτ=10
  {128, 128, 0  },  // wino m=300 cτ=3
  {128, 0  , 0  },  // wino m=300 cτ=30
  {170, 100, 195},  // wino m=500 cτ=10
};

// Names of background, signal and data samples
const vector<string> backgroundTitle = {
  "Z#rightarrow#mu#mu + jets",
  "tt",
  "VV",
  "W#rightarrow#mu#nu + jets",
  "QCD",
};

const vector<string> signalTitle = {
  "Wino_M_300_cTau_10",
  "Wino_M_300_cTau_3",
  "Wino_M_300_cTau_30",
  "Wino_M_500_cTau_10"
};

const vector<string> dataTitle = {
  "Single electron (2017B)",
};

// Path to trees with background, signal and data samples (also determines which samples will be merged)
const vector<vector<string>> inFileNameBackground = {
  // DY + jets
  {
    "../SR_MC/DYJetsM50_HT100to200/tree.root",
    "../SR_MC/DYJetsM50_HT100to200e/tree.root",
    "../SR_MC/DYJetsM50_HT200to400/tree.root",
    "../SR_MC/DYJetsM50_HT200to400e/tree.root",
    "../SR_MC/DYJetsM50_HT400to600/tree.root",
    "../SR_MC/DYJetsM50_HT400to600e/tree.root",
    "../SR_MC/DYJetsM50_HT600to800/tree.root",
    "../SR_MC/DYJetsM50_HT800to1200/tree.root",
    "../SR_MC/DYJetsM50_HT1200to2500/tree.root",
    "../SR_MC/DYJetsM50_HT2500toInf/tree.root",
  },
  // tops
  {
    "../SR_MC/TTHad/tree.root",
    "../SR_MC/TTLep/tree.root",
    "../SR_MC/TTSemi/tree.root",
    "../SR_MC/T_tch/tree.root",
    "../SR_MC/T_tWch/tree.root",
    "../SR_MC/TBar_tch/tree.root",
    "../SR_MC/TBar_tWch/tree.root",
  },
  // VV
  {
    "../SR_MC/WW/treeProducerXtracks/tree.root",
    "../SR_MC/WZ/treeProducerXtracks/tree.root",
    "../SR_MC/ZZ/treeProducerXtracks/tree.root",
  },
  // W->μν + jets
  {
    "../SR_MC/WJets_HT100to200/tree.root",
    "../SR_MC/WJets_HT200to400/tree.root",
    "../SR_MC/WJets_HT400to600/tree.root",
    "../SR_MC/WJets_HT600to800/tree.root",
    "../SR_MC/WJets_HT800to1200/tree.root",
    "../SR_MC/WJets_HT1200to2500/tree.root",
    "../SR_MC/WJets_HT2500toInf/tree.root",
  },
  // QCD
  {
    "../SR_MC/QCD_HT100to200/tree.root",
    "../SR_MC/QCD_HT200to300/tree.root",
    "../SR_MC/QCD_HT300to500/tree.root",
    "../SR_MC/QCD_HT500to700/tree.root",
    "../SR_MC/QCD_HT700to1000/tree.root",
    "../SR_MC/QCD_HT1000to1500/tree.root",
    "../SR_MC/QCD_HT1500to2000/tree.root",
    "../SR_MC/QCD_HT2000toInf/tree.root",
  },
  
};

const vector<string> inFileNameSignal = {
  "../adish/Signal/Wino_M_300_cTau_10/treeProducerXtracks/tree.root",
  "../adish/Signal/Wino_M_300_cTau_3/treeProducerXtracks/tree.root",
  "../adish/Signal/Wino_M_300_cTau_30/treeProducerXtracks/tree.root",
  "../adish/Signal/Wino_M_500_cTau_10/treeProducerXtracks/tree.root"
};

const vector<string> inFileNameData = {
  "../adish/Data/SingleElectron_Run2017B_17Nov2017/treeProducerXtracks/tree.root"
};

// Constants for tracker layers
const int nLayers = 14;
const double layerR[nLayers] = { 29, 68, 109, 160, 250, 340, 430, 520, 610, 696, 782, 868, 965, 1080 };

enum EVar{
  kCustom,
  
  // per event variables
  kNvertices,
  kNisoTracks,
  kNjets,
  kNjets30,
  kNjets30a,
  kMetSumEt,
  kMetPt,
  kMetMass,
  kMetEta,
  kMetPhi,
  kMetJetDphi,
  
  // per track variables
  kTrackNclusters,
  kTrackTotalDedx,
  kTrackDedxPerCluster,
  kTrackPt,
  kTrackEta,
  kTrackPhi,
  kTrackCaloEm,
  kTrackCaloHad,
  kTrackDxy,
  kTrackDz,
  kTrackCharge,
  kTrackMass,
  kTrackPid,
  
  // per jet variables
  kJetPt,
  kJetEta,
  kJetPhi,
  
  // per track per layer variables
  kDedx,  ///< dE/dx per layer
  kSizeX, ///< X cluster size in each layer
  kSizeY  ///< Y cluster size in each layer
};

inline int BackColor(EBackground bck){
  return TColor::GetColor(backColors[bck][0],backColors[bck][1],backColors[bck][2]);
}

inline int SignalColor(ESignal sig){
  return TColor::GetColor(signalColors[sig][0],signalColors[sig][1],signalColors[sig][2]);
}

//,,,{255, 215, 180},{0, 0, 128},{128, 128, 128},{255, 255, 255},{0, 0, 0}

//{2,63,165},{125,135,185},{190,193,212},{214,188,192},{187,119,132},{142,6,59},{74,111,227},{133,149,225},{181,187,227},{230,175,185},{224,123,145},{211,63,106},{17,198,56},{141,213,147},{198,222,199},{234,211,198},{240,185,141},{239,151,8},{15,207,192},{156,222,214},{213,234,231},{243,225,235},{246,196,225},{247,156,212}

//,{157,204,0},{194,0,136},{0,51,128},{255,164,5},{255,168,187},{66,102,0},{255,0,16},{94,241,242},{0,153,143},{224,255,102},{116,10,255},{153,0,0},{255,255,128},{255,255,0},{255,80,5}, {240,163,255}, {0,117,220}, {153,63,0}, {76,0,92}, {25,25,25},{0,92,49}, {43,206,72}, {255,204,153}, {128,128,128}, {148,255,181}, {143,124,0}



#endif /* Helpers_h */
