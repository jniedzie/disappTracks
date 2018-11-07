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
const double fillOpacity = 1.0;
const int fillStyleBack = 1000;
const int fillStyleSignal = 1000;
const int fillStyleData = 1000;

const vector<int> signalMarkers = {
  20, // wino m=300 cτ=3
  20, // wino m=300 cτ=10
  20, // wino m=300 cτ=30
  21, // wino m=500 cτ=10
  21, // wino m=500 cτ=20
  22, // wino m=650 cτ=10
  22, // wino m=650 cτ=20
  23, // wino m=800 cτ=10
  23, // wino m=800 cτ=20
  29, // wino m=1000 cτ=10
  29, // wino m=1000 cτ=20
};

const vector<vector<int>> backColors = {
  {230, 25 , 75 },  // QCD
  {60 , 180, 75 },  // Z->μμ + jets
  {0  , 130, 200},  // tops
  {245, 130, 48 },  // VV
  {145, 30 , 180},  // W->μν + jets
  {70 , 240, 240},  // Z->νν + jets
};

const vector<vector<int>> signalColors = {
  {128, 128, 0  },  // wino m=300 cτ=3
  {170, 110, 40 },  // wino m=300 cτ=10
  {128, 0  , 0  },  // wino m=300 cτ=30
  {170, 110, 40 },  // wino m=500 cτ=10
  {128, 0  , 0  },  // wino m=500 cτ=20
  {170, 110, 40 },  // wino m=650 cτ=10
  {128, 0  , 0  },  // wino m=650 cτ=20
  {170, 110, 40 },  // wino m=800 cτ=10
  {128, 0  , 0  },  // wino m=800 cτ=20
  {170, 110, 40 },  // wino m=1000 cτ=10
  {128, 0  , 0  },  // wino m=1000 cτ=20
};

// {140, 50 , 230} {170, 200, 195} {240, 50 , 100} {255, 225, 25 } {0  , 128, 128} {230, 190, 255} {170, 100, 195} {240, 50 , 230}

const vector<vector<int>> dataColors = {
  {200 , 10, 10},  // single electron (2017B)
};

// Names of background, signal and data samples
const vector<string> backgroundTitle = {
  "QCD",
  "Z#rightarrow#mu#mu + jets",
  "tt",
  "VV",
  "W#rightarrow#mu#nu + jets",
  "Z#rightarrow#nu#nu + jets",
};

const vector<string> signalTitle = {
  "Wino m=300 c#tau=3",
  "Wino m=300 c#tau=10",
  "Wino m=300 c#tau=30",
  "Wino m=500 c#tau=10",
  "Wino m=500 c#tau=20",
  "Wino m=650 c#tau=10",
  "Wino m=650 c#tau=20",
  "Wino m=800 c#tau=10",
  "Wino m=800 c#tau=20",
  "Wino m=1000 c#tau=10",
  "Wino m=1000 c#tau=20",
};

const vector<string> dataTitle = {
  "2017B",
};

// Path to trees with background, signal and data samples (also determines which samples will be merged)
const vector<vector<string>> inFileNameBackground = {
  // QCD
  {
    "../SR_MC/QCD_HT100to200/",
    "../SR_MC/QCD_HT200to300/",
    "../SR_MC/QCD_HT300to500/",
    "../SR_MC/QCD_HT500to700/",
    "../SR_MC/QCD_HT700to1000/",
    "../SR_MC/QCD_HT1000to1500/",
    "../SR_MC/QCD_HT1500to2000/",
    "../SR_MC/QCD_HT2000toInf/",
  },
  // DY + jets
  {
//    "../SR_MC/DYJetsM50_HT100to200/",
    "../SR_MC/DYJetsM50_HT100to200e/",
//    "../SR_MC/DYJetsM50_HT200to400/",
    "../SR_MC/DYJetsM50_HT200to400e/",
//    "../SR_MC/DYJetsM50_HT400to600/",
    "../SR_MC/DYJetsM50_HT400to600e/",
    "../SR_MC/DYJetsM50_HT600to800/",
    "../SR_MC/DYJetsM50_HT800to1200/",
    "../SR_MC/DYJetsM50_HT1200to2500/",
    "../SR_MC/DYJetsM50_HT2500toInf/",
  },
  // tops
  {
    "../SR_MC/TTHad/",
    "../SR_MC/TTLep/",
    "../SR_MC/TTSemi/",
    "../SR_MC/T_tch/",
    "../SR_MC/T_tWch/",
    "../SR_MC/TBar_tch/",
    "../SR_MC/TBar_tWch/",
  },
  // VV
  {
    "../SR_MC/WW/",
    "../SR_MC/WZ/",
    "../SR_MC/ZZ/",
  },
  // W->μν + jets
  {
    "../SR_MC/WJets_HT100to200/",
    "../SR_MC/WJets_HT200to400/",
    "../SR_MC/WJets_HT400to600/",
    "../SR_MC/WJets_HT600to800/",
    "../SR_MC/WJets_HT800to1200/",
    "../SR_MC/WJets_HT1200to2500/",
    "../SR_MC/WJets_HT2500toInf/",
  },
  // Z->νν + jets
  {
    "../SR_MC/ZvvJets_HT100to200/",
    "../SR_MC/ZvvJets_HT200to400/",
    "../SR_MC/ZvvJets_HT400to600/",
    "../SR_MC/ZvvJets_HT600to800/",
    "../SR_MC/ZvvJets_HT800to1200/",
    "../SR_MC/ZvvJets_HT1200to2500/",
    "../SR_MC/ZvvJets_HT2500toInf/",
  },
};

const vector<string> inFileNameSignal = {
  "../Signal/Wino_M_300_cTau_3/",
  "../Signal/Wino_M_300_cTau_10/",
  "../Signal/Wino_M_300_cTau_30/",
  "../Signal/Wino_M_500_cTau_10/",
  "../Signal/Wino_M_500_cTau_20/",
  "../Signal/Wino_M_650_cTau_10/",
  "../Signal/Wino_M_650_cTau_20/",
  "../Signal/Wino_M_800_cTau_10/",
  "../Signal/Wino_M_800_cTau_20/",
  "../Signal/Wino_M_1000_cTau_10/",
  "../Signal/Wino_M_1000_cTau_20/",
};

const vector<string> inFileNameData = {
  "../SR_DATA/MET_Run2017B_31Mar2018/"
};

enum ESignal{
  kWino_M_300_cTau_3,
  kWino_M_300_cTau_10,
  kWino_M_300_cTau_30,
  kWino_M_500_cTau_10,
  kWino_M_500_cTau_20,
  kWino_M_650_cTau_10,
  kWino_M_650_cTau_20,
  kWino_M_800_cTau_10,
  kWino_M_800_cTau_20,
  kWino_M_1000_cTau_10,
  kWino_M_1000_cTau_20,
  kNsignals
};

const double signalCrossSectionTwoTracks[kNsignals] = { // (fb)
  190,  // wino m=300 cτ=3
  190,  // wino m=300 cτ=10
  190,  // wino m=300 cτ=30
  22,   // wino m=500 cτ=10
  22,   // wino m=500 cτ=20
  6.4,  // wino m=650 cτ=10
  6.4,  // wino m=650 cτ=20
  2.2,  // wino m=800 cτ=10
  2.2,  // wino m=800 cτ=20
  0.62, // wino m=1000 cτ=10
  0.62, // wino m=1000 cτ=20
};

const double signalCrossSectionOneTrack[kNsignals] = { // (fb)
  380,  // wino m=300 cτ=3
  380,  // wino m=300 cτ=10
  380,  // wino m=300 cτ=30
  45,   // wino m=500 cτ=10
  45,   // wino m=500 cτ=20
  13,  // wino m=650 cτ=10
  13,  // wino m=650 cτ=20
  4.6,  // wino m=800 cτ=10
  4.6,  // wino m=800 cτ=20
  1.3, // wino m=1000 cτ=10
  1.3, // wino m=1000 cτ=20
};

enum EBackground{
  kQCD,
  kZmumuJets,
  kTT,
  kVV,
  kWmunuJets,
  kZnunuJets,
  kNbackgrounds
};

enum EData{
  kElectron_Run2017B,
  kNdata
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
  kTrackMissingOuterTrackerHits,
  kTrackPixelHits,
  kTrackTrackerHits,
  kTrackRelativeIsolation,
  kTrackAbsoluteIsolation,
  kTrackMetDphi,
  kTrackDedxPerHit,
  
  // per jet variables
  kJetPt,
  kJetEta,
  kJetPhi,
  kJetTrackDr,
  
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

inline int DataColor(EData data){
  return TColor::GetColor(dataColors[data][0],dataColors[data][1],dataColors[data][2]);
}

//,,,{255, 215, 180},{0, 0, 128},{128, 128, 128},{255, 255, 255},{0, 0, 0}

//{2,63,165},{125,135,185},{190,193,212},{214,188,192},{187,119,132},{142,6,59},{74,111,227},{133,149,225},{181,187,227},{230,175,185},{224,123,145},{211,63,106},{17,198,56},{141,213,147},{198,222,199},{234,211,198},{240,185,141},{239,151,8},{15,207,192},{156,222,214},{213,234,231},{243,225,235},{246,196,225},{247,156,212}

//,{157,204,0},{194,0,136},{0,51,128},{255,164,5},{255,168,187},{66,102,0},{255,0,16},{94,241,242},{0,153,143},{224,255,102},{116,10,255},{153,0,0},{255,255,128},{255,255,0},{255,80,5}, {240,163,255}, {0,117,220}, {153,63,0}, {76,0,92}, {25,25,25},{0,92,49}, {43,206,72}, {255,204,153}, {128,128,128}, {148,255,181}, {143,124,0}



#endif /* Helpers_h */
