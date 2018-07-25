//
//  Helpers.hpp
//
//  Created by Jeremi Niedziela on 13/06/2018.
//

#ifndef Helpers_h
#define Helpers_h

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

const bool analyzeData = false;

const int nLayers = 14;
const double layerR[nLayers] = { 29, 68, 109, 160, 250, 340, 430, 520, 610, 696, 782, 868, 965, 1080 };


const double fillOpacity = 0.1;

inline TLegend* GetLegend(double legendW = 0.15, double legendH = 0.5, double legendX = 0.75, double legendY = 0.25,
                          const char* header="")
{
  
  TLegend *leg = new TLegend(legendX,legendY,legendX+legendW,legendY+legendH);
  leg->SetHeader(header);
  return leg;
}

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

enum EBackground{
  kDYJetsToLL_M50,
  kTBar_tWch_noFullyHad,
  kTBar_tch,
  kTTLep_pow,
  kTTSemi_pow,
  kT_tWch_noFullyHad,
  kT_tch,
  kWJetsToLNu_LO,
  kWW,
  kWZ,
  kZZ,
  kNbackgrounds
};

static const char* backgroundTitle[kNbackgrounds] = {
  "DYJetsToLL_M50",
  "TBar_tWch_noFullyHad",
  "TBar_tch",
  "TTLep_pow",
  "TTSemi_pow",
  "T_tWch_noFullyHad",
  "T_tch",
  "WJetsToLNu_LO",
  "WW",
  "WZ",
  "ZZ",
};

const string inFileNameBackground[kNbackgrounds] = {
  "../adish/Background/DYJetsToLL_M50/tree.root",
  "../adish/Background/TBar_tWch_noFullyHad/tree.root",
  "../adish/Background/TBar_tch/tree.root",
  "../adish/Background/TTLep_pow/tree.root",
  "../adish/Background/TTSemi_pow/tree.root",
  "../adish/Background/T_tWch_noFullyHad/tree.root",
  "../adish/Background/T_tch/tree.root",
  "../adish/Background/WJetsToLNu_LO/tree.root",
  "../adish/Background/WW/tree.root",
  "../adish/Background/WZ/tree.root",
  "../adish/Background/ZZ/tree.root"
};


const int backColors[kNbackgrounds][3] = {
  {240,163,255},
  {0,117,220},
  {153,63,0},
  {76,0,92},
  {25,25,25},
  {0,92,49},
  {43,206,72},
  {255,204,153},
  {128,128,128},
  {148,255,181},
  {143,124,0}
};

inline int BackColor(EBackground bck){
  return TColor::GetColor(backColors[bck][0],backColors[bck][1],backColors[bck][2]);
}

//,{157,204,0},{194,0,136},{0,51,128},{255,164,5},{255,168,187},{66,102,0},{255,0,16},{94,241,242},{0,153,143},{224,255,102},{116,10,255},{153,0,0},{255,255,128},{255,255,0},{255,80,5}

const string inFileNameSignal = "../adish/Signal/tree.root";
const string inFileNameData = "../adish/Data/tree.root";

#endif /* Helpers_h */
