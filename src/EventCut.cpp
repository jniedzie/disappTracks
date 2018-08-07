//
//  EventCut.cpp
//
//  Created by Jeremi Niedziela on 23/07/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#include "EventCut.hpp"

#include <iostream>

using namespace std;

EventCut::EventCut(ECut cutType) :
minMetPt(0.0),
minMetNoMuPt(0.0),
minNjets(0),
minNtracks(0),
minNleptons(0),
maxNleptons(999999),
minNmuons(0),
maxNmuons(999999),
maxNtau(999999),
metNoMuTrigger(false),
metJetPhi(false),
metNoMuJetPhi(false),
muJetR0p4(false),
highJetMinPt(0),
highJetMinChHEF(0),
highJetMaxNeHEF(999999),
highJetMaxEta(999999)
{
  if(cutType&kEmpty) return;
  if(cutType&kOneTrack)         minNtracks = 1;
  if(cutType&kOneJet)           minNjets = 1;
  if(cutType&kOneLepton)        minNleptons = 1;
  if(cutType&kOneMuon)          minNmuons = 1;
  if(cutType&kMet100GeV)        minMetPt = 100.0;
  if(cutType&kMetNoMu100GeV)    minMetNoMuPt = 100.0;
  if(cutType&kMet200GeV)        minMetPt = 200.0;
  if(cutType&kMetNoMu200GeV)    minMetNoMuPt = 200.0;
  if(cutType&kMetNoMuTrigger)   metNoMuTrigger = true;
  if(cutType&kNoLepton)         maxNleptons = 0;
  if(cutType&kNoTau)            maxNtau = 0;
  if(cutType&kTwoMuon)          minNmuons = maxNmuons = 2;
  if(cutType&kMuonsFromZ)       muonsFromZ = true;
  if(cutType&kMetJetPhi0p5)     metJetPhi = true;
  if(cutType&kMetNoMuJetPhi0p5) metNoMuJetPhi = true;
  if(cutType&kMuJetR0p4)        muJetR0p4 = true;
  if(cutType&kHighJetPt100GeV)  highJetMinPt = 100.0;
  if(cutType&kHighJetChHEF0p1)  highJetMinChHEF = 0.1;
  if(cutType&kHighJetNeHEF0p8)  highJetMaxNeHEF = 0.8;
  if(cutType&kHighJetEta2p4)    highJetMaxEta = 2.4;
}

EventCut::~EventCut()
{
  
}
