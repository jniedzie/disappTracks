//
//  EventCut.cpp
//
//  Created by Jeremi Niedziela on 23/07/2018.
//

#include "EventCut.hpp"

#include <iostream>

using namespace std;

EventCut::EventCut(ECut cutType) :
minNjets(0),
minNtracks(0),
maxNtracks(999999),
minNleptons(0),
maxNleptons(999999),
minNmuons(0),
maxNmuons(999999),
maxNtau(999999),
metNoMuTrigger(false),
muonsFromZ(false),
metJetPhi(false),
metNoMuJetPhi(false),
muJetR0p4(false),
muTrackR0p4(false),
highJet(false),
tightMuon(false),
twoOpositeMuons(false),
minMetPt(0.0),
minMetNoMuPt(0.0),
highJetMinPt(0.0),
highJetMinChHEF(0.0),
highJetMaxNeHEF(999999.9),
highJetMaxEta(999999.9),
minJetMetPhi(0),
requirePassAllFilters(false)
{
  if(cutType&kEmpty)            return;
  if(cutType&kOneTrack)         minNtracks = 1;
  if(cutType&kOneJet)           minNjets = 1;
  if(cutType&kOneLepton)        minNleptons = 1;
  if(cutType&kOneMuon)          minNmuons = 1;
  if(cutType&kNoLepton)         maxNleptons = 0;
  if(cutType&kNoTau)            maxNtau = 0;
  if(cutType&kTwoMuon)          minNmuons = maxNmuons = 2;
  if(cutType&kMetNoMuTrigger)   metNoMuTrigger = true;
  if(cutType&kMuonsFromZ)       muonsFromZ = true;
  if(cutType&kMetJetPhi0p5)     minJetMetPhi = 0.5;
  if(cutType&kMetNoMuJetPhi0p5) metNoMuJetPhi = true;
  if(cutType&kMuJetR0p4)        muJetR0p4 = true;
  if(cutType&kMuTrackR0p4)      muTrackR0p4 = true;
  if(cutType&kHighJet)          highJet = true;
  if(cutType&kTightMuon)        tightMuon = true;
  if(cutType&kTwoOppositeMuons) twoOpositeMuons = true;
  if(cutType&kMet100GeV)        minMetPt = 100.0;
  if(cutType&kMetNoMu100GeV)    minMetNoMuPt = 100.0;
  if(cutType&kMet200GeV)        minMetPt = 200.0;
  if(cutType&kMetNoMu200GeV)    minMetNoMuPt = 200.0;
  if(cutType&kHighJetPt100GeV)  highJetMinPt = 100.0;
  if(cutType&kHighJetChHEF0p1)  highJetMinChHEF = 0.1;
  if(cutType&kHighJetNeHEF0p8)  highJetMaxNeHEF = 0.8;
  if(cutType&kHighJetEta2p4)    highJetMaxEta = 2.4;
}

EventCut::~EventCut()
{
  
}
