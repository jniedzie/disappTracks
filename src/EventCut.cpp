//
//  EventCut.cpp
//
//  Created by Jeremi Niedziela on 23/07/2018.
//

#include "EventCut.hpp"

#include <iostream>

using namespace std;

EventCut::EventCut() :
nTracks(range<int>()),
nJets(range<int>()),
nLeptons(range<int>()),
nMuons(range<int>()),
nTaus(range<int>()),
metPt(range<double>()),
metNoMuPt(range<double>()),
jetMetDeltaPhi(range<double>()),
leadingJetPt(range<double>()),
leadingJetEta(range<double>()),
leadingJetChHEF(range<double>()),
leadingJetNeHEF(range<double>()),
metNoMuTrigger(false),
muonsFromZ(false),
metJetPhi(false),
metNoMuJetPhi(false),
muJetR0p4(false),
muTrackR0p4(false),
highJet(false),
tightMuon(false),
twoOpositeMuons(false),
requirePassAllFilters(false)
{

}

EventCut::EventCut(const EventCut &c)
{
  nTracks         = c.nTracks;
  nJets           = c.nJets;
  nLeptons        = c.nLeptons;
  nMuons          = c.nMuons;
  nTaus           = c.nTaus;
  
  metPt           = c.metPt;
  metNoMuPt       = c.metNoMuPt;
  jetMetDeltaPhi  = c.jetMetDeltaPhi;
  
  leadingJetPt    = c.leadingJetPt;
  leadingJetEta   = c.leadingJetEta;
  leadingJetChHEF = c.leadingJetChHEF;
  leadingJetNeHEF = c.leadingJetNeHEF;
  
  metNoMuTrigger  = c.metNoMuTrigger;
  muonsFromZ      = c.muonsFromZ;
  metJetPhi       = c.metJetPhi;
  metNoMuJetPhi   = c.metNoMuJetPhi;
  muJetR0p4       = c.muJetR0p4;
  muTrackR0p4     = c.muTrackR0p4;
  highJet         = c.highJet;
  tightMuon       = c.tightMuon;
  twoOpositeMuons = c.twoOpositeMuons;
}

EventCut::~EventCut()
{
  
}
