//  EventCut.cpp
//
//  Created by Jeremi Niedziela on 23/07/2018.

#include "EventCut.hpp"

EventCut::EventCut() :

nTracks(range<int>()),
nJets(range<int>()),
nLeptons(range<int>()),
nMuons(range<int>()),
nTaus(range<int>()),

metPt(range<double>()),
metNoMuPt(range<double>()),
jetMetDeltaPhi(range<double>()),
jetMuonDeltaPhi(range<double>()),
trackMuonDeltaPhi(range<double>()),

leadingJetPt(range<double>()),
leadingJetEta(range<double>()),
leadingJetChHEF(range<double>()),
leadingJetNeHEF(range<double>()),

requiresMetNoMuTrigger(false),
requiresMuonsFromZ(false),
requiresTightMuon(false),
requiresTwoOpositeMuons(false),
requiresPassingAllFilters(false)
{

}

EventCut::EventCut(const EventCut &c)
{
  nTracks                 = c.nTracks;
  nJets                   = c.nJets;
  nLeptons                = c.nLeptons;
  nMuons                  = c.nMuons;
  nTaus                   = c.nTaus;
  
  metPt                   = c.metPt;
  metNoMuPt               = c.metNoMuPt;
  jetMetDeltaPhi          = c.jetMetDeltaPhi;
  jetMuonDeltaPhi         = c.jetMuonDeltaPhi;
  trackMuonDeltaPhi       = c.trackMuonDeltaPhi;
  
  leadingJetPt            = c.leadingJetPt;
  leadingJetEta           = c.leadingJetEta;
  leadingJetChHEF         = c.leadingJetChHEF;
  leadingJetNeHEF         = c.leadingJetNeHEF;
  
  requiresMetNoMuTrigger  = c.requiresMetNoMuTrigger;
  requiresMuonsFromZ      = c.requiresMuonsFromZ;
  requiresTightMuon       = c.requiresTightMuon;
  requiresTwoOpositeMuons = c.requiresTwoOpositeMuons;
  requiresPassingAllFilters   = c.requiresPassingAllFilters;
}

EventCut::~EventCut()
{
  
}
