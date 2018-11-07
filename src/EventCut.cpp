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

EventCut::~EventCut()
{
  
}
