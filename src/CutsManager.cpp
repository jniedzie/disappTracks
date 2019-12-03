//
//  CutsManager.cpp
//
//  Created by Jeremi Niedziela on 11/04/2019.
//

#include "CutsManager.hpp"

CutsManager::CutsManager()
{
  
}

CutsManager::~CutsManager()
{
  
}

void CutsManager::GetCuts(EventCut &eventCut, TrackCut &trackCut, JetCut &jetCut, LeptonCut &leptonCut)
{
  //----------------------------------------------------------------------------
  // Level 0 cuts
  //----------------------------------------------------------------------------
  
  // Remove bad jets
  jetCut.SetPt(range<double>(30.0, inf));
  jetCut.SetChargedHadronEnergyFraction(range<double>(0.01,0.99));
  jetCut.SetNeutralHadronEnergyFraction(range<double>(0.01,0.99));
  
  // Remove bad tracks
  trackCut.SetNmissingInnerPixel(range<int>(0, 0));
  trackCut.SetNmissingMiddleTracker(range<int>(0, 0));
  trackCut.SetRelativeIsolation(range<double>(0.0, 0.5));
  trackCut.SetNlayers(range<int>(2, inf));
  trackCut.SetEta(range<double>(-2.1, 2.1));
  
  // Check MET properties
  eventCut.SetMetNoMuPt(range<double>(200,inf));
  eventCut.SetRequireMetNoMuTrigger(true);
  eventCut.SetRequirePassingAllFilters(true);
  eventCut.SetJetMetDeltaPhi(range<double>(0.5,inf));
  
  // Check leading jet properties
  eventCut.SetLeadingJetPt(range<double>(100,inf));
  eventCut.SetLeadingJetEta(range<double>(-2.4,2.4));
  eventCut.SetLeadingJetNeHEF(range<double>(0,0.8));
  eventCut.SetLeadingJetChHEF(range<double>(0.1,1.0));
  
  // Check number of objects after cuts
  eventCut.SetNtracks(range<int>(1,inf));
  eventCut.SetNjets(range<int>(1,inf));
//  eventCut.SetNmuons(range<int>(0,0)); // we don't need this one, as nLeptons already contains loose muon cuts
  eventCut.SetNtaus(range<int>(0,0));
  eventCut.SetNleptons(range<int>(0,0));
  
  if(config.params["cuts_level"] == 0) return;
  
  //----------------------------------------------------------------------------
  // Level 1 cuts
  //----------------------------------------------------------------------------
 
  
  trackCut.SetCaloEmEnergy(range<double>(0.0, 2.0));
//  trackCut.SetVetoCracks(true);
  
  if(config.category == "2-tracks"){
    eventCut.SetNtracks(range<int>(2,2));
    trackCut.SetNlayers(range<int>(2,6));
  }
  else if(config.category == "3-layers"){
    eventCut.SetNtracks(range<int>(1,1));
    trackCut.SetNlayers(range<int>(3, 3));
  }
  else if(config.category == "4-layers"){
    eventCut.SetNtracks(range<int>(1,1));
    trackCut.SetNlayers(range<int>(4, 4));
  }
  else if(config.category == "5-6-layers"){
    eventCut.SetNtracks(range<int>(1,1));
    trackCut.SetNlayers(range<int>(5, 6));
  }
  else if(config.category == "all"){
    eventCut.SetNtracks(range<int>(1, 2));
    trackCut.SetNlayers(range<int>(3, 6));
  }
  
  if(config.params["cuts_level"] == 1) return;
  
  //----------------------------------------------------------------------------
  // Level 2 cuts (gen info)
  //----------------------------------------------------------------------------
  
//  eventCut.SetNgenPions(range<int>(1));
//  eventCut.SetGetPionPt(range<double>(500, inf));
  eventCut.SetMetPt(range<double>(500, inf));
  trackCut.SetEta(range<double>(0, 1.2));
  trackCut.SetNlayers(range<int>(6, inf));
  
//  trackCut.SetRequireCorrectCharge(true);
//  trackCut.SetRequireCorrectNlayers(true);
  
  if(config.params["cuts_level"] == 2) return;
  
  //
  
  return;
}

void CutsManager::GetZmumuCuts(EventCut &eventCut, TrackCut &trackCut, JetCut &jetCut, LeptonCut &leptonCut)
{
  GetCuts(eventCut, trackCut, jetCut, leptonCut);
  eventCut.SetNmuons(range<int>(2,2));
  eventCut.SetNleptons(range<int>(2,inf));
  eventCut.SetRequireMuonsFromZ(true);
  eventCut.SetRequireTwoOppositeMuons(true);
}

void CutsManager::GetWmunuCuts(EventCut &eventCut, TrackCut &trackCut, JetCut &jetCut, LeptonCut &leptonCut)
{
  GetCuts(eventCut, trackCut, jetCut, leptonCut);
  eventCut.SetNmuons(range<int>(1,1));
  eventCut.SetNleptons(range<int>(1,inf));
}

void CutsManager::GetLowMetCuts(EventCut &eventCut, TrackCut &trackCut, JetCut &jetCut, LeptonCut &leptonCut)
{
  GetCuts(eventCut, trackCut, jetCut, leptonCut);
  eventCut.SetMetNoMuPt(range<double>(150, inf));
}
