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
  jetCut.SetChargedHadronEnergyFraction(range<double>(0.01,0.99));
  jetCut.SetNeutralHadronEnergyFraction(range<double>(0.01,0.99));
  jetCut.SetPt(range<double>(30.0, inf));
  
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
  eventCut.SetRequireHighJet(true);
//  eventCut.SetNmuons(range<int>(0,0));
  eventCut.SetNtaus(range<int>(0,0));
  eventCut.SetNleptons(range<int>(0,0));
  
  if(config->performCutsLevel == 0) return;
  
  //----------------------------------------------------------------------------
  // Level 1 cuts
  //----------------------------------------------------------------------------
 
  if(config->category == "2-tracks"){
    eventCut.SetNtracks(range<int>(2,2));
    trackCut.SetNlayers(range<int>(2,10));
    
    trackCut.SetCaloEmEnergy(range<double>(0.0,2.0));
  }
  else if(config->category == "3-layers"){
    eventCut.SetNtracks(range<int>(1,1));
    trackCut.SetNlayers(range<int>(3, 3));
    
    trackCut.SetCaloEmEnergy(range<double>(0.0,2.0));
    trackCut.SetTrackMetDeltaPhi(range<double>(-TMath::Pi()/2,TMath::Pi()/2));
    trackCut.SetDedxPerCluster(range<double>(3.0,inf));
  }
  else if(config->category == "4-layers"){
    eventCut.SetNtracks(range<int>(1,1));
    trackCut.SetNlayers(range<int>(4, 4));
    
    trackCut.SetCaloEmEnergy(range<double>(0.0,2.0));
    trackCut.SetTrackMetDeltaPhi(range<double>(-TMath::Pi()/2,TMath::Pi()/2));
    trackCut.SetDedxPerCluster(range<double>(3.0,inf));
  }
  
  if(config->performCutsLevel == 1) return;
  
  //
  
  return;
}
