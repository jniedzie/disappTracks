//
//  TrackProcessor.cpp
//
//  Created by Jeremi Niedziela on 30/01/2019.
//

#include "TrackProcessor.hpp"

TrackProcessor::TrackProcessor()
{
  
}

TrackProcessor::~TrackProcessor()
{
  
}

unique_ptr<Track> TrackProcessor::GetRandomTrack(int nLayers, double maxEta)
{
  auto track = make_unique<Track>();
  
  track->eta = RandDouble(-maxEta, maxEta);
  track->phi = RandDouble(0, 2*TMath::Pi());
  
  track->nTrackerLayers = track->nPixelLayers = nLayers;
  
  double maxTheta = 2*atan(exp(-maxEta));
  
  double minL = layerR[nLayers-1]/sin(maxTheta);
  double maxL = layerR[nLayers]/sin(maxTheta);
  double decayR = RandDouble(minL, maxL);
  
  double theta = track->GetTheta();
  double decayX = decayR*sin(theta)*cos(track->phi);
  double decayY = decayR*sin(theta)*sin(track->phi);
  double decayZ = decayR*cos(theta);
  
  track->decayPoint = make_unique<Point>(decayX, decayY, decayZ);
  
  return track;
}

bool TrackProcessor::IsPassingCut(const shared_ptr<Track> track,
                                  const unique_ptr<TrackCut> &cut)
{
  // check number of hits in pixel, stip and tracker in general
  if(cut->GetRequireSameNpixelHitsLayers()){
    if(track->nPixelHits != track->nPixelLayers) return false;
  }
  
  if(cut->GetRequireSameNtrackerHitsLayers()){
    if(track->nTrackerHits != track->nTrackerLayers) return false;
  }
  
  if(cut->GetRequireMcMatch()){
    if(track->mcMatch == 0) return false;
  }
  
  if(cut->GetNpixelHits().IsOutside(track->nPixelHits))  return false;
  if(cut->GetNpixelLayers().IsOutside(track->nPixelLayers)) return false;
  if(cut->GetNmissingInnerPixel().IsOutside(track->nMissingInnerPixelHits)) return false;
  if(cut->GetNmissingMiddleTracker().IsOutside(track->nMissingMiddleTrackerHits)) return false;
  if(cut->GetNmissingOuterTracker().IsOutside(track->nMissingOuterTrackerHits))  return false;
  
  // check number of dedx, number of detectors, number of clusters
  if(cut->GetNdedxClusters().IsOutside(track->nDedxClusters)) return false;
  if(cut->GetNdetIDs().IsOutside(track->nDetIDs))  return false;
  
  
  if(cut->GetTotalDedx().IsOutside(track->GetTotalDedx())) return false;
  
  for(int iCluster=0;iCluster<track->nDedxClusters;iCluster++){
    if(cut->GetDedxPerCluster().IsOutside(track->dedx[iCluster])) return false;
  }
  
  // check basic kinematical variables
  if(cut->GetPt().IsOutside(track->pt)) return false;
  if(cut->GetEta().IsOutside(track->eta)) return false;
  
  // check calo energy
  if(cut->GetCaloEmEnergy().IsOutside(track->caloEmEnergy))  return false;
  if(cut->GetCaloHadEnergy().IsOutside(track->caloHadEnergy))  return false;
  
  // check isolation
  if(cut->GetRelativeIsolation().IsOutside(track->relativeIsolation))  return false;
  
  // Check track-met ΔΦ
  if(cut->GetTrackMetDeltaPhi().GetMin() > -1000){
    
    // We use TLorentzVector to automatically deal with shifting the angle to [-π,π]
    TLorentzVector metVector, trackVector;
    metVector.SetPtEtaPhiM(track->eventMetPt, track->eventMetEta, track->eventMetPhi, track->eventMetMass);
    trackVector.SetPtEtaPhiM(track->pt, track->eta, track->phi, track->mass);
    
    if(cut->GetTrackMetDeltaPhi().IsOutside(metVector.DeltaPhi(trackVector))) return false;
  }
  
  return true;
}
