//
//  Track.cpp
//
//  Created by Jeremi Niedziela on 16/07/2018.
//

#include "Track.hpp"

Track::Track() :
pt(inf),
eta(inf),
phi(inf),
mass(inf),
caloEmEnergy(-1.0),
caloHadEnergy(-1.0),
relativeIsolation(inf),
dxy(inf),
dxyErr(inf),
dz(inf),
dzErr(inf),
charge(inf),
pid(inf),
mcMatch(inf),
nTrackerLayers(inf),
nPixelLayers(inf),
nTrackerHits(inf),
nPixelHits(inf),
nMissingInnerPixelHits(inf),
nMissingOuterPixelHits(inf),
nMissingInnerStripHits(inf),
nMissingOuterStripHits(inf),
nMissingInnerTrackerHits(inf),
nMissingOuterTrackerHits(inf),
nMissingMiddleTrackerHits(inf),
nDetIDs(-1),
nDedxClusters(-1),
eventMetPt(inf),
eventMetEta(inf),
eventMetPhi(inf),
eventMetMass(inf)
{
  for(int iLayer=0;iLayer<nLayers;iLayer++){
    dedx.push_back(0.0);
    subDetId.push_back(-1);
    sizeX.push_back(-1);
    sizeY.push_back(-1);
  }
};

double Track::GetDedxInSubDet(int det)
{
  double dedxSum=0;
  
  for(int i=0;i<nLayers;i++){
    if(subDetId[i] == det){
      dedxSum += dedx[i];
    }
  }
  return dedxSum;
}

void Track::Print()
{
  cout<<"PID:"<<pid<<"\trel iso:"<<relativeIsolation<<endl;
  cout<<"eta:"<<eta<<"\tphi:"<<phi<<"\tpT:"<<pt<<endl;
  cout<<"Tracker layers:"<<nTrackerLayers<<"\tpixel layers:"<<nPixelLayers<<endl;
  cout<<"Missing outer tracker hits:"<<nMissingOuterTrackerHits<<endl;
  
//  for(int iLayer=0;iLayer<nLayers;iLayer++){
//    cout<<"Layer:"<<iLayer<<"\tsub-det ID:"<<subDetId[iLayer]<<"\tdEdx:"<<dedx[iLayer]<<endl;
//  }
}

bool Track::IsPassingCut(const unique_ptr<TrackCut> &cut)
{
  // check number of hits in pixel, stip and tracker in general
  if(cut->GetRequireSameNpixelHitsLayers()){
    if(nPixelHits != nPixelLayers) return false;
  }
  
  if(cut->GetRequireSameNtrackerHitsLayers()){
    if(nTrackerHits != nTrackerLayers) return false;
  }
  
  if(cut->GetRequireMcMatch()){
    if(mcMatch == 0) return false;
  }
  
  if(cut->GetNpixelHits().IsOutside(nPixelHits))  return false;
  if(cut->GetNpixelLayers().IsOutside(nPixelLayers)) return false;
  if(cut->GetNmissingInnerPixel().IsOutside(nMissingInnerPixelHits)) return false;
  if(cut->GetNmissingMiddleTracker().IsOutside(nMissingMiddleTrackerHits)) return false;
  if(cut->GetNmissingOuterTracker().IsOutside(nMissingOuterTrackerHits))  return false;
  
  // check number of dedx, number of detectors, number of clusters
  if(cut->GetNdedxClusters().IsOutside(nDedxClusters)) return false;
  if(cut->GetNdetIDs().IsOutside(nDetIDs))  return false;
  
  
  if(cut->GetTotalDedx().IsOutside(GetTotalDedx())) return false;
  
  for(int iCluster=0;iCluster<nDedxClusters;iCluster++){
    if(cut->GetDedxPerCluster().IsOutside(dedx[iCluster])) return false;
  }
  
  // check basic kinematical variables
  if(cut->GetPt().IsOutside(pt)) return false;
  if(cut->GetEta().IsOutside(eta)) return false;
  
  // check calo energy
  if(cut->GetCaloEmEnergy().IsOutside(caloEmEnergy))  return false;
  if(cut->GetCaloHadEnergy().IsOutside(caloHadEnergy))  return false;
  
  // check isolation
  if(cut->GetRelativeIsolation().IsOutside(relativeIsolation))  return false;
  
  // Check track-met ΔΦ
  
  if(cut->GetTrackMetDeltaPhi().GetMin() > -1000){
    
    // We use TLorentzVector to automatically deal with shifting the angle to [-π,π]
    TLorentzVector metVector, trackVector;
    metVector.SetPtEtaPhiM(eventMetPt, eventMetEta, eventMetPhi, eventMetMass);
    trackVector.SetPtEtaPhiM(pt, eta, phi, mass);
    
    if(cut->GetTrackMetDeltaPhi().IsOutside(metVector.DeltaPhi(trackVector))) return false;
  }
  
  return true;
}

void Track::SetDeDxInLayer(int layer, float value)
{
  dedx[layer] = value;
  CalculateInternals();
}

void Track::SetSubDetIdInLayer(int layer, int id)
{
  subDetId[layer] = id;
  CalculateInternals();
}

void Track::CalculateInternals()
{
  set<int> uniqueDets;
  for(int d : subDetId){
    uniqueDets.insert(d);
  }
  nDetIDs = (int)uniqueDets.size();
  
  nDedxClusters=0;
  for(float d : dedx){
    if(d > 0.000001) nDedxClusters++;
  }
}
