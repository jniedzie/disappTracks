//
//  Track.cpp
//
//  Created by Jeremi Niedziela on 16/07/2018.
//

#include "Track.hpp"

Track::Track() :
eta(99999),
phi(99999),
caloEmEnergy(-1.0),
caloHadEnergy(-1.0),
dxy(99999),
dxyErr(99999),
dz(99999),
dzErr(99999),
charge(99999),
mass(99999),
pt(99999),
pid(99999),
relIso03(99999),
nTrackerLayers(99999),
nPixelLayers(99999),
nTrackerHits(99999),
nPixelHits(99999),
nMissingInnerPixelHits(99999),
nMissingOuterPixelHits(99999),
nMissingInnerStripHits(99999),
nMissingOuterStripHits(99999),
nMissingInnerTrackerHits(99999),
nMissingOuterTrackerHits(99999),
nMissingMiddleTrackerHits(99999),
nDetIDs(-1),
nClusters(-1)
{
  for(int iLayer=0;iLayer<nLayers;iLayer++){
    dedx.push_back(0.0);
    subDetId.push_back(-1);
    sizeX.push_back(-1);
    sizeY.push_back(-1);
  }
};


void Track::CalculateInternals()
{
  set<int> uniqueDets;
  for(int d : subDetId){
    uniqueDets.insert(d);
  }
  nDetIDs = (int)uniqueDets.size();
  
  nClusters=0;
  for(float d : dedx){
    if(d > 0.000001) nClusters++;
  }
}

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
  cout<<"PID:"<<pid<<"\trel iso:"<<relIso03<<endl;
  cout<<"eta:"<<eta<<"\tphi:"<<phi<<"\tpT:"<<pt<<endl;
  cout<<"Tracker layers:"<<nTrackerLayers<<"\tpixel layers:"<<nPixelLayers<<endl;
  cout<<"Missing outer tracker hits:"<<nMissingOuterTrackerHits<<endl;
  
//  for(int iLayer=0;iLayer<nLayers;iLayer++){
//    cout<<"Layer:"<<iLayer<<"\tsub-det ID:"<<subDetId[iLayer]<<"\tdEdx:"<<dedx[iLayer]<<endl;
//  }
}

bool Track::IsPassingCut(TrackCut *cut)
{
  // check number of hits in pixel, stip and tracker in general
  if(cut->GetRequireSameNpixelHitsLayers()){
    if(nPixelHits != nPixelLayers) return false;
  }
  
  if(cut->GetRequireSameNtrackerHitsLayers()){
    if(nTrackerHits != nTrackerLayers) return false;
  }
  
  if(   nPixelHits < cut->GetMinNpixelHits()
     || nPixelHits > cut->GetMaxNpixelHits()) return false;
  
  if(   nPixelLayers < cut->GetMinNpixelLayers()
     || nPixelLayers > cut->GetMaxNpixelLayers()) return false;
  
  if(   nMissingInnerPixelHits < cut->GetMinNmissingInnerPixel()
     || nMissingInnerPixelHits > cut->GetMaxNmissingInnerPixel()) return false;
  
  if(   nMissingMiddleTrackerHits < cut->GetMinNmissingMiddleTracker()
     || nMissingMiddleTrackerHits > cut->GetMaxNmissingMiddleTracker()) return false;
  
  if(   nMissingOuterTrackerHits < cut->GetMinNmissingOuterTracker()
     || nMissingOuterTrackerHits > cut->GetMaxNmissingOuterTracker()) return false;
  
  // check number of dedx clusters
  if(GetNclusters() < cut->GetMinDedxClusters() || GetNclusters() > cut->GetMaxDedxClusters()){
    return false;
  }
  
  if(GetNdetIDs() < cut->GetMinDets() || GetNdetIDs() > cut->GetMaxDets()){
    return false;
  }
  
  // check values of dedx along the track
  if(GetTotalDedx() < cut->GetMinTotalDedx() || GetTotalDedx() > cut->GetMaxTotalDedx()){
    return false;
  }
  
  for(int iCluster=0;iCluster<GetNclusters();iCluster++){
    if(dedx[iCluster] < cut->GetMinDedxPerCluster()) return false;
  }
  
  // check pt
  if(pt < cut->GetMinPt()) return false;
  
  // check calo energy
  if(caloEmEnergy > cut->GetMaxEmCalo()) return false;
  if(caloHadEnergy > cut->GetMaxHadCalo()) return false;
  
  // check eta
  if(fabs(eta) > cut->GetMaxEta()) return false;
  
  // check isolation
  if(relIso03 > cut->GetMaxRelativeIsolation()) return false;
  
  return true;
}



