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
nMissingMiddleTrackerHits(99999)
{
  for(int iLayer=0;iLayer<nLayers;iLayer++){
    dedx.push_back(0.0);
    subDetId.push_back(-1);
    sizeX.push_back(-1);
    sizeY.push_back(-1);
  }
};


int Track::GetNclusters()
{
  int nClusters=0;
  for(float d : dedx){
    if(d > 0.000001) nClusters++;
  }
  return nClusters;
}

int Track::GetNdetIDs()
{
  set<int> uniqueDets;
  for(int d : subDetId){
    uniqueDets.insert(d);
  }
  return (int)uniqueDets.size();
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
  for(int iLayer=0;iLayer<nLayers;iLayer++){
    cout<<"Layer:"<<iLayer<<"\tsub-det ID:"<<subDetId[iLayer]<<"\tdEdx:"<<dedx[iLayer]<<endl;
  }
}

bool Track::IsPassingCut(TrackCut *cut)
{
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
  
  return true;
}



