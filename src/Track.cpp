//
//  Track.cpp
//  disappTracksTarget
//
//  Created by Jeremi Niedziela on 16/07/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
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
pid(99999)
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



