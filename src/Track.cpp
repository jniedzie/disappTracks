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
    detType.push_back(-1);
    layer.push_back(-1);
    ladder.push_back(-1);
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
  cout<<"Tracker layers:"<<GetLastBarrelLayer()+1<<endl;
  cout<<"Missing outer tracker hits:"<<nMissingOuterTrackerHits<<endl;
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

int Track::GetLastBarrelLayer()
{
  int lastBarrelLayer = -1;
  
  for(int iHit=0;iHit<GetNdEdxHits();iHit++){
    if(detType[iHit] == 2) continue; // Endcaps, we don't care
    
    if(layer[iHit] > lastBarrelLayer){
      lastBarrelLayer = layer[iHit];
    }
  }
  
  return lastBarrelLayer-1;
}

double Track::GetAverageDedx()
{
  return GetTotalDedx()/nDedxClusters;
}
