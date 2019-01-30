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
