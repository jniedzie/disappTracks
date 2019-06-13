//  Track.cpp
//
//  Created by Jeremi Niedziela on 16/07/2018.

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
eventMetMass(inf),
decayPoint(Point(0,0,0))
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

Track::Track(double _eta, double _phi, int _charge, int _nTrackerLayers, double _pt) :
eta(_eta),
phi(_phi),
charge(_charge),
nTrackerLayers(_nTrackerLayers),
pt(_pt),

mass(inf),
caloEmEnergy(-1.0),
caloHadEnergy(-1.0),
relativeIsolation(inf),
dxy(inf),
dxyErr(inf),
dz(inf),
dzErr(inf),
pid(inf),
mcMatch(inf),
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
eventMetMass(inf),
decayPoint(Point(0,0,0))
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
}

Track::Track(const Track &t) :
  dedx(t.dedx),
  subDetId(t.subDetId),
  sizeX(t.sizeX),
  sizeY(t.sizeY),
  detType(t.detType),
  layer(t.layer),
  ladder(t.ladder),
  pt(t.pt),
  eta(t.eta),
  phi(t.phi),
  mass(t.mass),
  caloEmEnergy(t.caloEmEnergy),
  caloHadEnergy(t.caloHadEnergy),
  relativeIsolation(t.relativeIsolation),
  dxy(t.dxy),
  dxyErr(t.dxyErr),
  dz(t.dz),
  dzErr(t.dzErr),
  charge(t.charge),
  pid(t.pid),
  mcMatch(t.mcMatch),
  nTrackerLayers(t.nTrackerLayers),
  nPixelLayers(t.nPixelLayers),
  nTrackerHits(t.nTrackerHits),
  nPixelHits(t.nPixelHits),
  nMissingInnerPixelHits(t.nMissingInnerPixelHits),
  nMissingOuterPixelHits(t.nMissingOuterPixelHits),
  nMissingInnerStripHits(t.nMissingInnerStripHits),
  nMissingOuterStripHits(t.nMissingOuterStripHits),
  nMissingInnerTrackerHits(t.nMissingInnerTrackerHits),
  nMissingOuterTrackerHits(t.nMissingOuterTrackerHits),
  nMissingMiddleTrackerHits(t.nMissingMiddleTrackerHits),
  nDetIDs(t.nDetIDs),
  nDedxClusters(t.nDedxClusters),
  eventMetPt(t.eventMetPt),
  eventMetEta(t.eventMetEta),
  eventMetPhi(t.eventMetPhi),
  eventMetMass(t.eventMetMass),
  decayPoint(t.decayPoint)
{

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
  cout<<"PID:"<<pid<<"\trel iso:"<<relativeIsolation<<endl;
  cout<<"eta:"<<eta<<"\tphi:"<<phi<<"\tpT:"<<pt<<endl;
  cout<<"Tracker layers:"<<nTrackerLayers<<endl;
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
  int lastBarrelLayer = -inf;
  int shift = 0; // shift layer index for strips
  
  for(int iHit=0;iHit<GetNdEdxHits();iHit++){
    if(detType[iHit] == 2) continue; // Endcaps, we don't care
    if(detType[iHit] == 0) shift = 4; // Strips, shift layer number by 4 pixel layers
    if(detType[iHit] == 1) shift = 0; // Pixel, don't shift
    
    if(layer[iHit]+shift > lastBarrelLayer){
      lastBarrelLayer = layer[iHit]+shift;
    }
  }
  
  return lastBarrelLayer-1;
}

double Track::GetAverageDedx()
{
  return GetTotalDedx()/nDedxClusters;
}
