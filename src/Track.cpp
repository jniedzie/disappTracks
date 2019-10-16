//  Track.cpp
//
//  Created by Jeremi Niedziela on 16/07/2018.

#include "Track.hpp"
#include "Logger.hpp"

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

double Track::GetDedxInBarrelLayer(int iLayer)
{
  double totalDedx=0;

  for(int iHit=0; iHit<dedx.size(); iHit++){
    
    if(detType[iHit] == 1){ // pixel barrel
      if(layer[iHit] != iLayer) continue; // check that it's the correct layer
      totalDedx += dedx[iHit];
    }
    else if(detType[iHit] == 0){
      if(layer[iHit]+4 != iLayer) continue;
      totalDedx += dedx[iHit];
    }
  }
  return totalDedx;
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
    if(detType[iHit] == 2 || layer[iHit]==0) continue; // Endcaps, we don't care
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

double Track::GetMinDedx()
{
  return *min(dedx.begin(), dedx.end());
}

double Track::GetMaxDedx()
{
  return *max(dedx.begin(), dedx.end());
}

double Track::GetDedxLikelihood()
{
  double likelihoodProduct = 1.;
  double logLikelihoodSum = 0.;
 
  int nHitsProcessed = 0;
  
  for(int iHit=0; iHit<GetNdEdxHits(); iHit++){
    double hitDedx = GetDeDxForHit(iHit);
    if(hitDedx==0) continue;
    
    double likelihood = 1.0;
    
    if(detType[iHit] == 2)      likelihood = GetLandgaus(100, hitDedx);                     // Endcaps
    else if(detType[iHit] == 0) likelihood = GetLandgaus(GetLayerForHit(iHit)+4, hitDedx);  // Strips
    else if(detType[iHit] == 1) likelihood = GetLandgaus(GetLayerForHit(iHit)  , hitDedx);  // Pixels
    else Log(0)<<"Error -- unknown detector type: "<<detType[iHit]<<"\n";
    
    likelihoodProduct *= likelihood;
    logLikelihoodSum  -= log(likelihood);
    
    if(likelihood != 1.0) nHitsProcessed++;
  }
  
  return logLikelihoodSum/nHitsProcessed;
}


double Track::GetLandgaus(int iLayer, double hitDedx)
{
  if (hitDedx<0) return 1.0;
  
  map<int, vector<double>> par = {
    {1,   {0.24, 1.60, 0.1, 0.55}}, // BPIX 1
    {2,   {0.22, 2.45, 0.1, 0.34}}, // BPIX 2
    {3,   {0.22, 2.45, 0.1, 0.34}}, // BPIX 3
    {4,   {0.23, 2.70, 0.1, 0.32}}, // BPIX 4
    {100, {0.24, 2.50, 0.1, 0.30}}, // FPIX
  };
  
  if(par.find(iLayer) == par.end()){
    Log(0)<<"Error -- parameters for likelihood in layer "<<iLayer<<" not found.\n";
    return 1.0;
  }
  
  // Numeric constants
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  double mpshift  = -0.22278298;       // Landau maximum location
  
  // Control constants
  double np = 100.0;      // number of convolution steps
  double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  
  // Variables
  double xx, mpc, fland, xlow,xupp, step, i;
  double sum = 0.0;
  
  
  // MP shift correction
  mpc = par[iLayer][1] - mpshift * par[iLayer][0];
  
  // Range of convolution integral
  xlow = hitDedx - sc * par[iLayer][3];
  xupp = hitDedx + sc * par[iLayer][3];
  
  step = (xupp-xlow) / np;
  
  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[iLayer][0]) / par[iLayer][0];
    sum += fland * TMath::Gaus(hitDedx, xx, par[iLayer][3]);
    
    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[iLayer][0]) / par[iLayer][0];
    sum += fland * TMath::Gaus(hitDedx, xx, par[iLayer][3]);
  }
  
  return (par[iLayer][2] * step * sum * invsq2pi / par[iLayer][3]);
  
}
