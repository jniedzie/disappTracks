//  TrackCut.cpp
//
//  Created by Jeremi Niedziela on 16/07/2018.

#include "TrackCut.hpp"

TrackCut::TrackCut() :

pt(range<double>()),
eta(range<double>()),
vetoCracks(false),
caloEmEnergy(range<double>()),
caloHadEnergy(range<double>()),
relativeIsolation(range<double>()),
trackMetDeltaPhi(range<double>()),

requiresMcMatch(false),

nLayers(range<int>()),
nPixelLayers(range<int>()),
nPixelHits(range<int>()),
nMissingInnerPixel(range<int>()),
nMissingOuterTracker(range<int>()),
nMissingMiddleTracker(range<int>()),

dedxPerCluster(range<double>()),
totalDeDx(range<double>()),
nDetIDs(range<int>()),
nDedxClusters(range<int>()),

requiresCorrectNlayers(false),
requiresCorrectCharge(false)
{

}

TrackCut::TrackCut(const TrackCut &c)
{
  pt                     = c.pt;
  eta                    = c.eta;
  vetoCracks             = c.vetoCracks;
  caloEmEnergy           = c.caloEmEnergy;
  caloHadEnergy          = c.caloHadEnergy;
  relativeIsolation      = c.relativeIsolation;
  trackMetDeltaPhi       = c.trackMetDeltaPhi;
  
  requiresMcMatch        = c.requiresMcMatch;
  
  nLayers                = c.nLayers;
  nPixelLayers           = c.nPixelLayers;
  nPixelHits             = c.nPixelHits;
  nMissingInnerPixel     = c.nMissingInnerPixel;
  nMissingOuterTracker   = c.nMissingOuterTracker;
  nMissingMiddleTracker  = c.nMissingMiddleTracker;
  
  dedxPerCluster         = c.dedxPerCluster;
  totalDeDx              = c.totalDeDx;
  nDetIDs                = c.nDetIDs;
  nDedxClusters          = c.nDedxClusters;
  
  requiresCorrectNlayers = c.requiresCorrectNlayers;
  requiresCorrectCharge  = c.requiresCorrectCharge;
}

TrackCut::~TrackCut()
{
  
}
