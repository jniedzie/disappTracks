//
//  TrackCut.cpp
//
//  Created by Jeremi Niedziela on 16/07/2018.
//

#include "TrackCut.hpp"

#include <iostream>

using namespace std;

TrackCut::TrackCut() :
pt(range<double>()),
eta(range<double>()),
caloEmEnergy(range<double>()),
caloHadEnergy(range<double>()),
relativeIsolation(range<double>()),
sameNpixelHitsLayers(false),
sameNtrackerHitsLayers(false),
requireMcMatch(false),
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
trackMetDeltaPhi(range<double>())
{

}

TrackCut::TrackCut(const TrackCut &c)
{
  pt = c.pt;
  eta = c.eta;
  caloEmEnergy = c.caloEmEnergy;
  caloHadEnergy = c.caloHadEnergy;
  relativeIsolation = c.relativeIsolation;
  sameNpixelHitsLayers = c.sameNpixelHitsLayers;
  sameNtrackerHitsLayers = c.sameNtrackerHitsLayers;
  requireMcMatch = c.requireMcMatch;
  nLayers = c.nLayers;
  nPixelLayers = c.nPixelLayers;
  nPixelHits = c.nPixelHits;
  nMissingInnerPixel = c.nMissingInnerPixel;
  nMissingOuterTracker = c.nMissingOuterTracker;
  nMissingMiddleTracker = c.nMissingMiddleTracker;
  dedxPerCluster = c.dedxPerCluster;
  totalDeDx = c.totalDeDx;
  nDetIDs = c.nDetIDs;
  nDedxClusters = c.nDedxClusters;
  trackMetDeltaPhi = c.trackMetDeltaPhi;
}

TrackCut::~TrackCut()
{
  
}
