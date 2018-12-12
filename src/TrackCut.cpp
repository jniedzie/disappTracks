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

TrackCut::~TrackCut()
{
  
}
