//
//  TrackCut.cpp
//
//  Created by Jeremi Niedziela on 16/07/2018.
//

#include "TrackCut.hpp"

#include <iostream>

using namespace std;

TrackCut::TrackCut(ECut cutType) :
nDedxClusters(range<int>()),
nDetIDs(range<int>()),
dedxPerCluster(range<double>()),
totalDeDx(range<double>()),
pt(range<double>()),
caloEmEnergy(range<double>()),
caloHadEnergy(range<double>()),
eta(range<double>()),
relativeIsolation(range<double>()),
sameNpixelHitsLayers(false),
sameNtrackerHitsLayers(false),
nPixelHits(range<int>()),
nPixelLayers(range<int>()),
nMissingInnerPixel(range<int>()),
nMissingMiddleTracker(range<int>()),
nMissingOuterTracker(range<int>())
{
  if(cutType&kEmpty) return;
  if(cutType&kPt50GeV)  pt = range<double>(50.0,  999999);
  if(cutType&kPt200GeV) pt = range<double>(200.0, 999999);
  if(cutType&kLowCalo){
      caloEmEnergy = range<double>(0.0, 0.5);
      caloEmEnergy = range<double>(0.0, 0.5);
  }
  if(cutType&kLowDEdx) totalDeDx = range<double>(0.0, 38.0);
  if(cutType&kShort)  nDedxClusters = range<int>(3, 4);
  if(cutType&kMedium) nDedxClusters = range<int>(3, 8);
  if(cutType&kEta2p4) eta = range<double>(-2.4, 2.4);
}

TrackCut::~TrackCut()
{
  
}
