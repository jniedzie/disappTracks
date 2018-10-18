//
//  TrackCut.cpp
//
//  Created by Jeremi Niedziela on 16/07/2018.
//

#include "TrackCut.hpp"

#include <iostream>

using namespace std;

TrackCut::TrackCut(ECut cutType) :
minDedxClusters(0),
maxDedxClusters(999999),
minDets(0),
maxDets(999999),
minDedxPerCluster(0.0),
minTotalDeDx(0.0),
maxTotalDeDx(999999),
minPt(0.0),
maxEmCalo(999999),
maxHadCalo(999999),
maxEta(999999),
maxRelIso(999999),
sameNpixelHitsLayers(false),
sameNtrackerHitsLayers(false),
minNpixelHits(0),
maxNpixelHits(999999),
minNpixelLayers(0),
maxNpixelLayers(999999),
minMissingInnerPixel(0),
maxMissingInnerPixel(999999),
minMissingMiddleTracker(0),
maxMissingMiddleTracker(999999),
minMissingOuterTracker(0),
maxMissingOuterTracker(999999)
{
  if(cutType&kEmpty) return;
  if(cutType&kPt50GeV)  minPt = 50.0;
  if(cutType&kPt200GeV) minPt = 200.0;
  if(cutType&kLowCalo){
      maxEmCalo = 0.5;
      maxHadCalo = 0.5;
  }
  if(cutType&kLowDEdx) maxTotalDeDx = 38.0;
  if(cutType&kShort){
      minDedxClusters = 3;
      maxDedxClusters = 4;
  }
  if(cutType&kMedium){
      minDedxClusters = 3;
      maxDedxClusters = 8;
  }
  if(cutType&kEta2p4) maxEta = 2.4;
}

TrackCut::~TrackCut()
{
  
}
