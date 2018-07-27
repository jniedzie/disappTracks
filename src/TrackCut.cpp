//
//  TrackCut.cpp
//
//  Created by Jeremi Niedziela on 16/07/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#include "TrackCut.hpp"

#include <iostream>

using namespace std;

TrackCut::TrackCut(ECut cutType) :
minDedxClusters(0),
maxDedxClusters(999999),
minDedxPerCluster(0.0),
minTotalDeDx(0.0),
maxTotalDeDx(999999),
minPt(0.0),
maxEmCalo(999999),
maxHadCalo(999999)
{
  if(cutType&kEmpty) return;
  if(cutType&kHighPt) minPt = 200.0;
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
}

TrackCut::~TrackCut()
{
  
}
